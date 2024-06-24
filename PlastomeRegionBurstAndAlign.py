#!/usr/bin/env python3
""" Extracts and aligns coding and non-coding regions across multiple plastid genomes
"""
__version__ = "m_gruenstaeudl@fhsu.edu|Wed 22 Nov 2023 04:35:09 PM CST"

# ------------------------------------------------------------------------------#
# IMPORTS
import argparse
import bisect
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Union, List
from Bio import SeqIO, Nexus, SeqRecord, AlignIO
from Bio.Align import Applications  # line necessary; see: https://www.biostars.org/p/13099/
from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition, SeqFeature
import coloredlogs
from collections import OrderedDict
from copy import deepcopy
from distutils.spawn import find_executable
from io import StringIO
import logging
import multiprocessing
import os
from re import sub
import sys
from Bio.Data.CodonTable import ambiguous_generic_by_id
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq

# ------------------------------------------------------------------------------#
# DEBUGGING HELP
# import ipdb
# ipdb.set_trace()
# -----------------------------------------------------------------#
# CLASSES AND FUNCTIONS


class ExtractAndCollect:
    def __init__(self, plastid_data: 'PlastidData', user_params: 'UserParameters'):
        """Parses all GenBank flatfiles of a given folder and extracts
        all sequence annotations of the type specified by the user for each flatfile
        INPUT: user specification on cds/int/igs
        """
        self.plastid_data = plastid_data
        self.user_params = user_params

    def extract(self):
        """Conduct extraction
        INPUT:  input folder, user specification on cds/int/igs
        OUTPUT: nucleotide and protein dictionaries
        """
        for f in self.plastid_data.files:
            self._extract_rec(f)

        if not self.plastid_data.nucleotides.items():
            log.critical(f"No items in main dictionary: {self.user_params.out_dir}")
            raise Exception()

    def _extract_rec(self, file: str):
        log.info(f"  parsing {file}")
        filepath = os.path.join(self.user_params.in_dir, file)
        record = SeqIO.read(filepath, "genbank")

        if self.user_params.select_mode == "cds":
            self._extract_cds(record)
        elif self.user_params.select_mode == "igs":
            self._extract_igs(record)
        elif self.user_params.select_mode == "int":
            self._extract_int(record)

    def _extract_cds(self, rec: SeqRecord):
        """Extracts all CDS (coding sequences = genes) from a given sequence record
        OUTPUT: saves to global main_odict_nucl and to global main_odict_prot
        """
        features = [
            f for f in rec.features if f.type == "CDS" and "gene" in f.qualifiers
        ]
        for feature in features:
            # Step 1. Extract nucleotide sequence of each gene and add to dictionary
            gene = GeneFeature(rec, feature)
            self.plastid_data.add_feature(gene)

            # Step 2. Translate nucleotide sequence to protein and add to dictionary
            protein = ProteinFeature(gene)
            self.plastid_data.add_protein(protein)

    def _extract_igs(self, rec: SeqRecord):
        """Extracts all IGS (intergenic spacers) from a given sequence record
        OUTPUT: saves to global main_odict_nucl
        """
        # Step 1. Extract all genes from record (i.e., cds, trna, rrna)
        # Resulting list contains adjacent features in order of appearance on genome
        # Note: No need to include "if feature.type=='tRNA'", because all tRNAs are also annotated as genes
        # Note: The statement "if feature.qualifier['gene'][0] not 'matK'" is necessary, as
        # matK is located inside trnK
        all_genes = [
            f for f in rec.features if f.type == "gene" and "gene" in f.qualifiers and f.qualifiers["gene"][0] != "matK"
        ]
        all_genes = self._split_compound(all_genes, rec)

        # Step 2. Loop through genes
        for count, idx in enumerate(range(0, len(all_genes) - 1), 1):
            cur_feat = all_genes[idx]
            adj_feat = all_genes[idx + 1]
            igs = IntergenicFeature(rec, cur_feat, adj_feat)

            # Only operate on genes that do not have compound locations (as it only messes things up)
            if not igs.compound_loc():
                # Step 3. Make IGS SeqFeature
                igs.set_seq_obj()

                # Step 4. Attach IGS to growing dictionary
                self.plastid_data.add_igs(igs)

            # Handle genes with compound locations
            else:
                log.warning(
                    f"{rec.name}: the IGS between `{igs.cur_name}` and `{igs.adj_name}` is "
                    f"currently not handled and would have to be extracted manually. "
                    f"Skipping this IGS ..."
                )
                continue

    def _extract_int(self, rec: SeqRecord):
        """Extracts all INT (introns) from a given sequence record
        OUTPUT: saves to global main_odict_nucl
        """
        # Step 1. Limiting the search to CDS containing introns
        features = [
            f for f in rec.features if (f.type == "CDS" or f.type == "tRNA") and "gene" in f.qualifiers
        ]
        for feature in features:
            # Step 1.a. If one intron in gene:
            if len(feature.location.parts) == 2:
                intron = IntronFeature(rec, feature)
                self.plastid_data.add_feature(intron)

            # Step 1.b. If two introns in gene:
            elif len(feature.location.parts) == 3:
                feature_copy = deepcopy(
                    feature
                )  # Important b/c feature is overwritten in extract_internal_intron()

                intron = IntronFeature(rec, feature)
                self.plastid_data.add_feature(intron)

                intron = IntronFeature(rec, feature_copy, 1)
                self.plastid_data.add_feature(intron)

    @staticmethod
    def _split_compound(genes: List[SeqFeature], record: SeqRecord):
        # find the compound features and remove them from the gene list
        compound_features = [
            f for f in genes if type(f.location) is CompoundLocation
        ]
        if len(compound_features) == 0:
            return genes
        genes = [
            f for f in genes if f not in compound_features
        ]

        log.info(f"  Resolving genes with compound locations for {record.name}")
        # create simple features from the compound features
        # sort the list by end location for proper handling of repeated annotations
        simple_features = []
        for f in compound_features:
            simple_features.extend(SeqFeature(p, f.type, f.id, f.qualifiers) for p in f.location.parts)
        simple_features.sort(key=lambda f: f.location.end, reverse=True)

        # find end locations of features for insertion index finding
        end_positions = [
            f.location.end for f in genes
        ]

        # insert the simple features at the correct indices in the gene list if applicable
        for simple_feature in simple_features:
            # extract feature location and find proper index
            insert_location = simple_feature.location
            insert_end = insert_location.end
            insertion_index = bisect.bisect_left(end_positions, insert_end)

            # only attempt insertion if there is no overlap with the previous gene
            # TODO: merge with overlapped previous gene if same gene?
            previous_feature = None if insertion_index == 0 else genes[insertion_index - 1]
            is_after_previous = not previous_feature or previous_feature.location.end < insert_location.start
            if not is_after_previous:
                continue

            # feature used for other insertion checks
            current_feature = None if insertion_index == len(genes) else genes[insertion_index]

            # directly insert feature if not overlapping with the current gene at this index
            is_before_current = not current_feature or insert_end < current_feature.location.start
            if is_before_current:
                genes.insert(insertion_index, simple_feature)
                end_positions.insert(insertion_index, insert_end)
                continue

            # if there is an overlap with the current feature which is this same gene
            # we assume it is a duplicate annotation;
            # in this case we decide to keep the longer of the two
            # TODO: merge with overlapped current gene if same gene instead?
            insert_gene = simple_feature.qualifiers["gene"][0]
            current_gene = current_feature.qualifiers["gene"][0]
            is_same_gene = current_gene == insert_gene
            if is_same_gene and len(simple_feature) > len(current_feature):
                genes[insertion_index] = simple_feature
                end_positions[insertion_index] = insert_end

        return genes


# -----------------------------------------------------------------#
class BiopythonExceptions(Exception):
    pass
# -----------------------------------------------------------------#


class DataCleaning:
    def __init__(self, plastid_data: 'PlastidData', user_params: 'UserParameters'):
        """Cleans the nucleotide and protein dictionaries
        INPUT:  nucleotide and protein dictionaries
        OUTPUT: nucleotide and protein dictionaries
        """
        self.plastid_data = plastid_data
        self.user_params = user_params
        log.info("cleaning extracted sequence annotations")

    def clean(self):
        self._dedup()
        self._remove_infreq()
        self._remove_short()
        self._remove_orfs()
        self._remove_excluded()

    def _dedup(self):
        log.info("  removing duplicate annotations")

        ### Inner Function - Start ###
        def remove_dups(my_dict: dict):
            """my_dict is modified in place"""
            for k, v in my_dict.items():
                unique_items = []
                seen_ids = set()
                for seqrec in v:
                    if seqrec.id not in seen_ids:
                        seen_ids.add(seqrec.id)
                        unique_items.append(seqrec)
                my_dict[k] = unique_items
        ### Inner Function - End ###

        remove_dups(self.plastid_data.nucleotides)
        if self.user_params.select_mode == "cds":
            remove_dups(self.plastid_data.proteins)

    def _remove_short(self):
        log.info(f"  removing annotations whose longest sequence is shorter than {self.user_params.min_seq_length} bp")
        for k, v in list(self.plastid_data.nucleotides.items()):
            longest_seq = max([len(s.seq) for s in v])
            if longest_seq < self.user_params.min_seq_length:
                log.info(f"    removing {k} due to minimum sequence length setting")
                del self.plastid_data.nucleotides[k]
                if self.plastid_data.proteins:
                    del self.plastid_data.proteins[k]

    def _remove_infreq(self):
        log.info(f"  removing annotations that occur in fewer than {self.user_params.min_num_taxa} taxa")
        for k, v in list(self.plastid_data.nucleotides.items()):
            if len(v) < self.user_params.min_num_taxa:
                log.info(f"    removing {k} due to minimum number of taxa setting")
                del self.plastid_data.nucleotides[k]
                if self.plastid_data.proteins:
                    del self.plastid_data.proteins[k]

    def _remove_orfs(self):
        log.info("  removing ORFs")
        list_of_orfs = [orf for orf in self.plastid_data.nucleotides.keys() if "orf" in orf]
        for orf in list_of_orfs:
            del self.plastid_data.nucleotides[orf]
            if self.plastid_data.proteins:
                del self.plastid_data.proteins[orf]

    def _remove_excluded(self):
        log.info("  removing user-defined genes")
        if self.user_params.exclude_list:
            for excluded in self.user_params.exclude_list:
                if excluded in self.plastid_data.nucleotides:
                    del self.plastid_data.nucleotides[excluded]
                    if self.user_params.select_mode == "cds" and self.plastid_data.proteins:
                        del self.plastid_data.proteins[excluded]
                else:
                    log.warning(f"    Region `{excluded}` to be excluded but not present in infile.")
                    pass

# -----------------------------------------------------------------#


class AlignmentCoordination:
    def __init__(self, plastid_data: 'PlastidData', user_params: 'UserParameters'):
        """Coordinates the alignment of nucleotide or protein sequences
        INPUT:  foo bar baz
        OUTPUT: foo bar baz
        """
        self.plastid_data = plastid_data
        self.user_params = user_params
        self.success_list = None
        log.info("conducting the alignment of extracted sequences")

    def save_unaligned(self):
        """Takes a dictionary of nucleotide sequences and saves all sequences of the same region
        into an unaligned nucleotide matrix
        INPUT: dictionary of sorted nucleotide sequences of all regions
        OUTPUT: unaligned nucleotide matrix for each region, saved to file
        """
        log.info("saving individual regions as unaligned nucleotide matrices")
        for k, v in self.plastid_data.nucleotides.items():
            # Define input and output names
            out_fn_unalign_nucl = os.path.join(self.user_params.out_dir, f"nucl_{k}.unalign.fasta")
            with open(out_fn_unalign_nucl, "w") as hndl:
                SeqIO.write(v, hndl, "fasta")

    def perform_MSA(self):
        if self.user_params.select_mode == "cds":
            self._prot_MSA()
        else:
            self._nuc_MSA()

    def _nuc_MSA(self):
        """
        Iterates over all unaligned nucleotide matrices and aligns each via a third-party software tool
        INPUT:  - dictionary of sorted nucleotide sequences of all regions (used only for region names!)
                - unaligned nucleotide matrices (present as files in FASTA format)
        OUTPUT: aligned nucleotide matrices (present as files in FASTA format)
        """
        log.info("conducting MSA based on nucleotide sequence data")
        log.info(f"  using {self.user_params.num_threads} CPUs")

        ### Inner Function - Start ###
        def single_nuc_MSA(k):
            # Define input and output names
            out_fn_unalign_nucl = os.path.join(self.user_params.out_dir, f"nucl_{k}.unalign.fasta")
            out_fn_aligned_nucl = os.path.join(self.user_params.out_dir, f"nucl_{k}.aligned.fasta")
            # Step 1. Align matrices via third-party alignment tool
            self._mafft_align(out_fn_unalign_nucl, out_fn_aligned_nucl)
        ### Inner Function - End ###

        # Step 2. Use ThreadPoolExecutor to parallelize alignment and back-translation
        if self.plastid_data.nucleotides.items():
            with ThreadPoolExecutor(max_workers=self.user_params.num_threads) as executor:
                future_to_nucleotide = {
                    executor.submit(single_nuc_MSA, k): k
                    for k in self.plastid_data.nucleotides.keys()
                }
                for future in as_completed(future_to_nucleotide):
                    k = future_to_nucleotide[future]
                    try:
                        future.result()  # If needed, you can handle results here
                    except Exception as e:
                        log.error("%r generated an exception: %s" % (k, e))
        else:
            log.critical("No items in nucleotide main dictionary to process")
            raise Exception()

    def _prot_MSA(self):
        """Iterates over all unaligned PROTEIN matrices, aligns them as proteins via
        third-party software, and back-translates each alignment to NUCLEOTIDES
        INPUT:  dictionary of sorted PROTEIN sequences of all regions
        OUTPUT: aligned nucleotide matrices (present as files in NEXUS format)
        """
        log.info("Conducting MSA based on protein sequence data, followed by back-translation to nucleotides")
        log.info(f"  using {self.user_params.num_threads} CPUs")

        ### Inner Function - Start ###
        def single_prot_MSA(k, v):
            # Define input and output names
            out_fn_unalign_prot = os.path.join(self.user_params.out_dir, f"prot_{k}.unalign.fasta")
            out_fn_aligned_prot = os.path.join(self.user_params.out_dir, f"prot_{k}.aligned.fasta")
            out_fn_unalign_nucl = os.path.join(self.user_params.out_dir, f"nucl_{k}.unalign.fasta")
            out_fn_aligned_nucl = os.path.join(self.user_params.out_dir, f"nucl_{k}.aligned.fasta")
            # Step 1. Write unaligned protein sequences to file
            with open(out_fn_unalign_prot, "w") as hndl:
                SeqIO.write(v, hndl, "fasta")
            # Step 2. Align matrices based on their PROTEIN sequences via third-party alignment tool
            self._mafft_align(out_fn_unalign_prot, out_fn_aligned_prot)
            # Step 3. Conduct actual back-translation from PROTEINS TO NUCLEOTIDES
            try:
                translator = BackTranslation(
                    "fasta", out_fn_aligned_prot,
                    out_fn_unalign_nucl, out_fn_aligned_nucl, 11
                )
                translator.backtranslate()
            except Exception as err:
                log.warning(
                    f"Unable to conduct back-translation of `{k}`. "
                    f"Error message: {err}."
                )
        ### Inner Function - End ###

        # Step 2. Use ThreadPoolExecutor to parallelize alignment and back-translation
        with ThreadPoolExecutor(max_workers=self.user_params.num_threads) as executor:
            future_to_protein = {
                executor.submit(single_prot_MSA, k, v): k
                for k, v in self.plastid_data.proteins.items()
            }
            for future in as_completed(future_to_protein):
                k = future_to_protein[future]
                try:
                    future.result()  # If needed, you can handle results here
                except Exception as e:
                    log.error(f"{k} generated an exception: {e}")

    def _mafft_align(self, input_file, output_file):
        """Perform sequence alignment using MAFFT"""
        mafft_cline = Applications.MafftCommandline(
            input=input_file, adjustdirection=True, thread=self.user_params.num_threads
        )
        stdout, stderr = mafft_cline()
        with open(output_file, "w") as hndl:
            hndl.write(stdout)

    def collect_MSAs(self):
        """Converts alignments to NEXUS format; then collect all successfully generated alignments
        INPUT:  dictionary of region names
        OUTPUT: list of alignments
        """
        log.info("collecting all successful alignments")
        success_list = []
        for k in self.plastid_data.nucleotides.keys():
            # Step 1. Define input and output names
            aligned_nucl_fasta = os.path.join(self.user_params.out_dir, f"nucl_{k}.aligned.fasta")
            aligned_nucl_nexus = os.path.join(self.user_params.out_dir, f"nucl_{k}.aligned.nexus")
            # Step 2. Convert FASTA alignment to NEXUS alignment
            try:
                AlignIO.convert(
                    aligned_nucl_fasta,
                    "fasta",
                    aligned_nucl_nexus,
                    "nexus",
                    molecule_type="DNA",
                )
            except Exception as e:
                log.warning(
                    f"Unable to convert alignment of `{k}` from FASTA to NEXUS.\n"
                    f"Error message: {e}"
                )
                continue  # skip to next k in loop, so that k is not included in success_list
            # Step 3. Import NEXUS files and append to list for concatenation
            try:
                alignm_nexus = AlignIO.read(aligned_nucl_nexus, "nexus")
                hndl = StringIO()
                AlignIO.write(alignm_nexus, hndl, "nexus")
                nexus_string = hndl.getvalue()
                # The following line replaces the gene name of sequence name with 'concat_'
                nexus_string = nexus_string.replace("\n" + k + "_", "\nconcat_")
                alignm_nexus = Nexus.Nexus.Nexus(nexus_string)
                success_list.append(
                    (k, alignm_nexus)
                )  # Function 'Nexus.Nexus.combine' needs a tuple.
            except Exception as e:
                log.warning(
                    f"Unable to add alignment of `{k}` to concatenation.\n"
                    f"Error message: {e}"
                )
                pass
        self.success_list = success_list

    def concat_MSAs(self):
        log.info("concatenate all successful alignments (in no particular order)")

        # Step 1. Define output names
        out_fn_nucl_concat_fasta = os.path.join(
            self.user_params.out_dir, "nucl_" + str(len(self.success_list)) + "concat.aligned.fasta"
        )
        out_fn_nucl_concat_nexus = os.path.join(
            self.user_params.out_dir, "nucl_" + str(len(self.success_list)) + "concat.aligned.nexus"
        )
        # Step 2. Do concatenation
        try:
            alignm_concat = Nexus.Nexus.combine(
                self.success_list
            )  # Function 'Nexus.Nexus.combine' needs a tuple
        except Exception as e:
            log.critical("Unable to concatenate alignments.\n" f"Error message: {e}")
            raise Exception()
        # Step 3. Write concatenated alignments to file in NEXUS format
        alignm_concat.write_nexus_data(filename=open(out_fn_nucl_concat_nexus, "w"))
        # Step 4. Convert the NEXUS file just generated to FASTA format
        AlignIO.convert(
            out_fn_nucl_concat_nexus, "nexus", out_fn_nucl_concat_fasta, "fasta"
        )

# -----------------------------------------------------------------#


class BackTranslation:
    def __init__(self, align_format: str,
                 prot_align_file: str, nuc_fasta_file: str, nuc_align_file: str,
                 table: int = 0, gap: str = "-"):
        """Back-translates protein sequences to nucleotide sequences
        Parameters:
        align_format: Format of the alignment file (e.g., 'fasta')
        prot_align_file: Path to the file containing the aligned protein sequences
        nuc_fasta_file: Path to the file containing the unaligned nucleotide sequences
        nuc_align_file: Path to the output file for the back-translated nucleotide sequences
        table: Genetic code table number (default is 0)
        gap: String that designates a nucleotide gap
        """
        if not gap or len(gap) != 1:
            raise ValueError("Please supply a single gap character")

        self.align_format = align_format
        self.prot_align_file = prot_align_file
        self.nuc_fasta_file = nuc_fasta_file
        self.nuc_align_file = nuc_align_file
        self.table = table
        self.gap = gap

    def _evaluate_nuc(self, identifier, nuc, prot):
        """Returns nucleotide sequence if works (can remove trailing stop)"""
        if len(nuc) % 3:
            log.warning(
                f"Nucleotide sequence for {identifier} is length {len(nuc)} (not a multiple of three)"
            )

        p = str(prot).upper().replace("*", "X")
        t = str(nuc.translate(self.table)).upper().replace("*", "X")
        if len(t) == len(p) + 1:
            if str(nuc)[-3:].upper() in ambiguous_generic_by_id[self.table].stop_codons:
                # Allow this...
                t = t[:-1]
                nuc = nuc[:-3]  # edit return value
        if len(t) != len(p):
            err = (
                f"Inconsistent lengths for {identifier}, ungapped protein {len(p)}, "
                f"tripled {len(p) * 3} vs ungapped nucleotide {len(nuc)}."
            )
            if t.endswith(p):
                err += f"\nThere are {len(t) - len(p)} extra nucleotides at the start."
            elif t.startswith(p):
                err += f"\nThere are {len(t) - len(p)} extra nucleotides at the end."
            elif p in t:
                err += "\nHowever, protein sequence found within translated nucleotides."
            elif p[1:] in t:
                err += "\nHowever, ignoring first amino acid, protein sequence found within translated nucleotides."
            log.warning(err)

        if t == p:
            return nuc
        elif p.startswith("M") and "M" + t[1:] == p:
            if str(nuc[0:3]).upper() in ambiguous_generic_by_id[self.table].start_codons:
                return nuc
            else:
                log.warning(
                    f"Translation check failed for {identifier}\n"
                    f"Would match if {nuc[0:3].upper()} was a start codon (check correct table used)"
                )

        else:
            m = "".join("." if x == y else "!" for (x, y) in zip(p, t))
            if len(prot) < 70:
                sys.stderr.write(f"Protein:     {p}\n")
                sys.stderr.write(f"             {m}\n")
                sys.stderr.write(f"Translation: {t}\n")
            else:
                for offset in range(0, len(p), 60):
                    sys.stderr.write(f"Protein:     {p[offset:offset + 60]}\n")
                    sys.stderr.write(f"             {m[offset:offset + 60]}\n")
                    sys.stderr.write(f"Translation: {t[offset:offset + 60]}\n\n")
            log.warning(f"Translation check failed for {identifier}\n")

    def _backtrans_seq(self, aligned_protein_record, unaligned_nucleotide_record):
        ######
        # Modification on 09-Sep-2022 by M. Gruenstaeudl
        # alpha = unaligned_nucleotide_record.seq.alphabet
        # if hasattr(alpha, "gap_char"):
        #    gap_codon = alpha.gap_char * 3
        #    assert len(gap_codon) == 3
        # else:
        #    from Bio.Alphabet import Gapped
        #    alpha = Gapped(alpha, gap)
        #    gap_codon = gap * 3
        ######


        # Per https://biopython.org/docs/1.81/api/Bio.Seq.html this is proper replacement for depreciated method ungap.
        #ungapped_protein = aligned_protein_record.seq.ungap(gap)
        ungapped_protein = aligned_protein_record.seq.replace(self.gap, "")

        ungapped_nucleotide = unaligned_nucleotide_record.seq
        if self.table:
            ungapped_nucleotide = self._evaluate_nuc(
                aligned_protein_record.id, ungapped_nucleotide, ungapped_protein
            )
        elif len(ungapped_protein) * 3 != len(ungapped_nucleotide):
            log.warning(
                f"Inconsistent lengths for {aligned_protein_record.id}, "
                f"ungapped protein {len(ungapped_protein)}, "
                f"tripled {len(ungapped_protein) * 3} vs "
                f"ungapped nucleotide {len(ungapped_nucleotide)}"
            )
        if ungapped_nucleotide is None:
            return None

        seq = []
        nuc = str(ungapped_nucleotide)
        gap_codon = self.gap * 3
        for amino_acid in aligned_protein_record.seq:
            if amino_acid == self.gap:
                seq.append(gap_codon)
            else:
                seq.append(nuc[:3])
                nuc = nuc[3:]
        if len(nuc) > 0:
            log.warning(
                f"Nucleotide sequence for {unaligned_nucleotide_record.id} "
                f"longer than protein {aligned_protein_record.id}"
            )
            return None

        aligned_nuc = unaligned_nucleotide_record[:]  # copy for most annotation
        aligned_nuc.letter_annotation = {}  # clear this
        aligned_nuc.seq = Seq("".join(seq))  # , alpha)  # Modification on 09-Sep-2022 by M. Gruenstaeudl
        if len(aligned_protein_record.seq) * 3 != len(aligned_nuc):
            return None

        return aligned_nuc

    def _backtrans_seqs(self, protein_alignment, nucleotide_records, key_function=lambda x: x):
        """Thread nucleotide sequences onto a protein alignment."""
        aligned = []
        for protein in protein_alignment:
            try:
                nucleotide = nucleotide_records[key_function(protein.id)]
            except KeyError:
                raise ValueError(
                    f"Could not find nucleotide sequence for protein {protein.id}"
                )
            sequence = self._backtrans_seq(protein, nucleotide)
            if sequence is not None:
                aligned.append(sequence)
        return MultipleSeqAlignment(aligned)

    def backtranslate(self):
        """Perform back-translation of a protein alignment to nucleotides."""

        # Step 1. Load the protein alignment
        prot_align = AlignIO.read(self.prot_align_file, self.align_format)
        # Step 2. Index the unaligned nucleotide sequences
        nuc_dict = SeqIO.index(self.nuc_fasta_file, "fasta")
        # Step 3. Perform back-translation
        nuc_align = self._backtrans_seqs(prot_align, nuc_dict)
        # Step 4. Write the back-translated nucleotide alignment to a file
        with open(self.nuc_align_file, "w") as output_handle:
            AlignIO.write(nuc_align, output_handle, self.align_format)

# -----------------------------------------------------------------#


class UserParameters:
    # constructor
    def __init__(self, pars: argparse.ArgumentParser):
        args = pars.parse_args()
        self._set_select_mode(args)
        self._set_in_dir(args)
        self._set_out_dir(args)
        self._set_fileext(args)
        self._set_exclude_list(args)
        self._set_min_seq_length(args)
        self._set_min_num_taxa(args)
        self._set_num_threads(args)
        self._set_verbose(args)

    # mutators
    def _set_select_mode(self, args: argparse.Namespace):
        self.select_mode = args.selectmode.lower()

    def _set_in_dir(self, args: argparse.Namespace):
        in_dir = args.inpd
        if not os.path.exists(in_dir):
            logging.critical(f"Input directory `{in_dir}` does not exist.")
            raise Exception()
        self.in_dir: str = in_dir

    def _set_out_dir(self, args: argparse.Namespace):
        out_dir = args.outd
        if not os.path.exists(out_dir):
            logging.critical(f"Output directory `{out_dir}` does not exist.")
            raise Exception()
        self.out_dir: str = out_dir

    def _set_fileext(self, args: argparse.Namespace):
        self.fileext = args.fileext

    def _set_exclude_list(self, args: argparse.Namespace):
        exclude_list = args.excllist
        if self.select_mode == "int":
            self.exclude_list = [i + "_intron1" for i in exclude_list] + \
                                [i + "_intron2" for i in exclude_list]
        else:
            self.exclude_list = exclude_list

    def _set_min_seq_length(self, args: argparse.Namespace):
        self.min_seq_length = args.minseqlength

    def _set_min_num_taxa(self, args: argparse.Namespace):
        self.min_num_taxa = args.minnumtaxa

    def _set_num_threads(self, args: argparse.Namespace):
        num_threads = args.numthreads
        if num_threads == "auto":
            try:
                num_threads = os.cpu_count()
            except NotImplementedError:
                num_threads = multiprocessing.cpu_count()
        else:
            num_threads = int(num_threads)
        self.num_threads = num_threads
        
    def _set_verbose(self, args: argparse.Namespace):
        self.verbose = args.verbose

# -----------------------------------------------------------------#


class PlastidData:
    @staticmethod
    def save_seq_to_dict(feature: Union['GeneFeature', 'ProteinFeature', 'IntronFeature'], odict: OrderedDict):
        record = SeqRecord.SeqRecord(
            feature.seq_obj, id=feature.seq_name, name="", description=""
        )
        if feature.gene_name in odict.keys():
            odict[feature.gene_name].append(record)
        else:
            odict[feature.gene_name] = [record]

    def __init__(self, user_params: UserParameters):
        self._set_mode(user_params)
        self._set_nucleotides()
        self._set_proteins()
        self._set_files(user_params)

    def _set_mode(self, user_params: UserParameters):
        self.mode = user_params.select_mode

    def _set_nucleotides(self):
        self.nucleotides = OrderedDict()

    def _set_proteins(self):
        self.proteins = OrderedDict() if self.mode == "cds" else None

    def _set_files(self, user_params: UserParameters):
        self.files = [
            f for f in os.listdir(user_params.in_dir) if f.endswith(user_params.fileext)
        ]

    def add_feature(self, feature: Union['GeneFeature', 'IntronFeature']):
        if feature.seq_obj is None:
            log.warning(f"{feature.seq_name} does not have a clear reading frame. Skipping this gene.")
        else:
            self.save_seq_to_dict(feature, self.nucleotides)

    def add_protein(self, protein: 'ProteinFeature'):
        if protein.seq_obj is not None:
            self.save_seq_to_dict(protein, self.proteins)

    def add_igs(self, igs: 'IntergenicFeature'):
        if igs.seq_obj is None:
            return

        igs.set_igs_names()
        record = SeqRecord.SeqRecord(
            igs.seq_obj, id=igs.seq_name, name="", description=""
        )

        if igs.igs_name in self.nucleotides.keys():
            self.nucleotides[igs.igs_name].append(record)
        elif igs.inv_igs_name in self.nucleotides.keys():
            pass  # Don't count IGS in the IRs twice
        else:
            self.nucleotides[igs.igs_name] = [record]

# -----------------------------------------------------------------#


class GeneFeature:
    def __init__(self, record: SeqRecord, feature: SeqFeature):
        self._set_gene_name(feature)
        self.seq_name = f"{self.gene_name}_{record.name}"
        self._set_seq_obj(record, feature)

    def _set_gene_name(self, feature: SeqFeature):
        self.gene_name = sub(
            r"\W", "", feature.qualifiers["gene"][0].replace("-", "_")
        )

    def _set_seq_obj(self, record: SeqRecord, feature: SeqFeature):
        self.seq_obj = feature.extract(record).seq
        self._trim_mult_three()

    def _trim_mult_three(self):
        trim_char = len(self.seq_obj) % 3
        if trim_char > 0 and self.seq_obj[:3] == "ATG":
            self.seq_obj = self.seq_obj[:-trim_char]
        elif trim_char > 0:
            self.seq_obj = None


class ProteinFeature:
    def __init__(self, gene: 'GeneFeature'):
        self.gene_name = gene.gene_name
        self.seq_name = gene.seq_name
        self._set_prot_obj(gene.seq_obj)

    def _set_prot_obj(self, seq_obj: SeqRecord):
        if seq_obj is None:
            self.seq_obj = None
        else:
            self.seq_obj = seq_obj.translate(table=11)  # Getting error TTA is not stop codon.


class IntronFeature:
    def __init__(self, record: SeqRecord, feature: SeqFeature, offset: int = 0):
        self._set_gene_name(feature, offset)
        self.seq_name = f"{self.gene_name}_{record.name}"
        self._set_seq_obj(record, feature, offset)

    def _set_gene_name(self, feature: SeqFeature, offset: int):
        self.gene_name = sub(
            r"\W", "", feature.qualifiers["gene"][0].replace("-", "_")
        ) + "_intron" + str(offset + 1)

    def _set_seq_obj(self, record, feature, offset):
        try:
            feature.location = FeatureLocation(
                feature.location.parts[offset].end,
                feature.location.parts[offset + 1].start
            )
        except Exception:
            feature.location = FeatureLocation(
                feature.location.parts[offset + 1].start,
                feature.location.parts[offset].end
            )
        try:
            self.seq_obj = feature.extract(record).seq
        except Exception as e:
            log.critical(
                f"Unable to conduct intron extraction for {feature.qualifiers['gene']}.\n"
                f"Error message: {e}"
            )
            raise Exception()


class IntergenicFeature:
    def __init__(self, record: SeqRecord, cur_feat: SeqFeature, adj_feat: SeqFeature):
        self.record = record
        self.cur_feat = cur_feat
        self.adj_feat = adj_feat
        self._set_gene_names()

        self.seq_obj = None
        self.igs_name = None
        self.inv_igs_name = None
        self.seq_name = None

    def _set_gene_names(self):
        self.cur_name = sub(
            r"\W", "", self.cur_feat.qualifiers["gene"][0].replace("-", "_")
        )
        self.adj_name = sub(
            r"\W", "", self.adj_feat.qualifiers["gene"][0].replace("-", "_")
        )

    def compound_loc(self):
        return type(self.cur_feat.location) is CompoundLocation or type(self.adj_feat.location) is CompoundLocation

    def set_seq_obj(self):
        # Note: It's unclear if +1 is needed here.
        start_pos = ExactPosition(self.cur_feat.location.end)  # +1)
        end_pos = ExactPosition(self.adj_feat.location.start)

        if int(start_pos) < int(end_pos):
            try:
                exact_location = FeatureLocation(start_pos, end_pos)
                self.seq_obj = exact_location.extract(self.record).seq
            except Exception as e:
                log.info(f"Start: {start_pos}, End: {end_pos}")
                log.warning(
                    f"\t{self.record.name}: Exception occurred for IGS between "
                    f"`{self.cur_name}` (start pos: {start_pos}) and "
                    f"`{self.adj_name}` (end pos:{end_pos}). "
                    f"Skipping this IGS ...\n"
                    f"Error message: {e}"
                )

    def set_igs_names(self):
        self.igs_name = f"{self.cur_name}_{self.adj_name}"
        self.inv_igs_name = f"{self.adj_name}_{self.cur_name}"
        self.seq_name = self.igs_name + "_" + self.record.name


# ------------------------------------------------------------------------------#
# MAIN
# ------------------------------------------------------------------------------#
def setup_logger(user_params: UserParameters) -> logging.Logger:
    logger = logging.getLogger(__name__)
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    log_level = logging.DEBUG if user_params.verbose else logging.INFO
    coloredlogs.install(fmt=log_format, level=log_level, logger=logger)
    return logger


def check_dependency(software: str = "mafft"):
    if find_executable(software) is None:
        log.critical(f"Unable to find alignment software `{software}`")
        raise Exception()


def main(user_params: UserParameters):
    plastid_data = PlastidData(user_params)

    extractor = ExtractAndCollect(plastid_data, user_params)
    extractor.extract()

    cleaner = DataCleaning(plastid_data, user_params)
    cleaner.clean()

    aligncoord = AlignmentCoordination(plastid_data, user_params)
    aligncoord.save_unaligned()
    aligncoord.perform_MSA()
    aligncoord.collect_MSAs()
    aligncoord.concat_MSAs()

    log.info("end of script\n")
    quit()


# ------------------------------------------------------------------------------#
# ARGPARSE
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Author|Version: " + __version__)
    # Required
    parser.add_argument(
        "--inpd",
        "-i",
        type=str,
        required=True,
        help="path to input directory (which contains the GenBank files)",
        default="./input",
    )
    # Optional
    parser.add_argument(
        "--outd",
        "-o",
        type=str,
        required=False,
        help="(Optional) Path to output directory",
        default="./output",
    )
    parser.add_argument(
        "--selectmode",
        "-s",
        type=str,
        required=False,
        help="(Optional) Type of regions to be extracted (i.e. `cds`, `int`, or `igs`)",
        default="cds",
    )
    parser.add_argument(
        "--fileext",
        "-f",
        type=str,
        required=False,
        help="(Optional) File extension of input files",
        default=".gb",
    )
    parser.add_argument(
        "--excllist",
        "-e",
        type=list,
        required=False,
        default=["rps12"],
        help="(Optional) List of genes to be excluded",
    )
    parser.add_argument(
        "--minseqlength",
        "-l",
        type=int,
        required=False,
        help="(Optional) Minimal sequence length (in bp) below which regions will not be extracted",
        default=3,
    )
    parser.add_argument(
        "--minnumtaxa",
        "-t",
        type=int,
        required=False,
        help="(Optional) Minimum number of taxa in which a region must be present to be extracted",
        default=2,
    )
    parser.add_argument(
        "--numthreads",
        "-n",
        type=str,
        required=False,
        help="(Optional) Number of CPUs to use; can be any positive integer or 'auto' (default)",
        default="auto",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="version",
        version="%(prog)s " + __version__,
        help="(Optional) Enable verbose logging",
        default=True,
    )
    params = UserParameters(parser)
    log = setup_logger(params)
    check_dependency()
    main(params)
# ------------------------------------------------------------------------------#
# EOF
# ------------------------------------------------------------------------------#
