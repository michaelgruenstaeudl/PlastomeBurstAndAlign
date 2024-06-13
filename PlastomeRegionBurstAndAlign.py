#!/usr/bin/env python3
""" Extracts and aligns coding and non-coding regions across multiple plastid genomes
"""
__version__ = "m_gruenstaeudl@fhsu.edu|Wed 22 Nov 2023 04:35:09 PM CST"

# ------------------------------------------------------------------------------#
# IMPORTS
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO, Nexus, SeqRecord, AlignIO
from Bio.Align import Applications  # line necessary; see: https://www.biostars.org/p/13099/
from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition
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
    def __init__(self, plastid_data: 'PlastidGenomeData', mainhelper: 'MainHelpers'):
        """Parses all GenBank flatfiles of a given folder and extracts
        all sequence annotations of the type specified by the user for each flatfile
        INPUT: user specification on cds/int/igs
        """
        self.plastid_data = plastid_data
        self.mainhelper = mainhelper

    def conduct_extraction(self):
        """Conduct extraction
        INPUT:  input folder, user specification on cds/int/igs
        OUTPUT: nucleotide and protein dictionaries
        """
        for f in self.plastid_data.files:
            self._extract_record(f)

        if not self.plastid_data.main_odict_nucl.items():
            log.critical(f"No items in main dictionary: {self.mainhelper.out_dir}")
            raise Exception()

    def _extract_record(self, file: str):
        log.info(f"  parsing {file}")
        filepath = os.path.join(self.mainhelper.in_dir, file)
        record = SeqIO.read(filepath, "genbank")

        if self.mainhelper.select_mode == "cds":
            self._extract_cds(record)
        elif self.mainhelper.select_mode == "igs":
            self._extract_igs(record)
        elif self.mainhelper.select_mode == "int":
            self._extract_int(record)

    def _extract_cds(self, rec: SeqRecord):
        """Extracts all CDS (coding sequences = genes) from a given sequence record
        OUTPUT: saves to global main_odict_nucl and to global main_odict_prot
        """

        def trim_mult_three(in_seq):
            trim_char = len(in_seq) % 3
            if trim_char > 0 and seq_obj[:3] == "ATG":
                in_seq = in_seq[:-trim_char]
            elif trim_char > 0:
                in_seq = None
            return in_seq

        def save_seq_to_dict(seq: SeqRecord, name: str, odict: OrderedDict, odict_key: str):
            record = SeqRecord.SeqRecord(
                seq, id=name, name="", description=""
            )
            if odict_key in odict.keys():
                odict[odict_key].append(record)
            else:
                odict[odict_key] = [record]

        features = [
            f for f in rec.features if f.type == "CDS" and "gene" in f.qualifiers
        ]
        for feature in features:
            gene_name = feature.qualifiers["gene"][0]
            seq_name = f"{gene_name}_{rec.name}"

            # Step 1. Extract nucleotide sequence of each gene
            seq_obj = feature.extract(rec).seq
            seq_obj = trim_mult_three(seq_obj)

            if seq_obj is None:
                log.warning(f"{seq_name} does not have a clear reading frame. Skipping this gene.")
                continue # Skip to next gene of the for loop.

            save_seq_to_dict(seq_obj, seq_name, self.plastid_data.main_odict_nucl, gene_name)

            # Step 2. Translate nucleotide sequence to amino acid sequence
            prot_obj = seq_obj.translate(
                table=11)  #, cds=True) # Getting error TTA is not stop codon.

            # Step 3. Save protein sequence to output dictionary
            save_seq_to_dict(prot_obj, seq_name, self.plastid_data.main_odict_prot, gene_name)

    def _extract_igs(self, rec: SeqRecord):
        """Extracts all IGS (intergenic spacers) from a given sequence record
        OUTPUT: saves to global main_odict_nucl
        """
        # Step 1. Extract all genes from record (i.e., cds, trna, rrna)
        # Resulting list contains adjacent features in order of appearance on genome
        # Note: No need to include "if feature.type=='tRNA'", because all tRNAs are also annotated as genes
        all_genes = [
            f for f in rec.features if (f.type == "gene" and "gene" in f.qualifiers)
        ]
        # Note: The statement "if feature.qualifier['gene'][0] not 'matK'" is necessary, as
        # matK is located inside trnK
        all_genes_minus_matk = [
            f for f in all_genes if f.qualifiers["gene"][0] != "matK"
        ]

        # Step 2. Loop through genes
        for count, idx in enumerate(range(0, len(all_genes_minus_matk) - 1), 1):
            cur_feat = all_genes_minus_matk[idx]
            cur_feat_name = cur_feat.qualifiers["gene"][0]
            adj_feat = all_genes_minus_matk[idx + 1]
            adj_feat_name = adj_feat.qualifiers["gene"][0]

            # Step 3. Define names of IGS
            if "gene" in cur_feat.qualifiers and "gene" in adj_feat.qualifiers:
                cur_feat_name_safe = sub(
                    r"\W", "", cur_feat_name.replace("-", "_")
                )
                adj_feat_name_safe = sub(
                    r"\W", "", adj_feat_name.replace("-", "_")
                )
                igs_name = cur_feat_name_safe + "_" + adj_feat_name_safe
                inv_igs_name = adj_feat_name_safe + "_" + cur_feat_name_safe

                # Only operate on genes that do not have compound locations (as it only messes things up)
                if (
                    type(cur_feat.location) is not CompoundLocation
                    and type(adj_feat.location) is not CompoundLocation
                ):
                    # Step 4. Make IGS SeqFeature
                    start_pos = ExactPosition(cur_feat.location.end)  # +1)
                    # Note: It's unclear if +1 is needed here.
                    end_pos = ExactPosition(adj_feat.location.start)
                    if int(start_pos) >= int(end_pos):
                        continue  # If there is no IGS, then simply skip this gene pair
                    else:
                        try:
                            exact_location = FeatureLocation(start_pos, end_pos)
                        except Exception as e:
                            log.warning(
                                f"\t{rec.name}: Exception occurred for IGS between "
                                f"`{cur_feat_name}` (start pos: {start_pos}) and "
                                f"`{adj_feat_name}` (end pos:{end_pos}). "
                                f"Skipping this IGS ...\n"
                                f"Error message: {e}"
                            )
                            continue
                    # Step 5. Make IGS SeqRecord
                    seq_obj = exact_location.extract(rec).seq
                    seq_name = igs_name + "_" + rec.name
                    seq_rec = SeqRecord.SeqRecord(
                        seq_obj, id=seq_name, name="", description=""
                    )
                    # Step 6. Attach seqrecord to growing dictionary
                    if (
                        igs_name in self.plastid_data.main_odict_nucl.keys()
                        or inv_igs_name in self.plastid_data.main_odict_nucl.keys()
                    ):
                        if igs_name in self.plastid_data.main_odict_nucl.keys():
                            tmp = self.plastid_data.main_odict_nucl[igs_name]
                            tmp.append(seq_rec)
                            self.plastid_data.main_odict_nucl[igs_name] = tmp
                        if inv_igs_name in self.plastid_data.main_odict_nucl.keys():
                            pass  # Don't count IGS in the IRs twice
                    else:
                        self.plastid_data.main_odict_nucl[igs_name] = [seq_rec]
                # Handle genes with compound locations
                else:
                    log.warning(
                        f"{rec.name}: the IGS between `{cur_feat_name}` and `{adj_feat_name}` is "
                        f"currently not handled and would have to be extracted manually. "
                        f"Skipping this IGS ..."
                    )
                    continue

    def _extract_int(self, rec: SeqRecord):
        """Extracts all INT (introns) from a given sequence record
        OUTPUT: saves to global main_odict_nucl
        """
        features = [
            f for f in rec.features if f.type == "CDS" or f.type == "tRNA"
        ]
        for feature in features:
            try:
                gene_name_base = feature.qualifiers["gene"][0]
                gene_name_base_safe = sub(
                    r"\W", "", gene_name_base.replace("-", "_")
                )
            except Exception as e:
                log.warning(
                    f"Unable to extract gene name for CDS starting "
                    f"at `{feature.location.start}` of `{rec.id}`. "
                    f"Skipping feature ...\n"
                    f"Error message: {e}"
                )
                continue
            ### Inner Function - Start ###
            def extract_internal_intron(rec, feature, gene_name, offset):
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
                    seq_name = gene_name + "_" + rec.name
                    seq_obj = feature.extract(rec).seq  # Here the actual extraction is conducted
                    seq_rec = SeqRecord.SeqRecord(
                        seq_obj, id=seq_name, name="", description=""
                    )
                    return seq_rec, gene_name
                except Exception as e:
                    log.critical(
                        f"Unable to conduct intron extraction for {feature.qualifiers['gene']}.\n"
                        f"Error message: {e}"
                    )
                    raise Exception()
            ### Inner Function - End ###
            # Step 1. Limiting the search to CDS containing introns
            # Step 1.a. If one intron in gene:
            if len(feature.location.parts) == 2:
                try:
                    gene_name = f"{gene_name_base_safe}_intron1"
                    seq_rec, gene_name = extract_internal_intron(
                        rec, feature, gene_name, 0
                    )
                    if gene_name not in self.plastid_data.main_odict_nucl.keys():
                        self.plastid_data.main_odict_nucl[gene_name] = [seq_rec]
                    else:
                        self.plastid_data.main_odict_nucl[gene_name].append(seq_rec)
                except Exception as e:
                    some_id = list(feature.qualifiers.values())[0]
                    log.warning(
                        f"An error for `{some_id}` occurred.\n"
                        f"Error message: {e}"
                    )
                    pass
            # Step 1.b. If two introns in gene:
            if len(feature.location.parts) == 3:
                copy_feature = deepcopy(
                    feature
                )  # Important b/c feature is overwritten in extract_internal_intron()
                try:
                    gene_name = f"{gene_name_base_safe}_intron1"
                    seq_rec, gene_name = extract_internal_intron(
                        rec, feature, gene_name, 0
                    )
                    if gene_name not in self.plastid_data.main_odict_nucl.keys():
                        self.plastid_data.main_odict_nucl[gene_name] = [seq_rec]
                    else:
                        self.plastid_data.main_odict_nucl[gene_name].append(seq_rec)
                except Exception as e:
                    some_id = list(feature.qualifiers.values())[0]
                    log.critical(
                        f"An error for `{some_id}` occurred.\n"
                        f"Error message: {e}"
                    )
                    raise Exception()
                    # pass
                feature = copy_feature
                try:
                    gene_name = f"{gene_name_base_safe}_intron2"
                    seq_rec, gene_name = extract_internal_intron(
                        rec, feature, gene_name, 1
                    )

                    if gene_name not in self.plastid_data.main_odict_intron2.keys():
                        self.plastid_data.main_odict_intron2[gene_name] = [seq_rec]
                    else:
                        self.plastid_data.main_odict_intron2[gene_name].append(seq_rec)
                except Exception as e:
                    some_id = list(feature.qualifiers.values())[0]
                    log.warning(
                        f"An issue occurred for gene `{some_id}`.\n"
                        f"Error message: {e}"
                    )
                    pass
        self.plastid_data.main_odict_nucl.update(self.plastid_data.main_odict_intron2)


# -----------------------------------------------------------------#
class BiopythonExceptions(Exception):
    pass
# -----------------------------------------------------------------#

class DataCleaning:
    def __init__(self, plastid_data: 'PlastidGenomeData', mainhelper: 'MainHelpers'):
        """Cleans the nucleotide and protein dictionaries
        INPUT:  nucleotide and protein dictionaries
        OUTPUT: nucleotide and protein dictionaries
        """
        self.plastid_data = plastid_data
        self.mainhelper = mainhelper
        log.info("cleaning extracted sequence annotations")

    def clean(self):
        self.remove_duplicate_annos()
        self.remove_annos_if_below_minnumtaxa()
        self.remove_annos_if_below_minseqlength()
        self.remove_orfs()
        self.remove_user_defined_genes()

    def remove_duplicate_annos(self):
        log.info("  removing duplicate annotations")
        ### Inner Function - Start ###
        def remove_duplicates(my_dict: dict):
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
        remove_duplicates(self.plastid_data.main_odict_nucl)
        if self.mainhelper.select_mode == "cds":
            remove_duplicates(self.plastid_data.main_odict_prot)

    def remove_annos_if_below_minseqlength(self):
        log.info(f"  removing annotations whose longest sequence is shorter than {self.mainhelper.min_seq_length} bp")
        for k, v in list(self.plastid_data.main_odict_nucl.items()):
            longest_seq = max([len(s.seq) for s in v])
            if longest_seq < self.mainhelper.min_seq_length:
                log.info(f"    removing {k} due to minimum sequence length setting")
                del self.plastid_data.main_odict_nucl[k]
                if self.plastid_data.main_odict_prot:
                    del self.plastid_data.main_odict_prot[k]

    def remove_annos_if_below_minnumtaxa(self):
        log.info(f"  removing annotations that occur in fewer than {self.mainhelper.min_num_taxa} taxa")
        for k, v in list(self.plastid_data.main_odict_nucl.items()):
            if len(v) < self.mainhelper.min_num_taxa:
                log.info(f"    removing {k} due to minimum number of taxa setting")
                del self.plastid_data.main_odict_nucl[k]
                if self.plastid_data.main_odict_prot:
                    del self.plastid_data.main_odict_prot[k]

    def remove_orfs(self):
        log.info("  removing ORFs")
        list_of_orfs = [orf for orf in self.plastid_data.main_odict_nucl.keys() if "orf" in orf]
        for orf in list_of_orfs:
            del self.plastid_data.main_odict_nucl[orf]
            if self.plastid_data.main_odict_prot:
                del self.plastid_data.main_odict_prot[orf]

    def remove_user_defined_genes(self):
        log.info("  removing user-defined genes")
        if self.mainhelper.exclude_list:
            for excluded in self.mainhelper.exclude_list:
                if excluded in self.plastid_data.main_odict_nucl:
                    del self.plastid_data.main_odict_nucl[excluded]
                    if self.mainhelper.select_mode == "cds" and self.plastid_data.main_odict_prot:
                        del self.plastid_data.main_odict_prot[excluded]
                else:
                    log.warning(f"    Region `{excluded}` to be excluded but not present in infile.")
                    pass

# -----------------------------------------------------------------#

class AlignmentCoordination:
    def __init__(self, plastid_data: 'PlastidGenomeData', mainhelper: 'MainHelpers'):
        """Coordinates the alignment of nucleotide or protein sequences
        INPUT:  foo bar baz
        OUTPUT: foo bar baz
        """
        self.plastid_data = plastid_data
        self.mainhelper = mainhelper
        self.success_list = None
        log.info("conducting the alignment of extracted sequences")

    def save_regions_as_unaligned_matrices(self):
        """Takes a dictionary of nucleotide sequences and saves all sequences of the same region
        into an unaligned nucleotide matrix
        INPUT: dictionary of sorted nucleotide sequences of all regions
        OUTPUT: unaligned nucleotide matrix for each region, saved to file
        """
        log.info("saving individual regions as unaligned nucleotide matrices")
        for k, v in self.plastid_data.main_odict_nucl.items():
            # Define input and output names
            out_fn_unalign_nucl = os.path.join(self.mainhelper.out_dir, f"nucl_{k}.unalign.fasta")
            with open(out_fn_unalign_nucl, "w") as hndl:
                SeqIO.write(v, hndl, "fasta")

    def align(self):
        if self.mainhelper.select_mode == "cds":
            self.conduct_protein_MSA_and_backtranslate()
        else:
            self.conduct_nucleotide_MSA()

    def conduct_nucleotide_MSA(self):
        """
        Iterates over all unaligned nucleotide matrices and aligns each via a third-party software tool
        INPUT:  - dictionary of sorted nucleotide sequences of all regions (used only for region names!)
                - unaligned nucleotide matrices (present as files in FASTA format)
        OUTPUT: aligned nucleotide matrices (present as files in FASTA format)
        """
        log.info("conducting MSA based on nucleotide sequence data")
        log.info(f"  using {self.mainhelper.num_threads} CPUs")

        ### Inner Function - Start ###
        def process_single_nucleotide_MSA(k, num_threads):
            # Define input and output names
            out_fn_unalign_nucl = os.path.join(self.mainhelper.out_dir, f"nucl_{k}.unalign.fasta")
            out_fn_aligned_nucl = os.path.join(self.mainhelper.out_dir, f"nucl_{k}.aligned.fasta")
            # Step 1. Align matrices via third-party alignment tool
            self._mafft_align(out_fn_unalign_nucl, out_fn_aligned_nucl)

        ### Inner Function - End ###
        # Step 2. Use ThreadPoolExecutor to parallelize alignment and back-translation
        if self.plastid_data.main_odict_nucl.items():
            with ThreadPoolExecutor(max_workers=self.mainhelper.num_threads) as executor:
                future_to_nucleotide = {
                    executor.submit(process_single_nucleotide_MSA, k, self.mainhelper.num_threads): k
                    for k in self.plastid_data.main_odict_nucl.keys()
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

    def conduct_protein_MSA_and_backtranslate(self):
        """Iterates over all unaligned PROTEIN matrices, aligns them as proteins via
        third-party software, and back-translates each alignment to NUCLEOTIDES
        INPUT:  dictionary of sorted PROTEIN sequences of all regions
        OUTPUT: aligned nucleotide matrices (present as files in NEXUS format)
        """
        log.info("Conducting MSA based on protein sequence data, followed by back-translation to nucleotides")
        log.info(f"  using {self.mainhelper.num_threads} CPUs")

        ### Inner Function - Start ###
        def process_single_protein_MSA(k, v, num_threads):
            # Define input and output names
            out_fn_unalign_prot = os.path.join(self.mainhelper.out_dir, f"prot_{k}.unalign.fasta")
            out_fn_aligned_prot = os.path.join(self.mainhelper.out_dir, f"prot_{k}.aligned.fasta")
            out_fn_unalign_nucl = os.path.join(self.mainhelper.out_dir, f"nucl_{k}.unalign.fasta")
            out_fn_aligned_nucl = os.path.join(self.mainhelper.out_dir, f"nucl_{k}.aligned.fasta")
            # Step 1. Write unaligned protein sequences to file
            with open(out_fn_unalign_prot, "w") as hndl:
                SeqIO.write(v, hndl, "fasta")
            # Step 2. Align matrices based on their PROTEIN sequences via third-party alignment tool
            self._mafft_align(out_fn_unalign_prot, out_fn_aligned_prot)
            # Step 3. Conduct actual back-translation from PROTEINS TO NUCLEOTIDES
            try:
                backtranslator = BackTranslation(
                    self.plastid_data
                )
                backtranslator.perform_back_translation(
                    "fasta", out_fn_aligned_prot,
                    out_fn_unalign_nucl, out_fn_aligned_nucl, 11
                )
            except Exception as e:
                log.warning(
                    f"Unable to conduct back-translation of `{k}`. "
                    f"Error message: {e}."
                )
        ### Inner Function - End ###
        # Step 2. Use ThreadPoolExecutor to parallelize alignment and back-translation
        with ThreadPoolExecutor(max_workers=self.mainhelper.num_threads) as executor:
            future_to_protein = {
                executor.submit(process_single_protein_MSA, k, v, self.mainhelper.num_threads): k
                for k, v in self.plastid_data.main_odict_prot.items()
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
            input=input_file, adjustdirection=True, thread=self.mainhelper.num_threads
        )
        stdout, stderr = mafft_cline()
        with open(output_file, "w") as hndl:
            hndl.write(stdout)

    def collect_successful_MSAs(self):
        """Converts alignments to NEXUS format; then collect all successfully generated alignments
        INPUT:  dictionary of region names
        OUTPUT: list of alignments
        """
        log.info("collecting all successful alignments")
        success_list = []
        for k in self.plastid_data.main_odict_nucl.keys():
            # Step 1. Define input and output names
            aligned_nucl_fasta = os.path.join(self.mainhelper.out_dir, f"nucl_{k}.aligned.fasta")
            aligned_nucl_nexus = os.path.join(self.mainhelper.out_dir, f"nucl_{k}.aligned.nexus")
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

    def concatenate_successful_MSAs(self):
        log.info("concatenate all successful alignments (in no particular order)")

        # Step 1. Define output names
        out_fn_nucl_concat_fasta = os.path.join(
            self.mainhelper.out_dir, "nucl_" + str(len(self.success_list)) + "concat.aligned.fasta"
        )
        out_fn_nucl_concat_nexus = os.path.join(
            self.mainhelper.out_dir, "nucl_" + str(len(self.success_list)) + "concat.aligned.nexus"
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
    def __init__(self, plastid_data: 'PlastidGenomeData'):
        """Back-translates protein sequences to nucleotide sequences
        INPUT:  foo bar baz
        OUTPUT: foo bar baz
        """
        self.main_odict_nucl = plastid_data.main_odict_nucl
        self.main_odict_prot = plastid_data.main_odict_prot

    def translate_and_evaluate(self, identifier, nuc, prot, table):
        """Returns nucleotide sequence if works (can remove trailing stop)"""
        if len(nuc) % 3:
            log.warning(
                f"Nucleotide sequence for {identifier} is length {len(nuc)} (not a multiple of three)"
            )

        p = str(prot).upper().replace("*", "X")
        t = str(nuc.translate(table)).upper().replace("*", "X")
        if len(t) == len(p) + 1:
            if str(nuc)[-3:].upper() in ambiguous_generic_by_id[table].stop_codons:
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
            if str(nuc[0:3]).upper() in ambiguous_generic_by_id[table].start_codons:
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

    def backtranslate_individual_sequence(self, aligned_protein_record, unaligned_nucleotide_record, gap, table=0):
        if not gap or len(gap) != 1:
            raise ValueError("Please supply a single gap character")

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
        gap_codon = "-" * 3
        ######


        # Per https://biopython.org/docs/1.81/api/Bio.Seq.html this is proper replacement for depreciated method ungap.
        #ungapped_protein = aligned_protein_record.seq.ungap(gap)
        ungapped_protein = aligned_protein_record.seq.replace(gap,"")

        ungapped_nucleotide = unaligned_nucleotide_record.seq
        if table:
            ungapped_nucleotide = self.translate_and_evaluate(
                aligned_protein_record.id, ungapped_nucleotide, ungapped_protein, table
            )
        elif len(ungapped_protein) * 3 != len(ungapped_nucleotide):
            log.warning(
                f"Inconsistent lengths for {aligned_protein_record.id}, "
                f"ungapped protein {len(ungapped_protein)}, "
                f"tripled {len(ungapped_protein) * 3} vs "
                f"ungapped nucleotide {len(ungapped_nucleotide)}"
            )

        seq = []
        nuc = str(ungapped_nucleotide)
        for amino_acid in aligned_protein_record.seq:
            if amino_acid == gap:
                seq.append(gap_codon)
            else:
                seq.append(nuc[:3])
                nuc = nuc[3:]
        assert (not nuc), (f"Nucleotide sequence for {unaligned_nucleotide_record.id} "
                           f"longer than protein {aligned_protein_record.id}")

        aligned_nuc = unaligned_nucleotide_record[:]  # copy for most annotation
        aligned_nuc.letter_annotation = {}  # clear this
        aligned_nuc.seq = Seq("".join(seq))  # , alpha)  # Modification on 09-Sep-2022 by M. Gruenstaeudl
        assert len(aligned_protein_record.seq) * 3 == len(aligned_nuc)
        return aligned_nuc

    def backtranslate_coordinator(self, protein_alignment, nucleotide_records, key_function=None, gap=None, table=0):
        """Thread nucleotide sequences onto a protein alignment."""
        if key_function is None:
            key_function = lambda x: x
        if gap is None:
            gap = "-"

        aligned = []
        for protein in protein_alignment:
            try:
                nucleotide = nucleotide_records[key_function(protein.id)]
            except KeyError:
                raise ValueError(
                    f"Could not find nucleotide sequence for protein {protein.id}"
                )
            aligned.append(self.backtranslate_individual_sequence(protein, nucleotide, gap, table))
        return MultipleSeqAlignment(aligned)

    def perform_back_translation(self, align_format, prot_align_file, nuc_fasta_file, nuc_align_file, table=0):
        """Perform back-translation of a protein alignment to nucleotides.
        Parameters:
        align_format: Format of the alignment file (e.g., 'fasta')
        prot_align_file: Path to the file containing the aligned protein sequences
        nuc_fasta_file: Path to the file containing the unaligned nucleotide sequences
        nuc_align_file: Path to the output file for the back-translated nucleotide sequences
        table: Genetic code table number (default is 0)
        """

        # Step 1. Load the protein alignment
        prot_align = AlignIO.read(prot_align_file, align_format)
        # Step 2. Index the unaligned nucleotide sequences
        nuc_dict = SeqIO.index(nuc_fasta_file, "fasta")
        # Step 3. Perform back-translation
        nuc_align = self.backtranslate_coordinator(prot_align, nuc_dict, gap="-", table=table)
        # Step 4. Write the back-translated nucleotide alignment to a file
        with open(nuc_align_file, "w") as output_handle:
            AlignIO.write(nuc_align, output_handle, align_format)

# -----------------------------------------------------------------#


class MainHelpers:
    # class methods
    @classmethod
    def setup_logger(cls, args: argparse.Namespace) -> logging.Logger:
        logger = logging.getLogger(__name__)
        log_format = "%(asctime)s [%(levelname)s] %(message)s"
        log_level = logging.DEBUG if args.verbose else logging.INFO
        coloredlogs.install(fmt=log_format, level=log_level, logger=logger)
        return logger

    @classmethod
    def test_if_alignsoftw_present(cls, softw: str = "mafft"):
        if find_executable(softw) is None:
            log.critical(f"Unable to find alignment software `{softw}`")
            raise Exception()

    # constructor
    def __init__(self, args: argparse.Namespace):
        self._set_select_mode(args)
        self._set_in_dir(args)
        self._set_out_dir(args)
        self._set_fileext(args)
        self._set_exclude_list(args)
        self._set_min_seq_length(args)
        self._set_min_num_taxa(args)
        self._set_num_threads(args)

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


class PlastidGenomeData:
    def __init__(self, mainhelper: MainHelpers):
        self._set_mode(mainhelper)
        self._set_main_odict_nucl()
        self._set_main_odict_prot()
        self._set_main_odict_intron2()
        self._set_files(mainhelper)

    def _set_mode(self, mainhelper: MainHelpers):
        self.mode = mainhelper.select_mode

    def _set_main_odict_nucl(self):
        self.main_odict_nucl = OrderedDict()

    def _set_main_odict_prot(self):
        self.main_odict_prot = OrderedDict() if self.mode == "cds" else None

    def _set_main_odict_intron2(self):
        self.main_odict_intron2 = OrderedDict() if self.mode == "int" else None

    def _set_files(self, mainhelper: MainHelpers):
        self.files = [
            f for f in os.listdir(mainhelper.in_dir) if f.endswith(mainhelper.fileext)
        ]


# ------------------------------------------------------------------------------#
# MAIN
# ------------------------------------------------------------------------------#
def main(args: argparse.Namespace):
    mainhelper = MainHelpers(args)
    plastid_data = PlastidGenomeData(mainhelper)

    extractor = ExtractAndCollect(plastid_data, mainhelper)
    extractor.conduct_extraction()

    cleaner = DataCleaning(plastid_data, mainhelper)
    cleaner.clean()

    aligncoord = AlignmentCoordination(plastid_data, mainhelper)
    aligncoord.save_regions_as_unaligned_matrices()
    aligncoord.align()
    aligncoord.collect_successful_MSAs()
    aligncoord.concatenate_successful_MSAs()

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
    par_args = parser.parse_args()
    log = MainHelpers.setup_logger(par_args)
    MainHelpers.test_if_alignsoftw_present()
    main(par_args)
# ------------------------------------------------------------------------------#
# EOF
# ------------------------------------------------------------------------------#
