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

from Bio import BiopythonWarning
import warnings

# ------------------------------------------------------------------------------#
# DEBUGGING HELP
# import ipdb
# ipdb.set_trace()
# -----------------------------------------------------------------#
# CLASSES AND FUNCTIONS


class ExtractAndCollect:
    def __init__(self, select_mode):
        """Parses all GenBank flatfiles of a given folder and extracts
        all sequence annotations of the type specified by the user for each flatfile
        INPUT: user specification on cds/int/igs
        """
        self.select_mode = select_mode

        # Setup warnings for non-multiple of 3.
        self.mult_err_flag = False

        def handle_biopython_warning(message, category, filename, lineno, file=None, line=None):
            if "Partial codon, len(sequence) not a multiple of three" in str(message):
                self.mult_err_flag = True
            else:
            # Let other warnings pass through
                warnings.showwarning(message, category, filename, lineno, file, line)

        # Register the warning handler function.
        warnings.showwarning = handle_biopython_warning

        log.info("parsing GenBank flatfiles and extracting their sequence annotations")

    def conduct_extraction(self, in_dir, fileext):
        """Conduct extraction
        INPUT:  input folder, user specification on cds/int/igs
        OUTPUT: nucleotide and protein dictionaries
        """
        self.main_odict_nucl = OrderedDict()
        self.main_odict_prot = OrderedDict() if self.select_mode == "cds" else None
        main_odict_intron2 = OrderedDict() if self.select_mode == "int" else None

        files = [f for f in os.listdir(in_dir) if f.endswith(fileext)]
        
        for f in files:
            log.info(f"  parsing {f}")
            filename = os.path.join(in_dir, f)

            rec, mult_err_flag = self.read_seq(filename)
                
            if self.select_mode == "cds":
                self._extract_cds(rec, mult_err_flag)
            if self.select_mode == "igs":
                self._extract_igs(rec)
            if self.select_mode == "int":
                self._extract_int(rec, main_odict_intron2)
                self.main_odict_nucl.update(main_odict_intron2)

            if not self.main_odict_nucl.items():
                log.critical(f"No items in main dictionary: {out_dir}")
                raise Exception()
        return (self.main_odict_nucl, self.main_odict_prot)

    def _extract_cds(self, rec, notMult3_warning_flag):
        """Extracts all CDS (coding sequences = genes) from a given sequence record
        OUTPUT: saves to global main_odict_nucl and to global main_odict_prot
        """

        def trim_mult_three(in_seq):
            trim_char = len(in_seq) % 3
            if trim_char > 0:
                in_seq = in_seq[:-trim_char]
            return in_seq
        
        for feature in rec.features:
            if feature.type == "CDS":
                if "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]
                    seq_name = f"{gene_name}_{rec.name}"

                    # Step 1. Extract nucleotide sequence of each gene
                    seq_obj = feature.extract(rec).seq

                    if notMult3_warning_flag or len(seq_obj)%3 != 0:
                        if seq_obj[:3] == "ATG":    # If start codon 'ATG', trim sequence render multiple of 3.
                            seq_obj = trim_mult_three(seq_obj)
                        
                        else:
                            log.warning(f"{seq_name} does not have a clear reading frame. Skipping this gene.")
                            continue # Skip to next gene of the fore loop.

                    seq_rec = SeqRecord.SeqRecord(
                        seq_obj, id=seq_name, name="", description=""
                    )
                    if gene_name in self.main_odict_nucl.keys():
                        tmp = self.main_odict_nucl[gene_name]
                        tmp.append(seq_rec)
                        self.main_odict_nucl[gene_name] = tmp
                    else:
                        self.main_odict_nucl[gene_name] = [seq_rec]

                    # Step 2. Translate nucleotide sequence to amino acid sequence
                    
                    if notMult3_warning_flag:
                        if seq_obj is None:
                            log.warning(f"{seq_name} does not have a clear reading frame. Skipping this gene.")
                            continue    # Skip to next gene of the for loop.
                    
                    seq_obj = feature.extract(rec).seq.translate(table=11) #, cds=True) # Getting error TTA is not stop codon. 


                    # Step 3. Save protein sequence to output dictionary
                    seq_rec = SeqRecord.SeqRecord(
                        seq_obj, id=seq_name, name="", description=""
                    )
                    if gene_name in self.main_odict_prot.keys():
                        tmp = self.main_odict_prot[gene_name]
                        tmp.append(seq_rec)
                        self.main_odict_prot[gene_name] = tmp
                    else:
                        self.main_odict_prot[gene_name] = [seq_rec]

    def _extract_igs(self, rec):
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
                        igs_name in self.main_odict_nucl.keys()
                        or inv_igs_name in self.main_odict_nucl.keys()
                    ):
                        if igs_name in self.main_odict_nucl.keys():
                            tmp = self.main_odict_nucl[igs_name]
                            tmp.append(seq_rec)
                            self.main_odict_nucl[igs_name] = tmp
                        if inv_igs_name in self.main_odict_nucl.keys():
                            pass  # Don't count IGS in the IRs twice
                    else:
                        self.main_odict_nucl[igs_name] = [seq_rec]
                # Handle genes with compound locations
                else:
                    log.warning(
                        f"{rec.name}: the IGS between `{cur_feat_name}` and `{adj_feat_name}` is "
                        f"currently not handled and would have to be extracted manually. "
                        f"Skipping this IGS ..."
                    )
                    continue

    def _extract_int(self, rec, main_odict_intron2):
        """Extracts all INT (introns) from a given sequence record
        OUTPUT: saves to global main_odict_nucl
        """
        for feature in rec.features:
            if feature.type == "CDS" or feature.type == "tRNA":
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
                        if gene_name not in self.main_odict_nucl.keys():
                            self.main_odict_nucl[gene_name] = [seq_rec]
                        else:
                            self.main_odict_nucl[gene_name].append(seq_rec)
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
                        if gene_name not in self.main_odict_nucl.keys():
                            self.main_odict_nucl[gene_name] = [seq_rec]
                        else:
                            self.main_odict_nucl[gene_name].append(seq_rec)
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

                        if gene_name not in main_odict_intron2.keys():
                            main_odict_intron2[gene_name] = [seq_rec]
                        else:
                            main_odict_intron2[gene_name].append(seq_rec)
                    except Exception as e:
                        some_id = list(feature.qualifiers.values())[0]
                        log.warning(
                            f"An issue occurred for gene `{some_id}`.\n"
                            f"Error message: {e}"
                        )
                        pass

       
    def read_seq(self, filename):
        rec = SeqIO.read(filename, "genbank")
        err_result = self.mult_err_flag
        self.mult_err_flag = False  # Reset the flag.

        return rec, err_result
     
      
# -----------------------------------------------------------------#
class BiopythonExceptions(Exception):
    pass
# -----------------------------------------------------------------#

class DataCleaning:
    def __init__(self, main_odict_nucl, main_odict_prot, select_mode):
        """Cleans the nucleotide and protein dictionaries
        INPUT:  nucleotide and protein dictionaries
        OUTPUT: nucleotide and protein dictionaries
        """
        self.main_odict_nucl = main_odict_nucl
        self.main_odict_prot = main_odict_prot
        self.select_mode = select_mode
        log.info("cleaning extracted sequence annotations")

    def remove_duplicate_annos(self):
        log.info("  removing duplicate annotations")
        ### Inner Function - Start ###
        def remove_duplicates(my_dict):
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
        remove_duplicates(self.main_odict_nucl)
        if self.select_mode == "cds":
            remove_duplicates(self.main_odict_prot)

    def remove_annos_if_below_minseqlength(self, min_seq_length):
        log.info(f"  removing annotations whose longest sequence is shorter than {min_seq_length} bp")
        for k, v in list(self.main_odict_nucl.items()):
            longest_seq = max([len(s.seq) for s in v])
            if longest_seq < min_seq_length:
                log.info(f"    removing {k} due to minimum sequence length setting")
                del self.main_odict_nucl[k]
                if self.main_odict_prot:
                    del self.main_odict_prot[k]

    def remove_annos_if_below_minnumtaxa(self, min_num_taxa):
        log.info(f"  removing annotations that occur in fewer than {min_num_taxa} taxa")
        for k, v in list(self.main_odict_nucl.items()):
            if len(v) < min_num_taxa:
                log.info(f"    removing {k} due to minimum number of taxa setting")
                del self.main_odict_nucl[k]
                if self.main_odict_prot:
                    del self.main_odict_prot[k]

    def remove_orfs(self):
        log.info("  removing ORFs")
        list_of_orfs = [orf for orf in self.main_odict_nucl.keys() if "orf" in orf]
        for orf in list_of_orfs:
            del self.main_odict_nucl[orf]
            if self.main_odict_prot:
                del self.main_odict_prot[orf]

    def remove_user_defined_genes(self, exclude_list):
        log.info("  removing user-defined genes")
        if exclude_list:
            if self.select_mode == "int":
                exclude_list = [i + "_intron1" for i in exclude_list] + \
                               [i + "_intron2" for i in exclude_list]
            for excluded in exclude_list:
                if excluded in self.main_odict_nucl:
                    del self.main_odict_nucl[excluded]
                    if self.select_mode == "cds" and self.main_odict_prot:
                        del self.main_odict_prot[excluded]
                else:
                    log.warning(f"    Region `{excluded}` to be excluded but not present in infile.")
                    pass
        return (self.main_odict_nucl, self.main_odict_prot)

# -----------------------------------------------------------------#

class AlignmentCoordination:
    def __init__(self, main_odict_nucl, main_odict_prot):
        """Coordinates the alignment of nucleotide or protein sequences
        INPUT:  foo bar baz
        OUTPUT: foo bar baz
        """
        self.main_odict_nucl = main_odict_nucl
        self.main_odict_prot = main_odict_prot
        log.info("conducting the alignment of extracted sequences")

    def save_regions_as_unaligned_matrices(self):
        """Takes a dictionary of nucleotide sequences and saves all sequences of the same region
        into an unaligned nucleotide matrix
        INPUT: dictionary of sorted nucleotide sequences of all regions
        OUTPUT: unaligned nucleotide matrix for each region, saved to file
        """
        log.info("saving individual regions as unaligned nucleotide matrices")
        for k, v in self.main_odict_nucl.items():
            # Define input and output names
            out_fn_unalign_nucl = os.path.join(out_dir, f"nucl_{k}.unalign.fasta")
            with open(out_fn_unalign_nucl, "w") as hndl:
                SeqIO.write(v, hndl, "fasta")

    def conduct_nucleotide_MSA(self):
        """
        Iterates over all unaligned nucleotide matrices and aligns each via a third-party software tool
        INPUT:  - dictionary of sorted nucleotide sequences of all regions (used only for region names!)
                - unaligned nucleotide matrices (present as files in FASTA format)
        OUTPUT: aligned nucleotide matrices (present as files in FASTA format)
        """
        log.info("conducting MSA based on nucleotide sequence data")
        # Step 1. Determine number of CPU core available
        try:
            num_threads = os.cpu_count()
        except NotImplementedError:
            num_threads = multiprocessing.cpu_count()
        log.info(f"  using {num_threads} CPUs")

        ### Inner Function - Start ###
        def process_single_nucleotide_MSA(k, num_threads):
            # Define input and output names
            out_fn_unalign_nucl = os.path.join(out_dir, f"nucl_{k}.unalign.fasta")
            out_fn_aligned_nucl = os.path.join(out_dir, f"nucl_{k}.aligned.fasta")
            # Step 1. Align matrices via third-party alignment tool
            self._mafft_align(out_fn_unalign_nucl, out_fn_aligned_nucl, num_threads)
        ### Inner Function - End ###
        # Step 2. Use ThreadPoolExecutor to parallelize alignment and back-translation
        if self.main_odict_nucl.items():
            with ThreadPoolExecutor(max_workers=num_threads) as executor:
                future_to_nucleotide = {
                    executor.submit(process_single_nucleotide_MSA, k, num_threads): k
                    for k in self.main_odict_nucl.keys()
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

        # Step 1. Determine number of CPU core available
        try:
            num_threads = os.cpu_count()
        except NotImplementedError:
            num_threads = multiprocessing.cpu_count()
        log.info(f"  using {num_threads} CPUs")

        ### Inner Function - Start ###
        def process_single_protein_MSA(k, v, num_threads):
            # Define input and output names
            out_fn_unalign_prot = os.path.join(out_dir, f"prot_{k}.unalign.fasta")
            out_fn_aligned_prot = os.path.join(out_dir, f"prot_{k}.aligned.fasta")
            out_fn_unalign_nucl = os.path.join(out_dir, f"nucl_{k}.unalign.fasta")
            out_fn_aligned_nucl = os.path.join(out_dir, f"nucl_{k}.aligned.fasta")
            # Step 1. Write unaligned protein sequences to file
            with open(out_fn_unalign_prot, "w") as hndl:
                SeqIO.write(v, hndl, "fasta")
            # Step 2. Align matrices based on their PROTEIN sequences via third-party alignment tool
            self._mafft_align(out_fn_unalign_prot, out_fn_aligned_prot, num_threads)
            # Step 3. Conduct actual back-translation from PROTEINS TO NUCLEOTIDES
            try:
                backtranslator = BackTranslation(
                    self.main_odict_nucl, self.main_odict_prot
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
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            future_to_protein = {
                executor.submit(process_single_protein_MSA, k, v, num_threads): k
                for k, v in self.main_odict_prot.items()
            }
            for future in as_completed(future_to_protein):
                k = future_to_protein[future]
                try:
                    future.result()  # If needed, you can handle results here
                except Exception as e:
                    log.error(f"{k} generated an exception: {e}")

    def _mafft_align(self, input_file, output_file, num_threads):
        """Perform sequence alignment using MAFFT"""
        mafft_cline = Applications.MafftCommandline(
            input=input_file, adjustdirection=True, thread=num_threads
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
        for k in self.main_odict_nucl.keys():
            # Step 1. Define input and output names
            aligned_nucl_fasta = os.path.join(out_dir, f"nucl_{k}.aligned.fasta")
            aligned_nucl_nexus = os.path.join(out_dir, f"nucl_{k}.aligned.nexus")
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
        return success_list

    def concatenate_successful_MSAs(self, success_list):
        log.info("concatenate all successful alignments (in no particular order)")

        # Step 1. Define output names
        out_fn_nucl_concat_fasta = os.path.join(
            out_dir, "nucl_" + str(len(success_list)) + "concat.aligned.fasta"
        )
        out_fn_nucl_concat_nexus = os.path.join(
            out_dir, "nucl_" + str(len(success_list)) + "concat.aligned.nexus"
        )
        # Step 2. Do concatenation
        try:
            alignm_concat = Nexus.Nexus.combine(
                success_list
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
    def __init__(self, main_odict_nucl, main_odict_prot):
        """Back-translates protein sequences to nucleotide sequences
        INPUT:  foo bar baz
        OUTPUT: foo bar baz
        """
        self.main_odict_nucl = main_odict_nucl
        self.main_odict_prot = main_odict_prot
        #log.info("conducting a backtranslation")

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


class MainHelperFunctions:
    def __init__(self):
        """Helper functions for main"""
        pass

    def setup_logger(self, verbose):
        global log
        log = logging.getLogger(__name__)
        log_format = "%(asctime)s [%(levelname)s] %(message)s"
        log_level = logging.DEBUG if verbose else logging.INFO
        coloredlogs.install(fmt=log_format, level=log_level, logger=log)

    def unpack_input_parameters(self, args):
        in_dir = args.inpd
        if not os.path.exists(in_dir):
            logging.critical(f"Input directory `{in_dir}` does not exist.")
            raise Exception()
        global out_dir
        out_dir = args.outd
        if not os.path.exists(out_dir):
            logging.critical(f"Output directory `{out_dir}` does not exist.")
            raise Exception()
        fileext = args.fileext
        exclude_list = args.excllist
        min_seq_length = args.minseqlength
        min_num_taxa = args.minnumtaxa
        select_mode = args.selectmode.lower()
        verbose = args.verbose
        return (
            in_dir,
            out_dir,
            fileext,
            exclude_list,
            min_seq_length,
            min_num_taxa,
            select_mode,
            verbose,
        )

    def test_if_alignsoftw_present(self, softw):
        if find_executable(softw) is not None:
            pass
        else:
            log.critical(f"Unable to find alignment software `{softw}`")
            raise Exception()

# ------------------------------------------------------------------------------#
# MAIN
# ------------------------------------------------------------------------------#
def main(args):
    mainhelper = MainHelperFunctions()
    (
        in_dir,
        out_dir,
        fileext,
        exclude_list,
        min_seq_length,
        min_num_taxa,
        select_mode,
        verbose,
    ) = mainhelper.unpack_input_parameters(args)
    mainhelper.setup_logger(verbose)
    mainhelper.test_if_alignsoftw_present("mafft")

    extractor = ExtractAndCollect(select_mode)
    main_odict_nucl, main_odict_prot = extractor.conduct_extraction(in_dir, fileext)

    cleaner = DataCleaning(main_odict_nucl, main_odict_prot, select_mode)
    cleaner.remove_duplicate_annos()
    cleaner.remove_annos_if_below_minnumtaxa(min_num_taxa)
    cleaner.remove_annos_if_below_minseqlength(min_seq_length)
    cleaner.remove_orfs()
    main_odict_nucl, main_odict_prot = cleaner.remove_user_defined_genes(exclude_list)

    aligncoord = AlignmentCoordination(main_odict_nucl, main_odict_prot)
    aligncoord.save_regions_as_unaligned_matrices()
    if not select_mode == "cds":
        aligncoord.conduct_nucleotide_MSA()
    if select_mode == "cds":
        aligncoord.conduct_protein_MSA_and_backtranslate()
    success_list = aligncoord.collect_successful_MSAs()
    aligncoord.concatenate_successful_MSAs(success_list)

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
        "--verbose",
        "-v",
        action="version",
        version="%(prog)s " + __version__,
        help="(Optional) Enable verbose logging",
        default=True,
    )
    args = parser.parse_args()
    main(args)
# ------------------------------------------------------------------------------#
# EOF
# ------------------------------------------------------------------------------#
