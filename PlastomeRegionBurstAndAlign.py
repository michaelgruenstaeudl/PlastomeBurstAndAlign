#!/usr/bin/env python3
""" A Python tool to extract and align genes, introns, and intergenic spacers across hundreds of plastid genomes using associative arrays 
"""
__version__ = "m_gruenstaeudl@fhsu.edu|Thu Jul 18 12:44:16 PM CEST 2024"

# ------------------------------------------------------------------------------#
# IMPORTS
import argparse
import bisect
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from functools import partial
from time import sleep
from typing import Union, List, Callable, Tuple, Mapping, Optional
from Bio import SeqIO, Nexus, SeqRecord, AlignIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition, SeqFeature
import coloredlogs
from collections import OrderedDict, defaultdict
from copy import deepcopy
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
        log.info(f"parsing GenBank flatfiles and extracting their sequence annotations using {self.user_params.num_threads} CPUs")

        # Step 0. Extract first genome in list for feature ordering
        if self.user_params.order == "seq":
            first_genome = self.plastid_data.files.pop(0)
            nuc_dict, prot_dict = self._extract_recs([first_genome])
            self.plastid_data.add_nucleotides(nuc_dict)
            self.plastid_data.add_proteins(prot_dict)

        # Step 1. Create the data for each worker
        file_lists = split_list(self.plastid_data.files, self.user_params.num_threads * 2)

        # Step 2. Use ProcessPoolExecutor to parallelize extraction
        mp_context = multiprocessing.get_context("fork")  # same method on all platforms
        with ProcessPoolExecutor(max_workers=self.user_params.num_threads, mp_context=mp_context) as executor:
            future_to_nuc = [
                executor.submit(self._extract_recs, file_list) for file_list in file_lists
            ]
            for future in as_completed(future_to_nuc):
                nuc_dict, prot_dict = future.result()
                self.plastid_data.add_nucleotides(nuc_dict)
                self.plastid_data.add_proteins(prot_dict)

        # Step 3. Stop execution if no nucleotides were extracted
        if not self.plastid_data.nucleotides.items():
            log.critical(f"No items in main dictionary: {self.user_params.out_dir}")
            raise Exception()

    def _extract_recs(self, files: List[str]) -> Tuple['PlastidDict', 'PlastidDict']:
        nuc_dict = IntergenicDict() if self.user_params.select_mode == 'igs' else PlastidDict()
        prot_dict = PlastidDict()
        extract_rec = self._extract_rec_gen(nuc_dict, prot_dict)

        for f in files:
            extract_rec(f)
        return nuc_dict, prot_dict

    def _extract_rec_gen(self, nuc_dict: 'PlastidDict', prot_dict: 'PlastidDict') -> Callable[[str], None]:
        """
        Creates the partial function to be used in extract_rec and then
        returns extract_rec
        :param nuc_dict: a PlastidDict that nucleotides will be saved to
        :param prot_dict: a PlastidDict that proteins will be saved to for CDS
        :return: the extraction function that will be used on all files
        """
        if self.user_params.select_mode == "cds":
            extract_fun = partial(self._extract_cds, gene_dict=nuc_dict, protein_dict=prot_dict)
        elif self.user_params.select_mode == "igs":
            extract_fun = partial(self._extract_igs, igs_dict=nuc_dict)
        elif self.user_params.select_mode == "int":
            extract_fun = partial(self._extract_int, int_dict=nuc_dict)
        else:
            extract_fun = None

        def extract_rec(file: str):
            try:
                log.info(f"  parsing {file}")
                filepath = os.path.join(self.user_params.in_dir, file)
                record = SeqIO.read(filepath, "genbank")
                extract_fun(record)
            except Exception as e:
                log.error("%r generated an exception: %s" % (file, e))

        return extract_rec

    def _extract_cds(self, rec: SeqRecord, gene_dict: 'PlastidDict', protein_dict: 'PlastidDict'):
        """Extracts all CDS (coding sequences = genes) from a given sequence record
        OUTPUT: saves to global main_odict_nucl and to global main_odict_prot
        """
        features = self._cds_features(rec)
        for feature in features:
            # Step 1. Extract nucleotide sequence of each gene and add to dictionary
            gene = GeneFeature(rec, feature)
            gene_dict.add_feature(gene)

            # Step 2. Translate nucleotide sequence to protein and add to dictionary
            protein = ProteinFeature(gene=gene)
            protein_dict.add_feature(protein)

    def _extract_igs(self, rec: SeqRecord, igs_dict: 'PlastidDict'):
        """Extracts all IGS (intergenic spacers) from a given sequence record
        OUTPUT: saves to global main_odict_nucl
        """
        # Step 1a. Extract all genes from record (i.e., cds, trna, rrna)
        # Resulting list contains adjacent features in order of appearance on genome
        all_genes = self._igs_features(rec)

        # Step 1b. Split all compound location features into simple location features
        splitter = CompoundSplitting(all_genes, rec)
        splitter.split()

        # Step 2. Loop through genes
        for count, idx in enumerate(range(0, len(all_genes) - 1), 1):
            current_feat = all_genes[idx]
            subsequent_feat = all_genes[idx + 1]

            # Step 3. Make IGS SeqFeature
            igs = IntergenicFeature(rec, current_feat, subsequent_feat)

            # Step 4. Attach IGS to growing dictionary
            igs_dict.add_feature(igs)

    def _extract_int(self, rec: SeqRecord, int_dict: 'PlastidDict'):
        """Extracts all INT (introns) from a given sequence record
        OUTPUT: saves to global main_odict_nucl
        """
        features = self._int_features(rec)
        for feature in features:
            # Step 1.a. If one intron in gene:
            if len(feature.location.parts) == 2:
                intron = IntronFeature(rec, feature)
                int_dict.add_feature(intron)

            # Step 1.b. If two introns in gene:
            elif len(feature.location.parts) == 3:
                feature_copy = deepcopy(
                    feature
                )  # Important b/c feature is overwritten in extract_internal_intron()

                intron = IntronFeature(rec, feature)
                int_dict.add_feature(intron)

                intron = IntronFeature(rec, feature_copy, 1)
                int_dict.add_feature(intron)

    def _cds_features(self, record: SeqRecord) -> List[SeqFeature]:
        return [
            f for f in record.features if f.type == "CDS" and self._not_exclude(f)
        ]

    def _igs_features(self, record: SeqRecord) -> List[SeqFeature]:
        # Note: No need to include "if feature.type=='tRNA'", because all tRNAs are also annotated as genes
        return [
            f for f in record.features if f.type == "gene" and self._not_exclude(f)
        ]

    def _int_features(self, record: SeqRecord) -> List[SeqFeature]:
        # Limiting the search to CDS containing introns
        return [
            f for f in record.features if (f.type == "CDS" or f.type == "tRNA") and self._not_exclude(f)
        ]

    def _not_exclude(self, feature: SeqFeature) -> bool:
        gene = PlastidFeature.get_gene(feature)
        return gene and gene not in self.user_params.exclude_list and "orf" not in gene


# -----------------------------------------------------------------#


class CompoundSplitting:
    def __init__(self, genes: List[SeqFeature], record: SeqRecord):
        self.genes = genes
        self.record = record

    def split(self):
        log.info(f"  resolving genes with compound locations in {self.record.name}")
        log.info(f"   resolving cis-spliced genes in {self.record.name}")
        merger = ExonSpliceMerger(self.genes, self.record)
        merger.merge()
        # reorder genes, even if we don't end up inserting anything
        self.genes.sort(key=lambda gene: gene.location.end)
        log.info(f"   resolving trans-spliced genes in {self.record.name}")
        insertor = ExonSpliceInsertor(self.genes, self.record, merger.trans_list)
        insertor.insert()


class ExonSpliceMerger:
    def __init__(self, genes: List[SeqFeature], record: SeqRecord):
        self.genes = genes
        self.rec_name: str = record.name

    def merge(self):
        self._setup()
        self._resolve_cis()

    def _setup(self):
        self._index = len(self.genes)
        self.trans_list: List[SeqFeature] = []
        self._current = None

    def _resolve_cis(self):
        for gene in reversed(self.genes):
            self._update_genes(gene)

            # we will not handle trans exons
            if self._is_trans():
                self._remove_trans()
                continue

            if type(gene.location) is CompoundLocation:
                self._merge_cis_exons()
                self._print_align()

            if self._is_same_gene():
                self._merge_adj_exons()
                self._print_align()

    def _is_trans(self):
        return self._current.qualifiers.get("trans_splicing") or self._current_gene == "rps12"

    def _is_same_gene(self):
        return self._current_gene == self._subsequent_gene

    def _remove_trans(self):
        self.trans_list.append(self._current)
        del self.genes[self._index]

    def _update_genes(self, current: SeqFeature):
        self._index -= 1
        self._set_subsequent()
        self._set_current(current)

    def _set_current(self, current: SeqFeature):
        self._current = current
        self._current_gene = PlastidFeature.get_safe_gene(current)

    def _set_subsequent(self, subsequent_index: Optional[int] = None):
        # if there is no subsequent feature
        if subsequent_index == len(self.genes) or not self._current:
            self._subsequent_gene: str = ""
            self._subsequent_loc: Union[str, FeatureLocation] = ""
        # typical behavior when proceeding to next iteration in `_resolve_cis`
        elif subsequent_index is None:
            self._subsequent_gene = self._current_gene
            self._subsequent_loc = self._current.location
        # typical behavior when performing non-complex exon merge in `_merge_adj_exons`
        else:
            subsequent = self.genes[subsequent_index]
            self._subsequent_gene = PlastidFeature.get_safe_gene(subsequent)
            self._subsequent_loc = subsequent.location

    def _merge_cis_exons(self):
        loc_parts = self._current.location.parts
        gene_start = min(p.start for p in loc_parts)
        gene_end = max(p.end for p in loc_parts)
        self._current.location = FeatureLocation(gene_start, gene_end)

    def _merge_adj_exons(self):
        gene_start = min(self._current.location.start, self._subsequent_loc.start)
        gene_end = max(self._current.location.end, self._subsequent_loc.end)
        self._current.location = FeatureLocation(gene_start, gene_end)

        # delete merged exon, and update new subsequent feature
        del self.genes[self._index + 1]
        self._set_subsequent(self._index + 1)

    def _print_align(self):
        log.debug(
            f"   Merging exons of {self._current_gene} within {self.rec_name}\n"
            "-----------------------------------------------------------\n"
            f"\t{self._current_gene}\t\t\t{self._subsequent_gene}\n"
            f"\t{self._current.location}\t\t{self._subsequent_loc}\n"
            "-----------------------------------------------------------\n"
        )


class ExonSpliceInsertor:
    def __init__(self, genes: List[SeqFeature], record: SeqRecord,
                 compound_features: Optional[List[SeqFeature]] = None):
        self.genes = genes
        self.rec_name: str = record.name
        self.compound_features = compound_features if compound_features is not None else None

    def insert(self):
        if self.compound_features is None:
            self._find_compounds()
        if len(self.compound_features) == 0:
            return
        self._create_simple()
        self._insert_simple()

    def _find_compounds(self):
        # find the compound features and remove them from the gene list
        self.compound_features = [
            f for f in self.genes if type(f.location) is CompoundLocation
        ]
        # directly remove elements from list;
        # list comprehension would create a new gene list
        for feature in self.compound_features:
            self.genes.remove(feature)
        # reorder genes, even if we don't end up inserting anything
        self.genes.sort(key=lambda gene: gene.location.end)

    def _create_simple(self):
        self.simple_features = []
        for f in self.compound_features:
            self.simple_features.extend(
                SeqFeature(location=p, type=f.type, id=f.id, qualifiers=f.qualifiers) for p in f.location.parts
            )
        self.simple_features.sort(key=lambda feat: feat.location.end, reverse=True)

    def _insert_simple(self):
        # find end locations of features for insertion index finding
        self._end_positions = [
            f.location.end for f in self.genes
        ]

        # insert the simple features at the correct indices in the gene list if applicable
        for insert_feature in self.simple_features:
            self._set_insert(insert_feature)
            self._set_adj()
            self._set_adj_tests()

            # using adjacency checks, attempt to insert
            self._try_repositioning()
            self._try_merging()

    def _set_insert(self, insert: SeqFeature):
        self._insert = insert
        self._is_repositioned = False
        self._insert_gene = PlastidFeature.get_safe_gene(self._insert)

        # extract feature location and find proper index
        insert_location = self._insert.location
        self._insert_start = insert_location.start
        self._insert_end = insert_location.end
        self._insert_index = bisect.bisect_left(self._end_positions, self._insert_end)

    def _set_adj(self):
        # set appropriate adjacent features
        self._previous = None if self._insert_index == 0 else self.genes[self._insert_index - 1]
        self._current = None if self._insert_index == len(self.genes) else self.genes[self._insert_index]

        self._previous_gene = "\t" if not self._previous else PlastidFeature.get_safe_gene(self._previous)
        self._current_gene = "" if not self._current else PlastidFeature.get_safe_gene(self._current)

        self._previous_loc = "\t\t\t\t\t" if not self._previous else self._previous.location
        self._current_loc = "" if not self._current else self._current.location

    def _set_adj_tests(self):
        # checks for how to handle the insert feature
        self._is_same_previous = False if not self._previous else self._previous_gene == self._insert_gene
        self._is_same_current = False if not self._current else self._current_gene == self._insert_gene

        self._is_after_previous = not self._previous or self._previous_loc.end < self._insert_start
        self._is_before_current = not self._current or self._insert_end < self._current_loc.start

    def _try_repositioning(self):
        # if insert feature does not overlap with adjacent features, and is a different gene from the others,
        # directly insert
        if self._is_after_previous and self._is_before_current and not self._is_same_previous and not self._is_same_current:
            self._message = f"Repositioning {self._insert_gene} within {self.rec_name}"
            self._insert_at_index()

    def _try_merging(self):
        if self._is_repositioned:
            return

        self._is_merged = False
        self._merge_right()
        self._merge_left()

        # perform merge if needed
        if self._is_merged:
            self._insert = SeqFeature(
                location=FeatureLocation(self._insert_start, self._insert_end, self._insert.location.strand),
                type=self._insert.type, id=self._insert.id, qualifiers=self._insert.qualifiers
            )
            # new adjacent features
            self._set_adj()
            self._message = f"Merging exons of {self._insert_gene} within {self.rec_name}"
            self._insert_at_index()

    def _merge_right(self):
        # if insert and current feature are the same gene, and insert starts before current,
        # remove current feature, and update ending location
        if self._is_same_current and self._insert_start < self._current_loc.start:
            self._insert_end = self._current_loc.end
            self._remove_at_index()
            self._is_merged = True

    def _merge_left(self):
        # if insert and previous feature are the same gene,
        # use the smaller start location, and remove previous feature
        if self._is_same_previous:
            self._insert_start = min(self._insert_start, self._previous_loc.start)
            # elements in list will shift to left, so update index
            self._insert_index -= 1
            self._remove_at_index()
            self._is_merged = True

    def _insert_at_index(self):
        self.genes.insert(self._insert_index, self._insert)
        self._end_positions.insert(self._insert_index, self._insert_end)
        self._print_align()
        self._is_repositioned = True

    def _remove_at_index(self):
        del self.genes[self._insert_index]
        del self._end_positions[self._insert_index]

    def _print_align(self):
        log.debug(
            f"   {self._message}\n"
            "-----------------------------------------------------------\n"
            f"\t\t{self._previous_gene}\t\t\t\t{self._insert_gene}\t\t\t\t{self._current_gene}\n"
            f"\t{self._previous_loc}\t{self._insert.location}\t{self._current_loc}\n"
            "-----------------------------------------------------------\n"
        )


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

    def clean(self):
        log.info("cleaning extracted sequence annotations")
        log.info(f"  removing annotations that occur in fewer than {self.user_params.min_num_taxa} taxa")
        log.info(f"  removing annotations whose longest sequence is shorter than {self.user_params.min_seq_length} bp")
        for k, v in list(self.plastid_data.nucleotides.items()):
            self._remove_infreq(k, v)
            self._remove_short(k, v)

    def _remove_short(self, k, v):
        longest_seq = max([len(s.seq) for s in v])
        if longest_seq < self.user_params.min_seq_length:
            log.info(f"    removing {k} for not reaching the minimum sequence length defined")
            self.plastid_data.remove_nuc(k)

    def _remove_infreq(self, k, v):
        if len(v) < self.user_params.min_num_taxa:
            log.info(f"    removing {k} for not reaching the minimum number of taxa defined")
            self.plastid_data.remove_nuc(k)


# -----------------------------------------------------------------#


class AlignmentCoordination:
    def __init__(self, plastid_data: 'PlastidData', user_params: 'UserParameters'):
        """Coordinates the alignment of nucleotide or protein sequences
        INPUT:  foo bar baz
        OUTPUT: foo bar baz
        """
        self.plastid_data = plastid_data
        self.user_params = user_params
        self.success_list = []

    def save_unaligned(self):
        """Takes a dictionary of nucleotide sequences and saves all sequences of the same region
        into an unaligned nucleotide matrix
        INPUT: dictionary of sorted nucleotide sequences of all regions
        OUTPUT: unaligned nucleotide matrix for each region, saved to file
        """
        log.info("saving individual regions as unaligned nucleotide matrices")
        for nuc in self.plastid_data.nucleotides.keys():
            out_fn_unalign_nucl = os.path.join(self.user_params.out_dir, f"nucl_{nuc}.unalign.fasta")
            with open(out_fn_unalign_nucl, "w") as hndl:
                SeqIO.write(self.plastid_data.nucleotides.get(nuc), hndl, "fasta")
            # Write unaligned protein sequences to file
            if self.user_params.select_mode == "cds":
                out_fn_unalign_prot = os.path.join(self.user_params.out_dir, f"prot_{nuc}.unalign.fasta")
                with open(out_fn_unalign_prot, "w") as hndl:
                    SeqIO.write(self.plastid_data.proteins.get(nuc), hndl, "fasta")

    def perform_MSA(self):
        log.info("conducting the alignment of extracted sequences")
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
        log.info(f"conducting multiple sequence alignments based on nucleotide sequence data using {self.user_params.num_threads} CPUs")

        ### Inner Function - Start ###
        def single_nuc_MSA(k: str):
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
        log.info(f"conducting multiple sequence alignments based on protein sequence data, followed by back-translation to nucleotides using {self.user_params.num_threads} CPUs")

        ### Inner Function - Start ###
        def single_prot_MSA(k: str):
            # Define input and output names
            out_fn_unalign_prot = os.path.join(self.user_params.out_dir, f"prot_{k}.unalign.fasta")
            out_fn_aligned_prot = os.path.join(self.user_params.out_dir, f"prot_{k}.aligned.fasta")
            out_fn_unalign_nucl = os.path.join(self.user_params.out_dir, f"nucl_{k}.unalign.fasta")
            out_fn_aligned_nucl = os.path.join(self.user_params.out_dir, f"nucl_{k}.aligned.fasta")
            # Step 1. Align matrices based on their PROTEIN sequences via third-party alignment tool
            self._mafft_align(out_fn_unalign_prot, out_fn_aligned_prot)
            # Step 2. Conduct actual back-translation from PROTEINS TO NUCLEOTIDES
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
                executor.submit(single_prot_MSA, k): k
                for k in self.plastid_data.nucleotides.keys()
            }
            for future in as_completed(future_to_protein):
                k = future_to_protein[future]
                try:
                    future.result()  # If needed, you can handle results here
                except Exception as e:
                    log.error(f"{k} generated an exception: {e}")

    def _mafft_align(self, input_file: str, output_file: str):
        """Perform sequence alignment using MAFFT"""
        mafft_cmd = ["mafft", "--thread", str(self.user_params.num_threads), "--adjustdirection", input_file]
        with open(output_file, 'w') as hndl, open(os.devnull, 'w') as devnull:
            process = subprocess.Popen(mafft_cmd, stdout=hndl, stderr=devnull, text=True)
            process.wait()

    def collect_MSAs(self):
        """Converts alignments to NEXUS format; then collect all successfully generated alignments
        INPUT:  dictionary of region names
        OUTPUT: list of alignments
        """
        log.info("collecting all successful alignments")

        nuc_lists = split_list(list(self.plastid_data.nucleotides.keys()), self.user_params.num_threads * 2)
        mp_context = multiprocessing.get_context("fork")  # same method on all platforms
        with ProcessPoolExecutor(max_workers=self.user_params.num_threads, mp_context=mp_context) as executor:
            future_to_success = [
                executor.submit(self._collect_MSA_list, msa_list)
                for msa_list in nuc_lists
            ]
            for future in as_completed(future_to_success):
                success_list = future.result()
                if len(success_list) > 0:
                    self.success_list.extend(success_list)

    def _collect_MSA_list(self, nuc_list: List[str]) -> List[Tuple[str, MultipleSeqAlignment]]:
        def collect_MSA(k: str) -> Optional[Tuple[str, MultipleSeqAlignment]]:
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
                return None
            # Step 3. Import NEXUS files and append to list for concatenation
            try:
                alignm_nexus = AlignIO.read(aligned_nucl_nexus, "nexus")
                hndl = StringIO()
                AlignIO.write(alignm_nexus, hndl, "nexus")
                nexus_string = hndl.getvalue()
                # The following line replaces the gene name of sequence name with 'concat_'
                nexus_string = nexus_string.replace("\n" + k + "_", "\nconcat_")
                alignm_nexus = Nexus.Nexus.Nexus(nexus_string)
                return k, alignm_nexus  # Function 'Nexus.Nexus.combine' needs a tuple.
            except Exception as e:
                log.warning(
                    f"Unable to add alignment of `{k}` to concatenation.\n"
                    f"Error message: {e}"
                )
                return None

        success_list = []
        for nuc in nuc_list:
            alignm_tup = collect_MSA(nuc)
            if alignm_tup is not None:
                success_list.append(alignm_tup)
        return success_list

    def concat_MSAs(self):
        def concat_sync():
            # Write concatenated alignments to file in NEXUS format
            mp_context = multiprocessing.get_context("spawn")
            nexus_write = mp_context.Process(target=alignm_concat.write_nexus_data,
                                             kwargs={"filename": out_fn_nucl_concat_nexus})
            nexus_write.start()

            # Write concatenated alignments to file in FASTA format
            fasta_write = mp_context.Process(target=alignm_concat.export_fasta,
                                             kwargs={"filename": out_fn_nucl_concat_fasta})
            fasta_write.start()

            # Wait for both files to be written before continuing
            while nexus_write.is_alive() or fasta_write.is_alive():
                sleep(0.5)

        def concat_seq():
            # Write concatenated alignments to file in NEXUS format
            alignm_concat.write_nexus_data(filename=out_fn_nucl_concat_nexus)
            log.info(" NEXUS written")

            # Write concatenated alignments to file in FASTA format
            AlignIO.convert(
                out_fn_nucl_concat_nexus, "nexus", out_fn_nucl_concat_fasta, "fasta"
            )
            log.info(" FASTA written")

        log.info(f"concatenate all successful alignments in `{self.user_params.order}` order")

        # sort alignments according to user specification
        self.plastid_data.set_order_map()
        self.success_list.sort(key=lambda t: self.plastid_data.order_map[t[0]])

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

        # Step 3. Write concatenated alignment to file,
        # either synchronously or sequentially, depending on user parameter
        if self.user_params.concat:
            log.info("writing concatenation to file sequentially")
            concat_seq()
        else:
            log.info("writing concatenation to file synchronously")
            concat_sync()

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

    def _evaluate_nuc(self, identifier: str, nuc: Seq, prot: Seq) -> Seq:
        """Returns nucleotide sequence if works (can remove trailing stop)"""
        if len(nuc) % 3:
            log.debug(
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
            debug = (
                f"Inconsistent lengths for {identifier}, ungapped protein {len(p)}, "
                f"tripled {len(p) * 3} vs ungapped nucleotide {len(nuc)}."
            )
            if t.endswith(p):
                debug += f"\nThere are {len(t) - len(p)} extra nucleotides at the start."
            elif t.startswith(p):
                debug += f"\nThere are {len(t) - len(p)} extra nucleotides at the end."
            elif p in t:
                debug += "\nHowever, protein sequence found within translated nucleotides."
            elif p[1:] in t:
                debug += "\nHowever, ignoring first amino acid, protein sequence found within translated nucleotides."
            log.debug(debug)

        if t == p:
            return nuc
        elif p.startswith("M") and "M" + t[1:] == p:
            if str(nuc[0:3]).upper() in ambiguous_generic_by_id[self.table].start_codons:
                return nuc
            else:
                log.debug(
                    f"Translation for {identifier} would match if {nuc[0:3].upper()} "
                    f"was a start codon (check correct table used)"
                )
                log.warning(f"Translation check failed for {identifier}")

        else:
            m = "".join("." if x == y else "!" for (x, y) in zip(p, t))
            if len(prot) < 70:
                log.debug(
                    f"Translation mismatch for {identifier} [0:{len(prot)}]\n"
                    f"Protein:     {p}\n"
                    f"             {m}\n"
                    f"Translation: {t}\n"
                )
            else:
                for offset in range(0, len(p), 60):
                    log.debug(
                        f"Translation mismatch for {identifier} [{offset}:{offset + 60}]\n"
                        f"Protein:     {p[offset:offset + 60]}\n"
                        f"             {m[offset:offset + 60]}\n"
                        f"Translation: {t[offset:offset + 60]}\n"
                    )
            log.warning(f"Translation check failed for {identifier}")

    def _backtrans_seq(self, aligned_protein_record: SeqRecord, unaligned_nucleotide_record: SeqRecord) -> SeqRecord:
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

        # Per https://biopython.org/docs/1.81/api/Bio.Seq.html,
        # `replace` is proper replacement for depreciated method `ungap`.
        try:
            ungapped_protein = aligned_protein_record.seq.replace(self.gap, "")
        except AttributeError:
            ungapped_protein = aligned_protein_record.seq.ungap(self.gap)

        ungapped_nucleotide = unaligned_nucleotide_record.seq
        if self.table:
            ungapped_nucleotide = self._evaluate_nuc(
                aligned_protein_record.id, ungapped_nucleotide, ungapped_protein
            )
        elif len(ungapped_protein) * 3 != len(ungapped_nucleotide):
            log.debug(
                f"Inconsistent lengths for {aligned_protein_record.id}, "
                f"ungapped protein {len(ungapped_protein)}, "
                f"tripled {len(ungapped_protein) * 3} vs "
                f"ungapped nucleotide {len(ungapped_nucleotide)}"
            )
            log.warning(
                f"Backtranslation failed for {aligned_protein_record.id} due to ungapped length mismatch"
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
            log.debug(
                f"Nucleotide sequence for {unaligned_nucleotide_record.id} "
                f"longer than protein {aligned_protein_record.id}"
            )
            log.warning(
                f"Backtranslation failed for {unaligned_nucleotide_record.id} due to unaligned length mismatch"
            )
            return None

        aligned_nuc = unaligned_nucleotide_record[:]  # copy for most annotation
        aligned_nuc.letter_annotation = {}  # clear this
        aligned_nuc.seq = Seq("".join(seq))  # , alpha)  # Modification on 09-Sep-2022 by M. Gruenstaeudl
        if len(aligned_protein_record.seq) * 3 != len(aligned_nuc):
            log.debug(
                f"Nucleotide sequence for {aligned_nuc.id} is length {len(aligned_nuc)} "
                f"but protein sequence {aligned_protein_record.seq} is length {len(aligned_protein_record.seq)} "
                f"(3 * {len(aligned_nuc)} != {len(aligned_protein_record.seq)})"
            )
            log.warning(
                f"Backtranslation failed for {aligned_nuc.id} due to aligned length mismatch"
            )
            return None

        return aligned_nuc

    def _backtrans_seqs(self, protein_alignment: MultipleSeqAlignment, nucleotide_records: Mapping,
                        key_function: Callable[[str], str] = lambda x: x):
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
        self._set_order(args)
        self._set_concat(args)

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
        self.exclude_list = args.excllist
        if self.select_mode == "igs":
            # Excluding matK is necessary, as matK is located inside trnK
            self.exclude_list.append("matK")

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

    def _set_order(self, args: argparse.Namespace):
        order = args.order.lower()
        if order == "alpha":
            self.order = order
        else:
            self.order = "seq"

    def _set_concat(self, args: argparse.Namespace):
        self.concat = args.concat


# -----------------------------------------------------------------#


class PlastidData:
    def __init__(self, user_params: UserParameters):
        self._set_mode(user_params)
        self._set_order_fun(user_params)
        self._set_nucleotides()
        self._set_proteins()
        self.order_map = None
        self._set_files(user_params)

    def _set_mode(self, user_params: UserParameters):
        self.mode = user_params.select_mode

    def _set_nucleotides(self):
        self.nucleotides = IntergenicDict() if self.mode == "igs" else PlastidDict()

    def _set_proteins(self):
        self.proteins = PlastidDict() if self.mode == "cds" else None

    def _set_files(self, user_params: UserParameters):
        self.files = [
            f for f in os.listdir(user_params.in_dir) if f.endswith(user_params.fileext)
        ]

    def _set_order_fun(self, user_params: UserParameters):
        self.order_fun = user_params.order

    @staticmethod
    def _add_plast_dict(pdict1: 'PlastidDict', pdict2: 'PlastidDict'):
        for key in pdict2.keys():
            if key in pdict1.keys():
                pdict1[key].extend(pdict2[key])
            else:
                pdict1[key] = pdict2[key]

    def _add_igs_dict(self, pdict: 'IntergenicDict'):
        self.nucleotides.nuc_inv_map.update(pdict.nuc_inv_map)
        for key in pdict.keys():
            if key in self.nucleotides.keys():
                self.nucleotides[key].extend(pdict[key])
            elif self.nucleotides.nuc_inv_map.get(key) in self.nucleotides.keys():
                continue  # Don't count IGS in the IRs twice
            else:
                self.nucleotides[key] = pdict[key]

    def add_nucleotides(self, pdict: 'PlastidDict'):
        if isinstance(pdict, IntergenicDict):
            self._add_igs_dict(pdict)
        else:
            self._add_plast_dict(self.nucleotides, pdict)

    def add_proteins(self, pdict: 'PlastidDict'):
        if self.proteins is not None:
            self._add_plast_dict(self.proteins, pdict)

    def set_order_map(self):
        order_list = list(self.nucleotides.keys())
        if self.order_fun == "alpha":
            order_list.sort()
        self.order_map = {
            nuc: index for index, nuc in enumerate(order_list)
        }

    def remove_nuc(self, nuc_name: str):
        if nuc_name not in self.nucleotides.keys():
            return

        del self.nucleotides[nuc_name]
        if self.proteins is not None:
            del self.proteins[nuc_name]


class PlastidDict(OrderedDict):
    def __init__(self):
        super().__init__()
        self.plastome_nucs = defaultdict(set)

    def _not_id(self, nuc_name: str, rec_name: str) -> bool:
        rec_set = self.plastome_nucs.get(rec_name)
        return True if not rec_set else nuc_name not in rec_set

    def _is_nuc(self, nuc_name: str) -> bool:
        return nuc_name in self.keys()

    def add_feature(self, feature: 'PlastidFeature'):
        if feature.seq_obj is None:
            log.warning(feature.status_str())
            return

        is_nuc = self._is_nuc(feature.nuc_name)
        not_id = self._not_id(feature.nuc_name, feature.rec_name)
        # if nuc list exists, and this nuc has not been added for the plastome
        if is_nuc and not_id:
            self[feature.nuc_name].append(feature.get_record())
        # if the nuc list does not exist
        elif not is_nuc:
            self[feature.nuc_name] = [feature.get_record()]
        # record that the feature for the plastome has been added
        self.plastome_nucs[feature.rec_name].add(feature.nuc_name)


class IntergenicDict(PlastidDict):
    def __init__(self):
        super().__init__()
        self.nuc_inv_map = {}

    def add_feature(self, igs: 'IntergenicFeature'):
        super().add_feature(igs)
        self.nuc_inv_map[igs.nuc_name] = igs.inv_nuc_name

# -----------------------------------------------------------------#


class PlastidFeature:
    # class fields
    type: str = "Plastid feature"
    default_exception: str = "exception"

    @staticmethod
    def get_gene(feat: SeqFeature) -> Optional[str]:
        return feat.qualifiers["gene"][0] if feat.qualifiers.get("gene") else None

    @staticmethod
    def safe_name(name: Optional[str]) -> Optional[str]:
        if name is None:
            return None

        return sub(
            r"\W", "", name.replace("-", "_")
        )

    @staticmethod
    def get_safe_gene(feat: SeqFeature) -> Optional[str]:
        return PlastidFeature.safe_name(PlastidFeature.get_gene(feat))

    def __init__(self, record: SeqRecord, feature: SeqFeature):
        self._exception = None
        self.seq_obj = None

        self._set_nuc_name(feature)
        self._set_rec_name(record)
        self._set_feature(feature)
        self._set_seq_obj(record)
        self._set_seq_name()

    def _set_nuc_name(self, feature: SeqFeature):
        self.nuc_name = self.get_safe_gene(feature)

    def _set_rec_name(self, record: SeqRecord):
        self.rec_name = record.name

    def _set_feature(self, feature: SeqFeature):
        self.feature = feature

    def _set_seq_obj(self, record: SeqRecord):
        try:
            self.seq_obj = self.feature.extract(record).seq
        except Exception as e:
            self._set_exception(e)

    def _set_seq_name(self):
        self.seq_name = f"{self.nuc_name}_{self.rec_name}"

    def _set_exception(self, exception: Exception = None):
        if exception is None:
            self._exception = self.default_exception
        else:
            self._exception = exception

    def status_str(self) -> str:
        message = f"skipped due to {self._exception}" if self._exception else "successfully extracted"
        return f"{self.type} '{self.nuc_name}' in {self.rec_name} {message}"

    def get_record(self) -> SeqRecord:
        if self.seq_obj:
            return SeqRecord.SeqRecord(
                self.seq_obj, id=self.seq_name, name="", description=""
            )
        return None


class GeneFeature(PlastidFeature):
    type = "Gene"
    default_exception = "ambiguous reading frame"

    def _set_seq_obj(self, record: SeqRecord):
        super()._set_seq_obj(record)
        self._trim_mult_three()

    def _trim_mult_three(self):
        trim_char = len(self.seq_obj) % 3
        if trim_char > 0 and self.seq_obj[:3] == "ATG":
            self.seq_obj = self.seq_obj[:-trim_char]
        elif trim_char > 0:
            self.seq_obj = None
            self._set_exception()


class ProteinFeature(GeneFeature):
    type = "Protein"
    default_exception = "ambiguous gene reading frame"

    def __init__(self, record: SeqRecord = None, feature: SeqFeature = None, gene: GeneFeature = None):
        # if we provide a GeneFeature to the constructor, we can just copy the attributes
        if gene is not None:
            self.__dict__.update(gene.__dict__)
        else:
            super().__init__(record, feature)
        self._set_prot_obj()

    def _set_prot_obj(self):
        if self.seq_obj is None:
            self._set_exception()
            return

        try:
            self.seq_obj = self.seq_obj.translate(table=11)  # Getting error TTA is not stop codon.
        except Exception as e:
            self._set_exception(e)


class IntronFeature(PlastidFeature):
    type = "Intron"

    def __init__(self, record: SeqRecord, feature: SeqFeature, offset: int = 0):
        self.offset = offset
        super().__init__(record, feature)

    def _set_nuc_name(self, feature: SeqFeature):
        super()._set_nuc_name(feature)
        self.nuc_name += "_intron" + str(self.offset + 1)

    def _set_feature(self, feature: SeqFeature):
        super()._set_feature(feature)
        exon_1 = self.feature.location.parts[self.offset]
        exon_2 = self.feature.location.parts[self.offset + 1]
        in_order = exon_2.start >= exon_1.end

        if in_order:
            self.feature.location = FeatureLocation(
                exon_1.end, exon_2.start
            )
        else:
            self.feature.location = FeatureLocation(
                exon_2.start, exon_1.end
            )


class IntergenicFeature(PlastidFeature):
    type = "Intergenic spacer"
    default_exception = "negative intergenic length"

    def __init__(self, record: SeqRecord, current_feat: SeqFeature, subsequent_feat: SeqFeature):
        self._set_sub_feat(subsequent_feat)
        super().__init__(record, current_feat)

    def _set_sub_feat(self, subsequent_feat: SeqFeature):
        self.subsequent_feat = subsequent_feat
        self.subsequent_name = self.get_safe_gene(subsequent_feat)

    def _set_feature(self, current_feat: SeqFeature):
        # Note: It's unclear if +1 is needed here.
        self.start_pos = ExactPosition(current_feat.location.end)  # +1)
        self.end_pos = ExactPosition(self.subsequent_feat.location.start)
        self.feature = FeatureLocation(self.start_pos, self.end_pos) if self.start_pos < self.end_pos else None
        self.default_exception = IntergenicFeature.default_exception + \
                                 f" (start pos: {self.start_pos}, end pos:{self.end_pos})"

    def _set_seq_obj(self, record: SeqRecord):
        if self.feature is None:
            self._set_exception()
        else:
            super()._set_seq_obj(record)

    def _set_seq_name(self):
        self.inv_nuc_name = f"{self.subsequent_name}_{self.nuc_name}"
        self.nuc_name += "_" + self.subsequent_name
        super()._set_seq_name()


# ------------------------------------------------------------------------------#
# MAIN
# ------------------------------------------------------------------------------#
def split_list(input_list: List, num_lists: int) -> List[List]:
    # find the desired number elements in each list
    list_len = len(input_list) // num_lists

    # create the first num_lists-1 lists
    nested_list = [
        input_list[index * list_len: (index + 1) * list_len] for index in range(num_lists - 1)
    ]
    # create the final list, including the division remainder number of elements
    nested_list.append(input_list[(num_lists - 1) * list_len:])
    return nested_list


def setup_logger(user_params: UserParameters) -> logging.Logger:
    logger = logging.getLogger(__name__)
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    log_level = logging.DEBUG if user_params.verbose else logging.INFO
    coloredlogs.install(fmt=log_format, level=log_level, logger=logger)
    logger.debug(f"{parser.prog} {__version__}")
    return logger


def check_dependency(software: str = "mafft"):
    if shutil.which(software) is None:
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
        action="store_true",
        required=False,
        help="(Optional) Enable verbose logging",
    )
    parser.add_argument(
        "--order",
        "-r",
        type=str,
        required=False,
        help="(Optional) Order that the alignments should be saved (`seq` or `alpha`)",
        default="seq",
    )
    parser.add_argument(
        "--concat",
        "-c",
        action="store_true",
        required=False,
        help="(Optional) Enable sequential writing of concatenation files",
    )
    params = UserParameters(parser)
    log = setup_logger(params)
    check_dependency()
    main(params)
# ------------------------------------------------------------------------------#
# EOF
# ------------------------------------------------------------------------------#
