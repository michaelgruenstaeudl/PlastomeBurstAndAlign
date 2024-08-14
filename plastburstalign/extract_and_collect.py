import bisect
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from typing import Union, List, Tuple, Optional, Any, Dict
from Bio import SeqIO, SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from copy import deepcopy
import multiprocessing
import os

# Package imports
from .plastid_data import (PlastidData, PlastidDict, IntergenicDict,
                           GeneFeature, ProteinFeature, IntergenicFeature,
                           IntronFeature, PlastidFeature)
from .logger import Logger, logger as log
from .helpers import split_list


class ExtractAndCollect:
    def __init__(self, plastid_data: 'PlastidData', user_params: Dict[str, Any]):
        """
        Coordinates the parsing of GenBank flatfiles, the extraction of the contained sequence annotations,
        and the collection of the resulting sequence records in a data structure.

        Args:
            plastid_data: Contains the locations of the files to be parsed,
                and where the extracted records will be stored.
            user_params: Specifications for how the extraction should be performed.
                These are `num_threads`, `out_dir`, `verbose`, and `exclude_list`.
        """
        self.plastid_data = plastid_data
        self.user_params = user_params
        self._set_extract_fun()

    def _set_extract_fun(self):
        if self.plastid_data.mode == "cds":
            self._extract_fun = self._extract_cds
        elif self.plastid_data.mode == "igs":
            self._extract_fun = self._extract_igs
        elif self.plastid_data.mode == "int":
            self._extract_fun = self._extract_int
        else:
            self._extract_fun = None

    def extract(self):
        """
        Parses all GenBank flatfiles of a given folder and extracts
        all sequence annotations of the type specified by the user for each flatfile.
        """
        log.info(
            f"parsing GenBank flatfiles and extracting their sequence annotations using "
            f"{self.user_params.get('num_threads')} processes"
        )

        # Step 0. Extract first genome in list for feature ordering
        if self.plastid_data.order == "seq":
            first_genome = self.plastid_data.files.pop(0)
            nuc_dict, prot_dict = self._extract_recs([first_genome])
            self.plastid_data.add_nucleotides(nuc_dict)
            self.plastid_data.add_proteins(prot_dict)

        # Step 1. Create the data for each worker
        file_lists = split_list(self.plastid_data.files, self.user_params.get("num_threads") * 2)

        # Step 2. Use ProcessPoolExecutor to parallelize extraction
        mp_context = multiprocessing.get_context("fork")  # same method on all platforms
        with ProcessPoolExecutor(max_workers=self.user_params.get("num_threads"), mp_context=mp_context) as executor:
            future_to_recs = [
                executor.submit(self._extract_recs, file_list) for file_list in file_lists
            ]
            for future in as_completed(future_to_recs):
                nuc_dict, prot_dict = future.result()
                self.plastid_data.add_nucleotides(nuc_dict)
                self.plastid_data.add_proteins(prot_dict)

        # Step 3. Stop execution if no nucleotides were extracted
        if not self.plastid_data.nucleotides.items():
            log.critical(f"No items in main dictionary: {self.user_params.get('out_dir')}")
            raise Exception()

    def _extract_recs(self, files: List[str]) -> Tuple['PlastidDict', 'PlastidDict']:
        def extract_rec(file: str):
            try:
                log.info(f" parsing {os.path.basename(file)}")
                record = SeqIO.read(file, "genbank")
                extract_fun(record)
            except Exception as e:
                log.error(" %r generated an exception: %s" % (os.path.basename(file), e))

        # new logger instance for this child process (multiprocessing assumed)
        Logger.reinitialize_logger(self.user_params.get("verbose"))

        # "bind" plastid dictionaries to the record extraction function
        nuc_dict = IntergenicDict() if self.plastid_data.mode == 'igs' else PlastidDict()
        extract_fun = partial(self._extract_fun, nuc_dict=nuc_dict)
        if self.plastid_data.mode.collects_proteins():
            prot_dict = PlastidDict()
            extract_fun = partial(extract_fun, protein_dict=prot_dict)
        else:
            prot_dict = None

        # extract the annotations and collect them in the plastid dictionaries
        for f in files:
            extract_rec(f)
        return nuc_dict, prot_dict

    def _extract_cds(self, rec: SeqRecord, nuc_dict: 'PlastidDict', protein_dict: 'PlastidDict'):
        """
        Extracts all CDS (coding sequences = genes) from a given sequence record
        """
        features = self._cds_features(rec)
        for feature in features:
            # Step 1. Extract nucleotide sequence of each gene and add to dictionary
            gene = GeneFeature(rec, feature)
            nuc_dict.add_feature(gene)

            # Step 2. Translate nucleotide sequence to protein and add to dictionary
            protein = ProteinFeature(gene=gene)
            protein_dict.add_feature(protein)

    def _extract_igs(self, rec: SeqRecord, nuc_dict: 'PlastidDict'):
        """
        Extracts all IGS (intergenic spacers) from a given sequence record
        """
        # Step 1. Extract all genes from record (i.e., cds, trna, rrna)
        # Resulting list contains adjacent features in order of appearance on genome
        all_genes = self._igs_features(rec)
        # Step 2. Loop through genes
        for count, idx in enumerate(range(0, len(all_genes) - 1), 1):
            current_feat = all_genes[idx]
            subsequent_feat = all_genes[idx + 1]

            # Step 3. Make IGS SeqFeature
            igs = IntergenicFeature(rec, current_feat, subsequent_feat)

            # Step 4. Attach IGS to growing dictionary
            nuc_dict.add_feature(igs)

    def _extract_int(self, rec: SeqRecord, nuc_dict: 'PlastidDict'):
        """
        Extracts all INT (introns) from a given sequence record
        """
        features = self._int_features(rec)
        for feature in features:
            # Step 1.a. If one intron in gene:
            if len(feature.location.parts) == 2:
                intron = IntronFeature(rec, feature)
                nuc_dict.add_feature(intron)

            # Step 1.b. If two introns in gene:
            elif len(feature.location.parts) == 3:
                feature_copy = deepcopy(
                    feature
                )  # Important b/c feature is overwritten in extract_internal_intron()

                intron = IntronFeature(rec, feature)
                nuc_dict.add_feature(intron)

                intron = IntronFeature(rec, feature_copy, 1)
                nuc_dict.add_feature(intron)

    def _cds_features(self, record: SeqRecord) -> List[SeqFeature]:
        return [
            f for f in record.features if f.type == "CDS" and self._not_exclude(f)
        ]

    def _igs_features(self, record: SeqRecord) -> List[SeqFeature]:
        # Note: No need to include "if feature.type=='tRNA'", because all tRNAs are also annotated as genes
        genes = [
            f for f in record.features if f.type == "gene" and self._not_exclude(f)
        ]
        # Handle spliced exons
        handler = ExonSpliceHandler(genes, record)
        handler.resolve()
        return genes

    def _int_features(self, record: SeqRecord) -> List[SeqFeature]:
        # Limiting the search to CDS containing introns
        return [
            f for f in record.features if (f.type == "CDS" or f.type == "tRNA") and self._not_exclude(f)
        ]

    def _not_exclude(self, feature: SeqFeature) -> bool:
        gene = PlastidFeature.get_gene(feature)
        return gene and gene not in self.user_params.get("exclude_list") and "orf" not in gene


# -----------------------------------------------------------------#


class ExonSpliceHandler:
    def __init__(self, genes: List[SeqFeature], record: SeqRecord):
        """
        Coordinates the handling of spliced exons.
        Cis and trans-spliced genes are considered separately.

        Args:
            genes: List of gene features.
            record: Record that contains the genes.
        """
        self.genes = genes
        self.record = record

    def resolve(self):
        """
        The handling of compound locations begins with merging the locations of
        cis-spliced genes. During this process, compound location genes qualified as "trans_splicing"
        or containing rps12 are identified and removed from the genes list, as they will be handled next.

        In the second step, the genes identified as trans-spliced are split into multiple simple location
        features. Each of these splits results in a simple location feature for each contiguous location
        of the compound location. All of these simple location features are inserted (or merged) into their
        expected location in the gene list.
        """
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
        """
        Coordinates the handling of cis-spliced exons by merging.

        Args:
            genes: List of gene features.
            record: Record that contains the genes.
        """
        self.genes = genes
        self.rec_name: str = record.name

    def merge(self):
        """
        Locations of compound cis-spliced genes are merged in-place as they occur within the source GenBank file
        (no repositioning of the annotations). Adjacent genes of the same name are assumed to be separately annotated
        cis-spliced genes and merged as well. During this process, compound location annotations qualified as
        "trans_splicing" or containing the gene rps12 (these are known to be trans but are not always qualified as such)
        are removed from the list of gene features, and retained in a separate list.
        """
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
        """
        Coordinates the handling of compound location features using insertion and merging.

        Args:
            genes: List of gene features.
            record: Record that contains the genes.
            compound_features: List of compound location features to be inserted, if already known (optional).
        """
        self.genes = genes
        self.rec_name: str = record.name
        self.compound_features = compound_features if compound_features is not None else None

    def insert(self):
        """
        If `compound_features` is not provided, all compound location features are identified and
        removed from the gene list. Otherwise, the provided `compound_feature` list is used,
        and it is assumed that these correspond to the same record as `genes` and have already been removed
        from this list.

        The genes identified as having compound locations (detailed above) are then split into simple location
        features. Each of these splits results in a simple location feature for each contiguous location
        of the compound location.

        We then examine each simple location and insert it into its expected location
        in the gene list. The expected location is determined by comparing the simple feature's end location
        with the end locations of the features in the gene list. If the expected location has no
        overlap with the proceeding and succeeding genes, and the feature is a different gene from those two,
        it is directly inserted into that expected location. Alternatively, if the expected location of the feature
        results in a flanking gene (strictly adjacent or overlapping) to be the same as the gene to be
        inserted they are merged. The merge can occur on the proceeding side. succeeding side, or both.
        When merging the gene location, the start location is minimized and the end location is maximized.
        """
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
        not_overlap = self._is_after_previous and self._is_before_current
        not_same = not self._is_same_previous and not self._is_same_current
        if not_overlap and not_same:
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


class DataCleaning:
    def __init__(self, plastid_data: 'PlastidData', user_params: Dict[str, Any]):
        """
        Coordinates the cleaning (removal) of dictionary regions based on properties of the regions.

        Args:
            plastid_data: Plastid data to be cleaned.
            user_params: Specifications for how the cleaning process should be performed.
                These are `min_seq_length` and `min_num_taxa`.
        """
        self.plastid_data = plastid_data
        self.user_params = user_params

    def clean(self):
        """
        Cleans the nucleotide and protein dictionaries according to user specifications.
        Specifically, this removes regions that are below the threshold of
        `min_seq_length` or `min_num_taxa`.
        """
        log.info("cleaning extracted sequence annotations")
        log.info(
            f"  removing annotations that occur in fewer than {self.user_params.get('min_num_taxa')} taxa"
        )
        log.info(
            f"  removing annotations whose longest sequence is shorter than {self.user_params.get('min_seq_length')} bp"
        )
        for feat_name, rec_list in list(self.plastid_data.nucleotides.items()):
            self._remove_infreq(feat_name, rec_list)
            self._remove_short(feat_name, rec_list)

    def _remove_short(self, feat_name: str, rec_list: List[SeqRecord]):
        longest_seq = max([len(s.seq) for s in rec_list])
        if longest_seq < self.user_params.get("min_seq_length"):
            log.info(f"    removing {feat_name} for not reaching the minimum sequence length defined")
            self.plastid_data.remove_nuc(feat_name)

    def _remove_infreq(self, feat_name: str, rec_list: List[SeqRecord]):
        if len(rec_list) < self.user_params.get("min_num_taxa"):
            log.info(f"    removing {feat_name} for not reaching the minimum number of taxa defined")
            self.plastid_data.remove_nuc(feat_name)
