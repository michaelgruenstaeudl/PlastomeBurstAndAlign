import bisect
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from typing import List, Tuple, Optional, Any, Dict, NamedTuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from copy import deepcopy
import multiprocessing
import os

# Package imports
from .seqfeature_ops import (PlastidData, PlastidDict, IntergenicDict,
                             GeneFeature, ProteinFeature, IntergenicFeature,
                             IntronFeature, PlastidFeature)
from .logging_ops import Logger, logger as log
from .helpers import split_list


class ExtractAndCollect:
    def __init__(self, plastid_data: PlastidData, user_params: Dict[str, Any]):
        """
        Coordinates the parsing of GenBank flatfiles, the extraction of the contained sequence annotations,
        and the collection of the resulting sequence records in a data structure.

        Args:
            plastid_data: Contains the locations of the files to be parsed,
                and where the extracted records will be stored.
            user_params: Specifications for how the extraction should be performed.
                These are `num_threads`, `out_dir`, `verbose`, and `exclude_cds`.
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
            log.critical("Undefined extraction mode")
            raise ValueError

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

    def _extract_recs(self, files: List[str]) -> Tuple[PlastidDict, PlastidDict]:
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

    def _extract_cds(self, rec: SeqRecord, nuc_dict: PlastidDict, protein_dict: PlastidDict):
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

    def _extract_igs(self, rec: SeqRecord, nuc_dict: PlastidDict):
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

    def _extract_int(self, rec: SeqRecord, nuc_dict: PlastidDict):
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
        return gene and gene not in self.user_params.get("exclude_cds") and "orf" not in gene


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
        log.info(f"   resolving trans-spliced genes in {self.record.name}")
        insertor = ExonSpliceInsertor(self.genes, self.record, merger.trans_list)
        insertor.insert()


class ExonSpliceMerger:
    class FeatureTuple(NamedTuple):
        feature: SeqFeature
        gene: Optional[str]
        is_trans: bool
        is_compound: bool

    def __init__(self, genes: List[SeqFeature], record: SeqRecord):
        """
        Coordinates the handling of cis-spliced exons by merging.

        Args:
            genes: List of gene features.
            record: Record that contains the genes.
        """
        self.genes = genes
        self.rec_name: str = record.name
        self.trans_list: List[SeqFeature] = []

    def merge(self):
        """
        Locations of compound cis-spliced genes are merged in-place as they occur within the source GenBank file
        (no repositioning of the annotations). Adjacent genes of the same name are assumed to be separately annotated
        cis-spliced genes and merged as well. During this process, compound location annotations qualified as
        "trans_splicing" or containing the gene rps12 (these are known to be trans but are not always qualified as such)
        are removed from the list of gene features, and retained in a separate list.
        """
        end_index = len(self.genes) - 1
        subsequent = self._get_feat_tuple(end_index)
        for index in range(end_index - 1, -1, -1):
            current = self._get_feat_tuple(index)

            if current.is_trans:
                self.trans_list.append(current.feature)
                del self.genes[index]
            elif current.is_compound:
                self._merge_cis_exons(current)
                # print new alignment
                self._print_align(current, subsequent)

            is_same_gene = current.gene == subsequent.gene
            if is_same_gene and not current.is_trans:
                self._merge_adj_exons(current, subsequent)
                # delete merged exon
                del self.genes[index + 1]
                # update new subsequent feature
                subsequent = self._get_feat_tuple(index + 1)
                # print new alignment
                self._print_align(current, subsequent)

            # update for next iteration
            subsequent = current

    def _get_feat_tuple(self, index: int) -> FeatureTuple:
        feature = self.genes[index]
        gene = PlastidFeature.get_safe_gene(feature)
        is_trans = self._get_is_trans(feature, gene)
        is_compound = type(feature.location) is CompoundLocation
        feat_tuple = self.FeatureTuple(feature, gene, is_trans, is_compound)
        return feat_tuple

    @staticmethod
    def _get_is_trans(feature: SeqFeature, gene: str) -> bool:
        return bool(feature.qualifiers.get("trans_splicing")) or gene == "rps12"

    @staticmethod
    def _merge_cis_exons(feature_tuple: FeatureTuple):
        loc_parts = feature_tuple.feature.location.parts
        gene_start = min(p.start for p in loc_parts)
        gene_end = max(p.end for p in loc_parts)
        feature_tuple.feature.location = FeatureLocation(gene_start, gene_end)

    @staticmethod
    def _merge_adj_exons(current: FeatureTuple, subsequent: FeatureTuple):
        gene_start = min(current.feature.location.start, subsequent.feature.location.start)
        gene_end = max(current.feature.location.end, subsequent.feature.location.end)
        current.feature.location = FeatureLocation(gene_start, gene_end)

    def _print_align(self, current: FeatureTuple, subsequent: FeatureTuple):
        log.debug(
            f"   Merging exons of {current.gene} within {self.rec_name}\n"
            "-----------------------------------------------------------\n"
            f"\t{current.gene}\t\t\t{subsequent.gene}\n"
            f"\t{current.feature.location}\t\t{subsequent.feature.location}\n"
            "-----------------------------------------------------------\n"
        )


class ExonSpliceInsertor:
    class FeatureTuple(NamedTuple):
        feature: Optional[SeqFeature]
        gene: Optional[str]
        location: Optional[FeatureLocation]

    class TestsTuple(NamedTuple):
        is_same_previous: bool
        is_same_current: bool
        is_after_previous: bool
        is_before_current: bool
        not_overlap: bool
        not_same: bool

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
        self.compound_features = compound_features
        self._end_positions = None

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

        self.genes.sort(key=lambda gene: gene.location.end)
        self._end_positions = [
            f.location.end for f in self.genes
        ]
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

    def _create_simple(self):
        self.simple_features = []
        for f in self.compound_features:
            self.simple_features.extend(
                SeqFeature(location=p, type=f.type, id=f.id, qualifiers=f.qualifiers) for p in f.location.parts
            )
        self.simple_features.sort(key=lambda feat: feat.location.end, reverse=True)

    def _insert_simple(self):
        """
        Insert the simple features at the correct indices in the gene list if applicable
        """
        for insert_feature in self.simple_features:
            # create a structure that can be compared with adjacent features
            insert = self.FeatureTuple(
                insert_feature,
                PlastidFeature.get_safe_gene(insert_feature),
                insert_feature.location
            )

            # find insertion index
            insert_index = bisect.bisect_left(self._end_positions, insert.location.end)

            # set appropriate adjacent features
            previous = self._get_feat_tuple(insert_index - 1)
            current = self._get_feat_tuple(insert_index)

            # checks for how to handle the insert feature
            tests = self._get_adj_tests(previous, insert, current)

            # directly reposition if possible
            if tests.not_overlap and tests.not_same:
                self._insert_at_index(insert_index, insert)
                message = f"Repositioning {insert.gene} within {self.rec_name}"
                self._print_align(previous, insert, current, message)
            # attempt to resolve by merging
            else:
                self._try_merging(insert_index, previous, insert, current, tests)

    def _get_feat_tuple(self, index: int) -> FeatureTuple:
        feature = None if index == 0 or index == len(self.genes) else self.genes[index]
        gene = None if not feature else PlastidFeature.get_safe_gene(feature)
        location = None if not feature else feature.location
        feat_tuple = self.FeatureTuple(
            feature,
            gene,
            location
        )
        return feat_tuple

    def _get_adj_tests(self, previous: FeatureTuple, insert: FeatureTuple, current: FeatureTuple) -> TestsTuple:
        is_same_previous = previous.gene == insert.gene
        is_same_current = current.gene == insert.gene
        is_after_previous = not previous.location or previous.location.end < insert.location.start
        is_before_current = not current.location or insert.location.end < current.location.start
        not_overlap = is_after_previous and is_before_current
        not_same = not is_same_previous and not is_same_current
        tests_tuple = self.TestsTuple(
            is_same_previous,
            is_same_current,
            is_after_previous,
            is_before_current,
            not_overlap,
            not_same
        )
        return tests_tuple

    def _try_merging(self, insert_index: int, previous: FeatureTuple, insert: FeatureTuple, current: FeatureTuple,
                     tests: TestsTuple):
        is_merged = False
        start = insert.location.start
        end = insert.location.end

        # if insert and current feature are the same gene, and insert starts before current,
        # remove current feature, and update ending location
        if tests.is_same_current and insert.location.start < current.location.start:
            end = current.location.end
            self._remove_at_index(insert_index)
            is_merged = True

        # if insert and previous feature are the same gene,
        # use the smaller start location, and remove previous feature
        if tests.is_same_previous:
            start = min(start, previous.location.start)
            # elements in list will shift to left, so update index
            insert_index -= 1
            self._remove_at_index(insert_index)
            is_merged = True

        # perform insertion if needed
        if is_merged:
            # updated feature to insert
            insert_feature = SeqFeature(
                location=FeatureLocation(start, end, insert.location.strand),
                type=insert.feature.type, id=insert.feature.id, qualifiers=insert.feature.qualifiers
            )
            insert = self.FeatureTuple(
                insert_feature,
                insert.gene,
                insert_feature.location
            )
            # new adjacent features
            previous = self._get_feat_tuple(insert_index - 1)
            current = self._get_feat_tuple(insert_index)
            # insert
            self._insert_at_index(insert_index, insert)
            message = f"Merging exons of {insert.gene} within {self.rec_name}"
            self._print_align(previous, insert, current, message)

    def _insert_at_index(self, insert_index: int, insert: FeatureTuple):
        self.genes.insert(insert_index, insert.feature)
        self._end_positions.insert(insert_index, insert.location.end)

    def _remove_at_index(self, index: int):
        del self.genes[index]
        del self._end_positions[index]

    @staticmethod
    def _print_align(previous: FeatureTuple, insert: FeatureTuple, current: FeatureTuple, message: str):
        previous_gene = "\t" if not previous.gene else previous.gene
        current_gene = "" if not current.gene else current.gene

        previous_loc = "\t\t\t\t\t" if not previous.location else previous.location
        current_loc = "" if not current.location else current.location

        log.debug(
            f"   {message}\n"
            "-----------------------------------------------------------\n"
            f"\t\t{previous_gene}\t\t\t\t{insert.gene}\t\t\t\t{current_gene}\n"
            f"\t{previous_loc}\t{insert.location}\t{current_loc}\n"
            "-----------------------------------------------------------\n"
        )


# -----------------------------------------------------------------#


class DataCleaning:
    def __init__(self, plastid_data: PlastidData, user_params: Dict[str, Any]):
        """
        Coordinates the cleaning (removal) of dictionary regions based on properties of the regions.

        Args:
            plastid_data: Plastid data to be cleaned.
            user_params: Specifications for how the cleaning process should be performed.
                These are `min_seq_length`, `min_num_taxa`, and `exclude_region`.
        """
        self.plastid_data = plastid_data
        self.min_num_taxa = user_params.get("min_num_taxa")
        self.min_seq_length = user_params.get("min_seq_length")
        self.exclude_region = user_params.get("exclude_region")

    def clean(self):
        """
        Cleans the nucleotide and protein dictionaries according to user specifications.
        Specifically, this removes regions that are below the threshold of
        `min_seq_length` or `min_num_taxa`.
        Additionally, any regions specified in `exclude_region` are removed.
        """
        log.info("cleaning extracted sequence annotations")
        if self.exclude_region:
            log.info(
                f"  removing excluded regions"
            )
        log.info(
            f"  removing annotations that occur in fewer than {self.min_num_taxa} taxa"
        )
        log.info(
            f"  removing annotations whose longest sequence is shorter than {self.min_seq_length} bp"
        )
        for feat_name, rec_list in list(self.plastid_data.nucleotides.items()):
            self._remove_excluded(feat_name)
            self._remove_infreq(feat_name, rec_list)
            self._remove_short(feat_name, rec_list)

    def _remove_short(self, feat_name: str, rec_list: List[SeqRecord]):
        longest_seq = max([len(s.seq) for s in rec_list])
        if longest_seq < self.min_seq_length:
            log.info(f"    removing {feat_name} for not reaching the minimum sequence length defined")
            self.plastid_data.remove_nuc(feat_name)

    def _remove_infreq(self, feat_name: str, rec_list: List[SeqRecord]):
        if len(rec_list) < self.min_num_taxa:
            log.info(f"    removing {feat_name} for not reaching the minimum number of taxa defined")
            self.plastid_data.remove_nuc(feat_name)

    def _remove_excluded(self, feat_name: str):
        if feat_name in self.exclude_region:
            log.info(f"    removing {feat_name} for being in exclusion list")
            self.plastid_data.remove_nuc(feat_name)
