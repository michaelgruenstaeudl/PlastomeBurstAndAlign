from typing import Optional, Dict, Any
from Bio import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition, SeqFeature
from collections import OrderedDict, defaultdict
import os
from re import sub, findall, match

# Package imports
from .logger import logger as log


class PlastidData:
    def __init__(self, user_params: Dict[str, Any]):
        """
        A data structure that primarily stores extracted sequence annotations.
        Additionally, the source files to be extracted and the concatenation order of regions are stored.

        Args:
            user_params: Specifications for how the data structure should be created and maintained.
                These are `select_mode`, `in_dir`, `fileext`, and `order`.
        """
        self._set_mode(user_params)
        self._set_order_fun(user_params)
        self._set_nucleotides()
        self._set_proteins()
        self.order_map = None
        self._set_files(user_params)

    def _set_mode(self, user_params: Dict[str, Any]):
        self.mode = user_params.get("select_mode")

    def _set_nucleotides(self):
        self.nucleotides = IntergenicDict() if self.mode == "igs" else PlastidDict()

    def _set_proteins(self):
        self.proteins = PlastidDict() if self.mode.collects_proteins() else None

    def _set_files(self, user_params: Dict[str, Any]):
        self.files = [
            os.path.join(user_params.get("in_dir"), f) for f in os.listdir(user_params.get("in_dir"))
            if f.endswith(user_params.get("fileext"))
        ]

    def _set_order_fun(self, user_params: Dict[str, Any]):
        self.order = user_params.get("order")

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
            # Don't count IGS in the IRs twice
            elif self.nucleotides.nuc_inv_map.get(key) in self.nucleotides.keys():
                continue
            else:
                self.nucleotides[key] = pdict[key]

    def add_nucleotides(self, pdict: 'PlastidDict'):
        """
        Merges a dictionary of region annotations into the dictionary of nucleotide region annotations
        maintained by the data structure.

        Args:
            pdict: Dictionary of nucleotide region annotations.
        """
        if isinstance(pdict, IntergenicDict):
            self._add_igs_dict(pdict)
        else:
            self._add_plast_dict(self.nucleotides, pdict)
        self.nucleotides.plastome_nucs.update(pdict.plastome_nucs)

    def add_proteins(self, pdict: 'PlastidDict'):
        """
        Merges a dictionary of region annotations into the dictionary of protein region annotations
        maintained by the data structure. If the extraction mode does not include proteins, calling
        this method will do nothing.

        Args:
            pdict: Dictionary of protein region annotations.
        """
        if self.mode.collects_proteins():
            self._add_plast_dict(self.proteins, pdict)
            self.proteins.plastome_nucs.update(pdict.plastome_nucs)

    def set_order_map(self):
        """
        Orders the keys of the nucleotide dictionary according to `order`, creates a map from region name to ordinal
        ordering, and sets this map as `order_map`.
        """
        order_list = list(self.nucleotides.keys())
        if self.order == "alpha":
            order_list.sort()
        self.order_map = {
            nuc: index for index, nuc in enumerate(order_list)
        }

    def remove_nuc(self, nuc_name: str):
        """
        Removes the specified region from the nucleotide dictionary. If the extraction mode includes proteins,
        this will also remove the region from the protein dictionary.

        Args:
            nuc_name: Name of the region to remove from the collected data.
        """
        if nuc_name not in self.nucleotides.keys():
            return

        del self.nucleotides[nuc_name]
        if self.mode.collects_proteins():
            del self.proteins[nuc_name]


class PlastidDict(OrderedDict):
    def __init__(self):
        """
        A dictionary used to store extracted sequence annotations.
        """
        super().__init__()
        self.plastome_nucs = defaultdict(set)

    def _not_id(self, nuc_name: str, rec_name: str) -> bool:
        rec_set = self.plastome_nucs.get(rec_name)
        return True if not rec_set else nuc_name not in rec_set

    def _is_nuc(self, nuc_name: str) -> bool:
        return nuc_name in self.keys()

    def add_feature(self, feature: 'PlastidFeature'):
        """
        Adds an annotation to the dictionary. If an annotation with the same name has already been added
        for the plastome (record name), it will not be added.

        Args:
            feature: A plastid feature annotation.
        """
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
        """
        A dictionary used to store extracted intergenic regions.
        """
        super().__init__()
        self.nuc_inv_map = {}

    def add_feature(self, igs: 'IntergenicFeature'):
        super().add_feature(igs)
        self.nuc_inv_map[igs.nuc_name] = igs.inv_nuc_name


# -----------------------------------------------------------------#


class PlastidFeature:
    _type: str = "Plastid feature"
    _default_exception: str = "exception"

    @staticmethod
    def get_gene(feat: SeqFeature) -> Optional[str]:
        """
        Finds and returns the gene name for a feature. If the feature is not a gene, `None` will be returned.

        Args:
            feat: `Biopython` feature.

        Returns:
            The gene name.

        """
        return feat.qualifiers["gene"][0] if feat.qualifiers.get("gene") else None

    @staticmethod
    def clean_gene(gene: Optional[str], cut_trna: bool = True) -> Optional[str]:
        """
        Standardizes a gene name by case and delimiter usage. By default, it also removes additional
        qualifiers in the name if tRNA. If there is no initial gene name (`None`), `None` will be returned.

        Args:
            gene: Gene name.
            cut_trna: Option to disregard additional qualifiers from the name (default is `TRUE`).

        Returns:
            Cleaned gene name.

        """
        if not gene:
            return None

        # look for standardized gene name
        gene_pattern = r"^[a-zA-Z]{3}[a-zA-Z0-9]+"
        gene_sub = match(gene_pattern, gene)
        # if non-standard, just return a "safe" version of the name
        if not gene_sub:
            return sub(
                r"\W+", "_", gene.replace("-", "_")
            )
        # the first match is the beginning of the cleaned name
        cleaned_gene = gene_sub.group()[0:3].lower() + gene_sub.group()[3].upper()
        if cut_trna:
            return cleaned_gene

        # look for tRNA qualifiers
        trna_pattern = r"(\b\w{3}\b)"
        qual_subs = findall(trna_pattern, gene)
        # if not tRNA, return the gene
        if not qual_subs:
            return cleaned_gene

        # handle each tRNA qualifier appropriately
        codon_pattern = r"\b[ACGTUacgtu]{3}\b"
        for qualifier in qual_subs:
            # insert delimiter
            cleaned_gene += "_"
            # qualifier is codon
            if match(codon_pattern, qualifier):
                cleaned_gene += qualifier.upper()
            # qualifier is amino acid
            else:
                cleaned_gene += qualifier[0].upper() + qualifier[1:4].lower()
        return cleaned_gene

    @staticmethod
    def get_safe_gene(feat: SeqFeature) -> Optional[str]:
        """
        Finds and returns the cleaned gene name for a feature. If the feature is not a gene, `None` will be returned.

        Args:
            feat: `Biopython` feature.

        Returns:
            Cleaned gene name.

        """
        return PlastidFeature.clean_gene(PlastidFeature.get_gene(feat))

    def __init__(self, record: SeqRecord, feature: SeqFeature):
        """
        Primarily a wrapper for a `SeqRecord` that can be accessed by `get_record()`.
        It contains other variables and methods that document the creation of the `SeqRecord`.

        Args:
            record: The record that the feature is being extracted from.
            feature: The feature to be extracted.
        """
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
            self._exception = self._default_exception
        else:
            self._exception = exception

    def status_str(self) -> str:
        """
        Retrieves a string that describes the status of the contained `SeqRecord`.

        Returns:
            A string describing the `SeqRecord` status.

        """
        message = f"skipped due to {self._exception}" if self._exception else "successfully extracted"
        return f"{self._type} '{self.nuc_name}' in {self.rec_name} {message}"

    def get_record(self) -> SeqRecord:
        """
        Instantiates the contained `SeqRecord` and returns it.
        If the `SeqRecord` could not be successfully extracted, `None` is returned.
        Additional information about the `None` return can be accessed by `status_str()`.

        Returns:
            The instantiated `SeqRecord`.

        """
        if self.seq_obj:
            return SeqRecord.SeqRecord(
                self.seq_obj, id=self.seq_name, name="", description=""
            )
        return None


class GeneFeature(PlastidFeature):
    """
    Primarily a wrapper for a gene `SeqRecord` that can be accessed by `get_record()`.
    It contains other variables and methods that document the creation of the `SeqRecord`.
    """

    _type = "Gene"
    _default_exception = "ambiguous reading frame"

    def _set_seq_obj(self, record: SeqRecord):
        super()._set_seq_obj(record)
        self._trim_mult_three()

    def _trim_mult_three(self):
        """
        Accesses if a proper reading frame exists for the gene.
        """
        trim_char = len(self.seq_obj) % 3
        if trim_char > 0 and self.seq_obj[:3] == "ATG":
            self.seq_obj = self.seq_obj[:-trim_char]
        elif trim_char > 0:
            self.seq_obj = None
            self._set_exception()


class ProteinFeature(GeneFeature):
    _type = "Protein"
    _default_exception = "ambiguous gene reading frame"

    def __init__(self, record: SeqRecord = None, feature: SeqFeature = None, gene: GeneFeature = None):
        """
        Primarily a wrapper for a protein `SeqRecord` that can be accessed by `get_record()`.
        It contains other variables and methods that document the creation of the `SeqRecord`.

        If an existing gene feature (`gene`) is provided instead of a `Biopython` `record` and `feature`,
        the object will be copied and the contained nucleotide sequence will be translated to a protein sequence.
        This is how the class is used within the module, but the standard `PlastidFeature` constructor
        exists for use cases where the nucleotide sequence of a feature is not used.

        Args:
            record: The record that the feature is being extracted from (default is `None`).
            feature: The feature to be extracted (default is `None`).
            gene: An existing gene feature (default is `None`).
        """
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
    _type = "Intron"

    def __init__(self, record: SeqRecord, feature: SeqFeature, offset: int = 0):
        """
        Primarily a wrapper for an intron `SeqRecord` that can be accessed by `get_record()`.
        It contains other variables and methods that document the creation of the `SeqRecord`.

        Args:
            record: The record that the feature is being extracted from.
            feature: The feature to be extracted.
            offset: An index offset for which intron is being considered (default is 0).
        """
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
    _type = "Intergenic spacer"
    _default_exception = "negative intergenic length"

    def __init__(self, record: SeqRecord, current_feat: SeqFeature, subsequent_feat: SeqFeature):
        """
        Primarily a wrapper for an intergenic spacer `SeqRecord` that can be accessed by `get_record()`.
        It contains other variables and methods that document the creation of the `SeqRecord`.

        Args:
            record: The record that the feature is being extracted from.
            current_feat: The feature that borders the beginning of the intergenic space.
            subsequent_feat: The feature that borders the end of the intergenic space.
        """
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
        self.default_exception = (
                f"{IntergenicFeature._default_exception}"
                f" (start pos: {self.start_pos}, end pos:{self.end_pos})"
        )

    def _set_seq_obj(self, record: SeqRecord):
        if self.feature is None:
            self._set_exception()
        else:
            super()._set_seq_obj(record)

    def _set_seq_name(self):
        self.inv_nuc_name = f"{self.subsequent_name}_{self.nuc_name}"
        self.nuc_name += "_" + self.subsequent_name
        super()._set_seq_name()
