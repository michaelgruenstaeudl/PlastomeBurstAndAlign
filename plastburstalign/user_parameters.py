import argparse
import multiprocessing
import os
from enum import auto, StrEnum
from typing import List, Union, Any, Dict, Optional

# Package imports
from .logger import Logger, logger as log


class SelectMode(StrEnum):
    """
    An enumeration class used to keep track of selection and alignment modes.
    """
    CDS = auto()
    INT = auto()
    IGS = auto()

    def collects_proteins(self):
        return self == type(self).CDS


class UserParameters(Dict[str, Any]):
    def __init__(
            self,
            in_dir: str = "./input",
            out_dir: str = "./output",
            select_mode: str = "cds",
            fileext: str = ".gb",
            exclude_list: Optional[List[str]] = None,
            min_seq_length: int = 3,
            min_num_taxa: int = 2,
            num_threads: Union[str, int] = "auto",
            verbose: bool = False,
            order: str = "seq",
            concat: bool = False
    ):
        """
        A dictionary that validates user parameters that are to be used within the package.


        Args:
            in_dir: Path to input directory (which contains the GenBank files
            out_dir: Path to output directory
            select_mode: Type of regions to be extracted (i.e. `cds`, `int`, or `igs`)
            fileext: File extension of input files
            exclude_list: List of genes to be excluded
            min_seq_length: Minimal sequence length (in bp) below which regions will not be extracted
            min_num_taxa: Minimum number of taxa in which a region must be present to be extracted
            num_threads: Number of CPUs to use; can be any positive integer or 'auto'
            verbose: Enable verbose logging
            order: Order that the alignments should be saved ('seq' or 'alpha')
            concat: Enable sequential writing of concatenation files
        """
        super().__init__()
        self._set_select_mode(select_mode)
        self._set_in_dir(in_dir)
        self._set_out_dir(out_dir)
        self._set_fileext(fileext)
        self._set_exclude_list(exclude_list)
        self._set_min_seq_length(min_seq_length)
        self._set_min_num_taxa(min_num_taxa)
        self._set_num_threads(num_threads)
        self._set_verbose(verbose)
        self._set_order(order)
        self._set_concat(concat)

    def _set_select_mode(self, select_mode: str):
        try:
            self["select_mode"] = SelectMode(select_mode.lower())
        except ValueError:
            log.critical(f"Invalid selection mode provided: `{select_mode}`")
            raise ValueError

    def _set_in_dir(self, in_dir: str):
        if not os.path.exists(in_dir) or not isinstance(in_dir, str):
            log.critical(f"Input directory `{in_dir}` does not exist.")
            raise FileNotFoundError
        self["in_dir"] = in_dir

    def _set_out_dir(self, out_dir: str):
        if not os.path.exists(out_dir) or not isinstance(out_dir, str):
            log.critical(f"Output directory `{out_dir}` does not exist.")
            raise FileNotFoundError
        self["out_dir"] = out_dir

    def _set_fileext(self, fileext: str):
        self._test_type("fileext", fileext, str)
        self["fileext"] = fileext

    def _set_exclude_list(self, exclude_list: Optional[List[str]]):
        if not exclude_list:
            exclude_list = ["rps12"]
        self._test_type("exclude_list", exclude_list, list)
        self["exclude_list"] = exclude_list
        if self["select_mode"] == "igs":
            # Excluding matK is necessary, as matK is located inside trnK
            self["exclude_list"].append("matK")

    def _set_min_seq_length(self, min_seq_length: int):
        self._test_type("min_seq_length", min_seq_length, int)
        self["min_seq_length"] = min_seq_length

    def _set_min_num_taxa(self, min_num_taxa: int):
        self._test_type("min_num_taxa", min_num_taxa, int)
        self["min_num_taxa"] = min_num_taxa

    def _set_num_threads(self, num_threads: Union[str, int]):
        if isinstance(num_threads, str) and not num_threads.isnumeric():
            try:
                num_threads = os.cpu_count()
            except NotImplementedError:
                num_threads = multiprocessing.cpu_count()
        elif isinstance(num_threads, str):
            num_threads = int(num_threads)
        self["num_threads"] = num_threads

    def _set_verbose(self, verbose: bool):
        self._test_type("verbose", verbose, bool)
        self["verbose"] = verbose
        Logger.set_logger(verbose)

    def _set_order(self, order: str):
        self._test_type("order", order, str)
        if order == "alpha":
            self["order"] = order
        else:
            self["order"] = "seq"

    def _set_concat(self, concat: bool):
        self._test_type("concat", concat, bool)
        self["concat"] = concat

    @staticmethod
    def _test_type(param_name: str, param_val: Any, expected_type: Any):
        if not isinstance(param_val, expected_type):
            log.critical(f"Invalid `{param_name}` type provided: {type(param_val)}")
            raise TypeError


class UserParametersScript(UserParameters):
    _construct_name_map: Dict[str, str] = {
        "inpd": "in_dir",
        "outd": "out_dir",
        "selectmode": "select_mode",
        "fileext": "fileext",
        "excllist": "exclude_list",
        "minseqlength": "min_seq_length",
        "minnumtaxa": "min_num_taxa",
        "numthreads": "num_threads",
        "verbose": "verbose",
        "order": "order",
        "concat": "concat"
    }

    def __init__(self):
        """
        An interface to create `UserParameters` within the context of a script execution.
        """
        parser = self._get_parser()
        args = parser.parse_args()
        args_dict = self._get_construct_dict(args)
        super().__init__(**args_dict)

    @staticmethod
    def _get_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--inpd",
            "-i",
            type=str,
            required=False,
            help="Path to input directory (which contains the GenBank files)",
        )
        parser.add_argument(
            "--outd",
            "-o",
            type=str,
            required=False,
            help="(Optional) Path to output directory",
        )
        parser.add_argument(
            "--selectmode",
            "-s",
            type=str,
            required=False,
            help="(Optional) Type of regions to be extracted (i.e. `cds`, `int`, or `igs`)",
        )
        parser.add_argument(
            "--fileext",
            "-f",
            type=str,
            required=False,
            help="(Optional) File extension of input files",
        )
        parser.add_argument(
            "--excllist",
            "-e",
            type=str,
            nargs="+",
            required=False,
            help="(Optional) List of genes to be excluded",
        )
        parser.add_argument(
            "--minseqlength",
            "-l",
            type=int,
            required=False,
            help="(Optional) Minimal sequence length (in bp) below which regions will not be extracted",
        )
        parser.add_argument(
            "--minnumtaxa",
            "-t",
            type=int,
            required=False,
            help="(Optional) Minimum number of taxa in which a region must be present to be extracted",
        )
        parser.add_argument(
            "--numthreads",
            "-n",
            type=int,
            required=False,
            help="(Optional) Number of CPUs to use; can be any positive integer or 'auto' (default)",
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
            help="(Optional) Order that the alignments should be saved ('seq' or 'alpha')",
        )
        parser.add_argument(
            "--concat",
            "-c",
            action="store_true",
            required=False,
            help="(Optional) Enable sequential writing of concatenation files",
        )
        return parser

    @classmethod
    def _get_construct_dict(cls, args: argparse.Namespace) -> Dict[str, Any]:
        # map the cmd argument names to class constructor parameters,
        # and remove empty argument values (`None`)
        args_dict = vars(args)
        args_dict = {
            cls._construct_name_map[k]: v for k, v in args_dict.items() if v
        }
        return args_dict
