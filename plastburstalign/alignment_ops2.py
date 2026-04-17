import glob
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from io import StringIO
import time 
from typing import List, Callable, Tuple, Mapping, Optional, Dict, Any

import requests
from Bio import SeqIO, Nexus, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import ambiguous_generic_by_id
from Bio.Seq import Seq
# from Bio.Align.Applications import MafftCommandline
# Package imports
from .helpers import split_list
from .logging_ops import Logger, logger as log
from .seqfeature_ops import PlastidData


# ===================================================================
# Module-level helper functions for multiprocessing
# ===================================================================

def _save_region_unaligned(args: tuple) -> None:
    """
    Helper function to save one region's unaligned sequences to file.
    Must be at module level to be pickleable for multiprocessing.
    
    Args:
        args: Tuple of (feat_name, nucleotides_data, proteins_data, out_dir, collect_proteins)
    """
    feat_name, nucleotides_data, proteins_data, out_dir, collect_proteins = args
    try:
        out_fn_unalign_nucl = os.path.join(out_dir, f"nucl_{feat_name}.unalign.fasta")
        with open(out_fn_unalign_nucl, "w") as hndl:
            SeqIO.write(nucleotides_data, hndl, "fasta")
        
        if collect_proteins and proteins_data:
            out_fn_unalign_prot = os.path.join(out_dir, f"prot_{feat_name}.unalign.fasta")
            with open(out_fn_unalign_prot, "w") as hndl:
                SeqIO.write(proteins_data, hndl, "fasta")
    except Exception as e:
        log.error(f"Error saving {feat_name}: {e}")


def _single_nuc_MSA(args: tuple) -> None:
    """
    Module-level helper for nucleotide MSA alignment (pickleable for multiprocessing).
    
    Args:
        args: Tuple of (feat_name, out_dir, mafft_exec_path, num_threads)
    """
    feat_name, out_dir, mafft_exec_path, num_threads = args
    log.debug(f" aligning {feat_name}")
    out_fn_unalign_nucl = os.path.join(out_dir, f"nucl_{feat_name}.unalign.fasta")
    fasta_dir = os.path.join(out_dir, "fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    out_fn_aligned_nucl = os.path.join(fasta_dir, f"nucl_{feat_name}.aligned.fasta")
    
    if not mafft_exec_path:
        log.error(f"  skipping alignment of {feat_name}; MAFFT not configured")
        return
    
    mafft_cmd = [mafft_exec_path, "--adjustdirection", out_fn_unalign_nucl]
    with open(out_fn_aligned_nucl, 'w') as hndl, open(os.devnull, 'w') as devnull:
        process = subprocess.Popen(mafft_cmd, stdout=hndl, stderr=devnull, text=True)
        process.wait()


def _single_prot_MSA(args: tuple) -> None:
    """
    Module-level helper for protein MSA alignment + back-translation (pickleable for multiprocessing).
    
    Args:
        args: Tuple of (feat_name, out_dir, mafft_exec_path, num_threads, backtrans_table)
    """
    feat_name, out_dir, mafft_exec_path, num_threads, backtrans_table = args
    log.debug(f" aligning {feat_name}")
    out_fn_unalign_prot = os.path.join(out_dir, f"prot_{feat_name}.unalign.fasta")
    out_fn_aligned_prot = os.path.join(out_dir, f"prot_{feat_name}.aligned.fasta")
    out_fn_unalign_nucl = os.path.join(out_dir, f"nucl_{feat_name}.unalign.fasta")
    fasta_dir = os.path.join(out_dir, "fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    out_fn_aligned_nucl = os.path.join(fasta_dir, f"nucl_{feat_name}.aligned.fasta")
    
    if not mafft_exec_path:
        log.error(f"  skipping alignment of {feat_name}; MAFFT not configured")
        return
    
    # Align proteins
    mafft_cmd = [mafft_exec_path, "--adjustdirection", out_fn_unalign_prot]
    with open(out_fn_aligned_prot, 'w') as hndl, open(os.devnull, 'w') as devnull:
        process = subprocess.Popen(mafft_cmd, stdout=hndl, stderr=devnull, text=True)
        process.wait()
    
    # Back-translate
    try:
        translator = BackTranslation(
            feat_name,
            out_fn_aligned_prot,
            out_fn_unalign_nucl,
            out_fn_aligned_nucl,
            backtrans_table
        )
        translator.backtranslate()
    except Exception as err:
        log.warning(f"Unable to conduct back-translation of `{feat_name}`. Error message: {err}.")


class AlignmentCoordination:
    def __init__(self, plastid_data: PlastidData, user_params: Dict[str, Any]):
        """
        Coordinates the alignment of nucleotide or protein sequences.

        Args:
            plastid_data: Plastid data to be aligned.
            user_params: Specifications for how the alignment should be performed.
                These are `num_threads`, `out_dir`, `verbose`, and `concat`.
        """
        self.plastid_data = plastid_data
        self.user_params = user_params
        self.success_list = []
        self.mafft = MAFFT(user_params)
        self._unaligned_saved = False
        self._aligned_saved = False



   
    def save_unaligned(self):
        """
        For each region in the plastid data, an unaligned nucleotide matrix of sequences is saved to file.
        These are saved in FASTA format, located in the path specified by `user_params.out_dir`.
        If the plastid data contains protein sequences, those are saved as unaligned protein matrices.
        """
        log.info("saving individual regions as unaligned nucleotide matrices")
        
        num_workers = self.user_params.get("num_threads")
        mp_context = multiprocessing.get_context()  # same on all platforms
        
        # Prepare arguments for each region
        args_list = []
        for feat_name in self.plastid_data.nucleotides.keys():
            nuc_data = self.plastid_data.nucleotides.get(feat_name)
            prot_data = self.plastid_data.proteins.get(feat_name) if self.plastid_data.mode.collects_proteins() else None
            args_list.append((feat_name, nuc_data, prot_data, self.user_params.get("out_dir"), self.plastid_data.mode.collects_proteins()))
        
        # Use ProcessPoolExecutor to parallelize
        with ProcessPoolExecutor(max_workers=num_workers, mp_context=mp_context) as executor:
            futures = [executor.submit(_save_region_unaligned, args) for args in args_list]
            for future in futures:
                future.result()  # Wait for completion and catch exceptions
        
        self._unaligned_saved = True

    def perform_MSA(self):
        """
        For each region in the plastid data, the sequences are aligned and saved to file as a nucleotide matrix.
        These are saved in FASTA format, located in the path specified by `user_params.out_dir`.

        In the nucleotide case, this is the direct process of aligning the unaligned matrix of each region
        using a third-party software tool. For the protein case, each unaligned matrix is aligned, and
        then back-translated to a nucleotide matrix before saving to file.
        """
        # if not self.mafft.exec_path:
        #     log.error("unable to align regions; alignment tool is not properly configured")
        #     return
        if not self._unaligned_saved:
            self.save_unaligned()

        log.info("conducting the alignment of extracted sequences")
        if self.plastid_data.mode.collects_proteins():
            self._prot_MSA()
        else:
            self._nuc_MSA()
        self._aligned_saved = True

    def _nuc_MSA(self):
        log.info(
            f" conducting multiple sequence alignments based on nucleotide sequence data using "
            f"{self.user_params.get('num_threads')} threads"
        )

        # Step 1. Prepare arguments for each region
        num_workers = self.user_params.get("num_threads")
        args_list = [
            (feat_name, self.user_params.get("out_dir"), self.mafft.exec_path, num_workers)
            for feat_name in self.plastid_data.nucleotides.keys()
        ]

        # Step 2. Use ProcessPoolExecutor to parallelize alignment
        if args_list:
            mp_context = multiprocessing.get_context("fork")
            with ProcessPoolExecutor(max_workers=num_workers, mp_context=mp_context) as executor:
                futures = [executor.submit(_single_nuc_MSA, args) for args in args_list]
                for future in as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        log.error(f"Exception during alignment: {e}")
        else:
            log.critical("No items in nucleotide main dictionary to process")
            raise Exception()

    def _prot_MSA(self):
        log.info(
            f" conducting multiple sequence alignments based on protein sequence data, "
            f"followed by back-translation to nucleotides using {self.user_params.get('num_threads')} threads"
        )

        # Step 1. Prepare arguments for each region
        num_workers = self.user_params.get("num_threads")
        args_list = [
            (feat_name, self.user_params.get("out_dir"), self.mafft.exec_path, num_workers, 11)
            for feat_name in self.plastid_data.nucleotides.keys()
        ]

        # Step 2. Use ProcessPoolExecutor to parallelize alignment and back-translation
        if args_list:
            mp_context = multiprocessing.get_context("fork")
            with ProcessPoolExecutor(max_workers=num_workers, mp_context=mp_context) as executor:
                futures = [executor.submit(_single_prot_MSA, args) for args in args_list]
                for future in as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        log.error(f"Exception during alignment: {e}")

    def perform_MSAs(self):
        """
        For each region in the plastid data, the sequences are aligned and saved to file in FASTA and NEXUS format
        in the path specified by `user_params.out_dir`.
        """
        if not self._aligned_saved:
            self.perform_MSA()

        # log.info(
        #     f"collecting all successful alignments using {self.user_params.get('num_threads')} processes"
        # )
        log.info(
            f" converting all successful alignments to NEXUS format using {self.user_params.get('num_threads')} threads"
        )
        msa_lists = split_list(list(self.plastid_data.nucleotides.keys()), self.user_params.get("num_threads") * 2)
        # Use ThreadPoolExecutor for I/O-bound file conversions
        with ThreadPoolExecutor(max_workers=self.user_params.get("num_threads")) as executor:
            future_to_success = [
                executor.submit(self.convert_MSAs, msa_list)
                for msa_list in msa_lists
            ]
            for future in as_completed(future_to_success):
                success_list = future.result()
                if len(success_list) > 0:
                    self.success_list.extend(success_list)

   
    def convert_MSAs(self, msa_list: List[str]) -> List[Tuple[str, str]]:
        """Convert each region's aligned FASTA file to NEXUS format and collect successful conversions."""
        def _convert_FASTA_to_NEXUS(msa_name: str) -> Optional[str]:
            fasta_dir = os.path.join(self.user_params.get("out_dir"), "fasta")
            aligned_nucl_fasta = os.path.join(fasta_dir, f"nucl_{msa_name}.aligned.fasta")
            nexus_dir = os.path.join(self.user_params.get("out_dir"), "nexus")
            os.makedirs(nexus_dir, exist_ok=True)
            aligned_nucl_nexus = os.path.join(nexus_dir, f"nucl_{msa_name}.aligned.nexus")
            try:
                AlignIO.convert(
                    aligned_nucl_fasta,
                    "fasta",
                    aligned_nucl_nexus,
                    "nexus",
                    molecule_type="DNA",
                )
                return aligned_nucl_nexus
            except Exception:
                log.warning(f"Unable to convert alignment of `{msa_name}` from FASTA to NEXUS.")
                return None

        success_list = []
        for msa_name in msa_list:
            nexus_path = _convert_FASTA_to_NEXUS(msa_name)
            if nexus_path is not None:
                success_list.append((msa_name, nexus_path))
        return success_list

    def concat_MSAs(self):
        """
        Concatenated region alignments are saved to file in both NEXUS and FASTA formats,
        located in the path specified by `user_params.out_dir`.
        """
        # if not self.success_list:
        #     self.collect_MSAs()

        def _concat_sync():
            # Write concatenated alignments to file in NEXUS format
            mp_context = multiprocessing.get_context()
            nexus_write = mp_context.Process(target=alignm_concat.write_nexus_data,
                                             kwargs={"filename": out_fn_nucl_concat_nexus})
            nexus_write.start()

            # Write concatenated alignments to file in FASTA format
            fasta_write = mp_context.Process(target=alignm_concat.export_fasta,
                                             kwargs={"filename": out_fn_nucl_concat_fasta})
            fasta_write.start()

            # Wait for both files to be written before continuing
            nexus_write.join()
            fasta_write.join()
            if nexus_write.exitcode != 0 or fasta_write.exitcode != 0:
                log.error("Child writer failed (nexus=%s, fasta=%s)", nexus_write.exitcode, fasta_write.exitcode)
                raise RuntimeError("Child process failure while writing concatenated alignments")
           
        def _load_nexus_alignment(msa_name: str, nexus_path: str) -> Nexus.Nexus.Nexus:
            # Load the Nexus alignment and rename sequence IDs for concatenation
            alignm_nexus = AlignIO.read(nexus_path, "nexus")
            hndl = StringIO()
            AlignIO.write(alignm_nexus, hndl, "nexus")
            nexus_string = hndl.getvalue()

            gene_name_pattern = fr"\n(_R_)?{msa_name}_"
            nexus_string = re.sub(gene_name_pattern, "\nconcat_", nexus_string)
            return Nexus.Nexus.Nexus(nexus_string)

        log.info(f" concatenating all successful alignments in `{self.plastid_data.order}` order")
        # sort alignments according to user specification
        self.plastid_data.set_order_map()
        self.success_list.sort(key=lambda t: self.plastid_data.order_map[t[0]])

        # Load Nexus alignments with renamed sequence IDs
        self.success_list = [
            (msa_name, _load_nexus_alignment(msa_name, nexus_path))
            for msa_name, nexus_path in self.success_list
        ]

        # Step 1. Define output names
        out_fn_nucl_concat_fasta = os.path.join(
            self.user_params.get("out_dir"), "nucl_" + str(len(self.success_list)) + "concat.aligned.fasta"
        )
        out_fn_nucl_concat_nexus = os.path.join(
            self.user_params.get("out_dir"), "nucl_" + str(len(self.success_list)) + "concat.aligned.nexus"
        )

        # Step 2. Do concatenation
        try:
            alignm_concat = Nexus.Nexus.combine(
                self.success_list
            )  
        except Exception as e:
            log.critical("Unable to concatenate alignments.\n" f"Error message: {e}")
            raise Exception()

        # Step 3. Write concatenated alignment to file,
        log.info(" writing concatenation to file in NEXUS and FASTA formats")
        _concat_sync()
       
# -----------------------------------------------------------------#


class BackTranslation:
    def __init__(
            self,
            prot_name: str,
            prot_align_file: str,
            nuc_fasta_file: str,
            nuc_align_file: str,
            table: int = 0,
            align_format: str = "fasta",
            gap: str = "-"
    ):
        """
        Coordinate the back-translation of protein sequences to nucleotide sequences.

        Args:
            prot_name: Name of the protein being back-translated
            prot_align_file: Path to the file containing the aligned protein sequences.
            nuc_fasta_file: Path to the file containing the unaligned nucleotide sequences.
            nuc_align_file: Path to the output file for the back-translated nucleotide sequences.
            table: Genetic code table number (default is 0).
            align_format: Format of the alignment file (default is 'fasta').
            gap: String that designates a nucleotide gap (default is '-').
        """
        self.prot_name = prot_name
        self.align_format = align_format
        self.prot_align_file = prot_align_file
        self.nuc_fasta_file = nuc_fasta_file
        self.nuc_align_file = nuc_align_file
        self.table = table
        self.gap = gap

    def _evaluate_nuc(self, identifier: str, nuc: Seq, prot: Seq) -> Seq:
        """
        Returns nucleotide sequence if works (can remove trailing stop).
        """
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
                    f"translation for {identifier} would match if {nuc[0:3].upper()} "
                    f"was a start codon (check correct table used)"
                )
                log.warning(f"translation check failed for {identifier}")

        else:
            m = "".join("." if x == y else "!" for (x, y) in zip(p, t))
            if len(prot) < 70:
                log.debug(
                    f"translation mismatch for {identifier} [0:{len(prot)}]\n"
                    f"protein:     {p}\n"
                    f"             {m}\n"
                    f"translation: {t}\n"
                )
            else:
                for offset in range(0, len(p), 60):
                    log.debug(
                        f"translation mismatch for {identifier} [{offset}:{offset + 60}]\n"
                        f"protein:     {p[offset:offset + 60]}\n"
                        f"             {m[offset:offset + 60]}\n"
                        f"translation: {t[offset:offset + 60]}\n"
                    )
            log.warning(f"translation check failed for {identifier}")

    def _backtrans_seq(
            self,
            identifier: str,
            aligned_protein_record: SeqRecord,
            unaligned_nucleotide_record: SeqRecord,
    ) -> SeqRecord:
        # Per https://biopython.org/docs/1.81/api/Bio.Seq.html,
        # `replace` is proper replacement for depreciated method `ungap`.
        try:
            ungapped_protein = aligned_protein_record.seq.replace(self.gap, "")
        except AttributeError:
            ungapped_protein = aligned_protein_record.seq.ungap(self.gap)

        ungapped_nucleotide = unaligned_nucleotide_record.seq
        if self.table:
            ungapped_nucleotide = self._evaluate_nuc(
                identifier, ungapped_nucleotide, ungapped_protein
            )
        elif len(ungapped_protein) * 3 != len(ungapped_nucleotide):
            log.debug(
                f"Inconsistent lengths for {aligned_protein_record.id}, "
                f"ungapped protein {len(ungapped_protein)}, "
                f"tripled {len(ungapped_protein) * 3} vs "
                f"ungapped nucleotide {len(ungapped_nucleotide)}"
            )
            log.warning(
                f"Backtranslation failed for {identifier} due to ungapped length mismatch"
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
                f"Backtranslation failed for {identifier} due to unaligned length mismatch"
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
                f"Backtranslation failed for {identifier} due to aligned length mismatch"
            )
            return None

        return aligned_nuc

    def _backtrans_seqs(self, protein_alignment: MultipleSeqAlignment, nucleotide_records: Mapping,
                        key_function: Callable[[str], str] = lambda x: x):
        """
        Thread nucleotide sequences onto a protein alignment.
        """
        aligned = []
        for protein in protein_alignment:
            identifier = self._get_identifier(protein.id)
            try:
                nucleotide = nucleotide_records[key_function(protein.id)]
            except KeyError:
                raise ValueError(
                    f"Could not find nucleotide sequence for {identifier}"
                )
            sequence = self._backtrans_seq(identifier, protein, nucleotide)
            if sequence is not None:
                aligned.append(sequence)
        return MultipleSeqAlignment(aligned)

    def _get_identifier(self, rec_id: str) -> str:
        rec_name = rec_id.replace(self.prot_name + "_", "")
        identifier = f"protein '{self.prot_name}' in {rec_name}"
        return identifier

    def backtranslate(self):
        """
        Perform back-translation of a protein alignment to nucleotides.
        """
        # Step 1. Load the protein alignment
        prot_align = AlignIO.read(self.prot_align_file, self.align_format)
        # Step 2. Index the unaligned nucleotide sequences
        nuc_dict = SeqIO.index(self.nuc_fasta_file, "fasta")
        # Step 3. Perform back-translation
        nuc_align = self._backtrans_seqs(prot_align, nuc_dict)
        # Step 4. Write the back-translated nucleotide alignment to a file
        with open(self.nuc_align_file, "w") as output_handle:
            AlignIO.write(nuc_align, output_handle, self.align_format)


class MAFFT:

    def __init__(self, user_params: Optional[Dict[str, Any]] = None):
        """
        An interface for executing MAFFT alignment.

        If the instance is able to detect an installed version of MAFFT, that will be used for alignment.
        Alternatively, downloading a version to be used by the package will be attempted.
        If this download is successful, this version will be used by the current instance of `MAFFT`
        and will be immediately available for subsequent initializations.

        Args:
            user_params: Specifications for how the alignment should be performed.
                This is `num_threads`.
        """
        log.info("configuring alignment tool")
        self._set_num_threads(user_params)
        self._set_exec_path()
        log.info(
            f" conducting alignments using mafft with "
            f"{self.num_threads} threads"
        )

    def _set_num_threads(self, user_params: Optional[Dict[str, Any]]):
        if user_params and isinstance(user_params, dict) and user_params.get("num_threads"):
            # self.num_threads = str(user_params.get("num_threads"))
            self.num_threads = "1"
        else:
            self.num_threads = "1"


    def _set_exec_path(self):
        self._check_mafft()


    def _check_mafft(self):
        if shutil.which("mafft"):
            log.info(f"  using installed MAFFT for alignment")
            self.exec_path = "mafft"
        else:
            log.info(f"  unable to find MAFFT installed")

   

    def align(self, input_file: str, output_file: str):
        """
        Performs sequence alignment using MAFFT.

        Args:
            input_file: Path of input sequence file.
            output_file: Path of output sequence file.
        """
        if not self.exec_path:
            log.error("  skipping attempted alignment; alignment tool is not properly configured")
        elif not os.path.exists(input_file) or not isinstance(input_file, str):
            log.error("  skipping attempted alignment; specified input file does not exist")
        elif not isinstance(output_file, str):
            log.error(f"  skipping attempted alignment; specified output file is incorrect type: {type(output_file)}")
        else:
            mafft_cmd = [self.exec_path, "--thread", self.num_threads, "--adjustdirection", input_file]
            with open(output_file, 'w') as hndl, open(os.devnull, 'w') as devnull:
                process = subprocess.Popen(mafft_cmd, stdout=hndl, stderr=devnull, text=True)
                process.wait()
