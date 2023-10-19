 #!/usr/bin/env python3
'''Extract and align coding and non-coding regions across multiple plastid genomes'''
__version__ = 'm_gruenstaeudl@fhsu.edu|Mon 09 Oct 2023 08:18:54 PM CDT'

# ------------------------------------------------------------------------------#
## IMPORTS
import argparse
from typing import Sequence
from Bio import SeqIO, AlignIO, Nexus, SeqRecord  # line is necessary for similar reason stated in: https://www.biostars.org/p/13099/
from Bio.Align import Applications  # line is necessary for similar reason stated in: https://www.biostars.org/p/13099/
from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition
import coloredlogs
import collections
import copy
import io
import logging
import os
import re
import subprocess
import sys
from Bio.Align.Applications import MafftCommandline

# ------------------------------------------------------------------------------#
## DEBUGGING HELP
import ipdb


# ipdb.set_trace()
# -----------------------------------------------------------------#
# CLASSES AND FUNCTIONS
class ExtractAndCollect:
    def __init__(self, main_d_nucl):
        self.main_d_nucl = main_d_nucl

    def do_CDS(self, main_d_prot, rec, len_cutoff, log):
        for feature in rec.features:
            if feature.type == 'CDS' and 'gene' in feature.qualifiers:
                gene_name = feature.qualifiers['gene'][0]
                seq_name = f"{gene_name}_{rec.name}"

                # Nucleotide sequences
                seq_obj = feature.extract(rec).seq
                seq_rec = SeqRecord.SeqRecord(seq_obj, id=seq_name, name="", description="")

                self.update_main_d_nucl(gene_name, seq_rec)

                # Protein sequences
                seq_obj = feature.extract(rec).seq.translate(table=11)
                seq_rec = SeqRecord.SeqRecord(seq_obj, id=seq_name, name="", description="")

                if len(seq_obj) >= len_cutoff:
                    self.update_main_d_prot(main_d_prot, gene_name, seq_rec)

    def do_IGS(self, rec, _fname, len_cutoff, log):
        all_genes_minus_matK = [feature for feature in rec.features if (
            feature.type == 'gene' and
            'gene' in feature.qualifiers and
            feature.qualifiers['gene'][0] != 'matK'
        )]

        for idx in range(len(all_genes_minus_matK) - 1):
            cur_feat = all_genes_minus_matK[idx]
            adj_feat = all_genes_minus_matK[idx + 1]

            if 'gene' in cur_feat.qualifiers and 'gene' in adj_feat.qualifiers:
                cur_feat_name = cur_feat.qualifiers['gene'][0]
                adj_feat_name = adj_feat.qualifiers['gene'][0]

                IGS_name = self.create_IGS_name(cur_feat_name, adj_feat_name, _fname)

                if type(cur_feat.location) is not CompoundLocation and type(adj_feat.location) is not CompoundLocation:
                    exact_location = FeatureLocation(cur_feat.location.end, adj_feat.location.start)

                    if int(exact_location.start) < int(exact_location.end):
                        seq_obj = exact_location.extract(rec).seq
                        seq_name = f"{IGS_name}_{rec.name}"
                        seq_rec = SeqRecord.SeqRecord(seq_obj, id=seq_name, name="", description="")

                        if len(seq_obj) >= len_cutoff:
                            self.update_main_d_nucl(IGS_name, seq_rec)

    def do_INT(self, main_d_intron2, rec, log):
        for feature in rec.features:
            if feature.type in ['CDS', 'tRNA']:
                gene_name_base = feature.qualifiers.get('gene', ['Unknown'])[0]
                gene_name_base_SAFE = self.sanitize_gene_name(gene_name_base)

                if len(feature.location.parts) == 2:
                    gene_name = f"{gene_name_base_SAFE}_intron1"
                    seq_rec, gene_name = self.extract_INT_internal(rec, feature, gene_name, 0, log)
                    self.update_main_d_nucl(gene_name, seq_rec)

                if len(feature.location.parts) == 3:
                    copy_feature = feature
                    gene_name = f"{gene_name_base_SAFE}_intron1"
                    seq_rec, gene_name = self.extract_INT_internal(rec, feature, gene_name, 0, log)
                    self.update_main_d_nucl(gene_name, seq_rec)

                    feature = copy_feature
                    gene_name = f"{gene_name_base_SAFE}_intron2"
                    seq_rec, gene_name = self.extract_INT_internal(rec, feature, gene_name, 1, log)
                    self.update_main_d_intron2(gene_name, seq_rec)

    def extract_INT_internal(self, rec, feature, gene_name, offset, log):
        try:
            feature.location = FeatureLocation(
                feature.location.parts[offset].end,
                feature.location.parts[offset + 1].start
            )
        except:
            feature.location = FeatureLocation(
                feature.location.parts[offset + 1].start,
                feature.location.parts[offset].end
            )

        seq_name = f"{gene_name}_{rec.name}"
        seq_obj = feature.extract(rec).seq
        seq_rec = SeqRecord.SeqRecord(seq_obj, id=seq_name, name="", description="")

        return seq_rec, gene_name

    def update_main_d_nucl(self, gene_name, seq_rec):
        if gene_name in self.main_d_nucl:
            self.main_d_nucl[gene_name].append(seq_rec)
        else:
            self.main_d_nucl[gene_name] = [seq_rec]

    def update_main_d_prot(self,main_d_prot, gene_name, seq_rec):
        if gene_name in main_d_prot:
            main_d_prot[gene_name].append(seq_rec)
        else:
            main_d_prot[gene_name] = [seq_rec]

    def create_IGS_name(self, cur_feat_name, adj_feat_name, _fname):
        cur_feat_name_SAFE = self.sanitize_gene_name(cur_feat_name)
        adj_feat_name_SAFE = self.sanitize_gene_name(adj_feat_name)
        IGS_name = f"{cur_feat_name_SAFE}_{adj_feat_name_SAFE}"
        inv_IGS_name = f"{adj_feat_name_SAFE}_{cur_feat_name_SAFE}"
        return IGS_name

    def sanitize_gene_name(self, gene_name):
        gene_name_SAFE = gene_name.replace('-', '_')
        gene_name_SAFE = re.sub(r'[^\w]', '', gene_name_SAFE)
        return gene_name_SAFE


# -----------------------------------------------------------------#

def proteinalign_and_backtranslate(main_d_prot, out_dir, log):
    # Step 1. Write unaligned protein sequences to file
    for k, v in main_d_prot.items():
        out_fn_unalign_prot = os.path.join(out_dir, f'prot_{k}.unalign.fasta')
        out_fn_aligned_prot = os.path.join(out_dir, f'prot_{k}.aligned.fasta')
        
        with open(out_fn_unalign_prot, 'w') as hndl:
            SeqIO.write(v, hndl, 'fasta')

        # Step 2. Align sequences with MAFFT
        # import subprocess
        # subprocess.call(['mafft', '--auto', out_fn_unalign_prot, '>', out_fn_aligned_prot])
        mafft_align(out_fn_unalign_prot, out_fn_aligned_prot)

    # Step 4. Check if back-translation script exists
    path_to_back_transl_helper = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'align_back_trans.py')
    if not os.path.isfile(path_to_back_transl_helper):
        log.critical("Unable to find `align_back_trans.py` alongside this script")
        raise Exception()

    # Step 5. Conduct actual back-translation
    for k, v in main_d_prot.items():
        out_fn_unalign_nucl = os.path.join(out_dir, f'nucl_{k}.unalign.fasta')
        out_fn_aligned_nucl = os.path.join(out_dir, f'nucl_{k}.aligned.fasta')
        out_fn_aligned_prot = os.path.join(out_dir, f'prot_{k}.aligned.fasta')

        ## TO DO - BUG! ##
        # For some reason, the path_to_back_transl_helper spits out only nexus files
        cmd = ['python3', path_to_back_transl_helper, 'fasta', out_fn_aligned_prot, out_fn_unalign_nucl,
               out_fn_aligned_nucl, '11']

        try:
            ## TO DO - BUG! ##
            # For some reason, stderr is not properly saved to log_msg in cases when the back-translation python package does not work; then log_msg is empty

            log_msg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            cmd_prt = ' '.join(cmd)
            ## TO DO - BUG! ##
            # The following lines do not work because log_msg is empty, as indicated above.

            # log_prt = str(log_msg, 'utf-8')
            # log.warning("Unable to conduct back-translation of `%s`. Command used: %s. Error message: %s" % (k, cmd_prt, log_prt))
           
            log.warning(f"Unable to conduct back-translation of `{k}`. Command used: {cmd_prt}. Error message: {e.output.decode('utf-8')}")

def mafft_align(input_file, output_file):
    # Perform sequence alignment using MAFFT
    mafft_cline = Applications.MafftCommandline(input=input_file)
    stdout, stderr = mafft_cline()

    with open(output_file, 'w') as hndl:
        hndl.write(stdout)

# ------------------------------------------------------------------------------#

def extract_INT_internal(rec, feature, gene_name, offset, log):
    # Extract intron region
    try:
        # Calculate the location of the intron
        intron_start = feature.location.parts[offset].end
        intron_end = feature.location.parts[offset + 1].start
        feature.location = FeatureLocation(intron_start, intron_end)
    except:
        # Swap the positions if unable to calculate the location
        intron_start = feature.location.parts[offset + 1].start
        intron_end = feature.location.parts[offset].end
        feature.location = FeatureLocation(intron_start, intron_end)

    # Extract the intron sequence
    try:
        seq_name = f"{gene_name}_{rec.name}"
        seq_obj = feature.extract(rec).seq
        seq_rec = SeqRecord.SeqRecord(seq_obj, id=seq_name, name="", description="")
        return seq_rec, gene_name
    except:
        log.critical(f"Unable to conduct intron extraction for {feature.qualifiers['gene']}")
        raise Exception()

# ------------------------------------------------------------------------------#
def remove_duplicates(my_dict):
    for k, v in my_dict.items():
        idtags = []
        for counter, seqrec in enumerate(v):
            if seqrec.id in idtags:
                # v.pop(counter)
                del v[counter]
            else:
                idtags.append(seqrec.id)
        my_dict[k] = v
    # return my_dict

# ------------------------------------------------------------------------------#
## MAIN HELPER FUNCTIONS
# ------------------------------------------------------------------------------#
def setup_logger(verbose):
    global log
    log = logging.getLogger(__name__)
    log_format = '%(asctime)s [%(levelname)s] %(message)s'
    log_level = logging.DEBUG if verbose else logging.INFO
    coloredlogs.install(fmt=log_format, level=log_level, logger=log)
    return log

def unpack_input_parameters(args):
    in_dir = args.inpd
    out_dir = args.outd
    fileext = args.fileext
    exclude_list = args.excllist
    len_cutoff = args.lencutoff
    tax_cutoff = args.taxcutoff
    select_mode = args.selectmode.lower()
    verbose = args.verbose
    return in_dir, out_dir, fileext, exclude_list, len_cutoff, tax_cutoff, select_mode, verbose

def extract_and_collect_annotations(in_dir, fileext, select_mode, len_cutoff, log):
    main_d_nucl = collections.OrderedDict()
    main_d_prot = collections.OrderedDict() if select_mode == 'cds' else None
    main_d_intron2 = collections.OrderedDict() if select_mode == 'int' else None

    action = "extracting annotations from genome records"
    log.info("%s" % action)

    files = [f for f in os.listdir(in_dir) if f.endswith(fileext)]
    for f in files:
        log.info('reading file `%s`' % (f))
        rec = SeqIO.read(os.path.join(in_dir, f), 'genbank')

        if select_mode == 'cds':
            ExtractAndCollect(main_d_nucl).do_CDS(main_d_prot, rec, len_cutoff, log)
        if select_mode == 'igs':
            ExtractAndCollect(main_d_nucl).do_IGS(rec, f, len_cutoff, log)
        if select_mode == 'int':
            ExtractAndCollect(main_d_nucl).do_INT(main_d_intron2, rec, log)
            main_d_nucl.update(main_d_intron2)

        if not main_d_nucl.items():
            log.critical("No items in main dictionary")
            raise Exception()

    return main_d_nucl, main_d_prot

def remove_annotations_below_taxa_cutoff(main_d_nucl, main_d_prot, tax_cutoff):
    for k, v in main_d_nucl.items():
        if len(v) < tax_cutoff:
            del main_d_nucl[k]
            if main_d_prot:
                del main_d_prot[k]

def remove_orfs(main_d_nucl, main_d_prot):
    list_of_orfs = [orf for orf in main_d_nucl.keys() if "orf" in orf]
    for orf in list_of_orfs:
        del main_d_nucl[orf]
        if main_d_prot:
            del main_d_prot[orf]

def remove_user_defined_genes(main_d_nucl, main_d_prot, exclude_list, select_mode, log):
    if exclude_list:
        if select_mode == 'int':
            to_be_excluded = [i + "_intron1" for i in exclude_list] + [i + "_intron2" for i in exclude_list]
            exclude_list = to_be_excluded
        for excluded in exclude_list:
            if excluded in main_d_nucl:
                del main_d_nucl[excluded]
                if select_mode == 'cds' and main_d_prot:
                    del main_d_prot[excluded]
            else:
                log.warning("Region `%s` to be excluded but unable to be found in infile." % excluded)
                pass

def extract_and_combine_nucleotide_sequences(main_d_nucl, out_dir):
    action = "extracting and combining nucleotide sequences genewise"
    log.info("%s" % action)

    for k, v in main_d_nucl.items():
        out_fn_unalign_nucl = os.path.join(out_dir, 'nucl_' + k + '.unalign.fasta')
        with open(out_fn_unalign_nucl, 'w') as hndl:
            SeqIO.write(v, hndl, 'fasta')

def multiple_sequence_alignment_nucleotide(main_d_nucl, out_dir):
    action = "conducting multiple sequence alignment based on nucleotide sequence data"
    log.info("%s" % action)

    if main_d_nucl.items():
        for k, v in main_d_nucl.items():
            out_fn_unalign_nucl = os.path.join(out_dir, 'nucl_' + k + '.unalign.fasta')
            out_fn_aligned_nucl = os.path.join(out_dir, 'nucl_' + k + '.aligned.fasta')
            num_threads = 1
            alignm_cmdlexec = Applications.MafftCommandline(input=out_fn_unalign_nucl, adjustdirection=True, thread=num_threads)
            stdout, stderr = alignm_cmdlexec()
            with open(out_fn_aligned_nucl, 'w') as hndl:
                hndl.write(stdout)
    else:
        log.critical("No items in nucleotide main dictionary to process")
        raise Exception()

def conduct_protein_alignment_and_back_translation(main_d_prot, out_dir, log):
    action = "conducting multiple sequence alignment based on protein sequence data, followed by back-translation to nucleotides"
    log.info("%s" % action)
    try:
        proteinalign_and_backtranslate(main_d_prot, out_dir, log)
    except:
        raise Exception()

def convert_fasta_alignment_to_nexus(main_d_nucl, out_dir):
    action = "converting FASTA alignment to NEXUS alignment"
    log.info("%s" % action)
    successList = []

    for k in main_d_nucl.keys():
        aligned_nucl_fasta = os.path.join(out_dir, 'nucl_' + k + '.aligned.fasta')
        aligned_nucl_nexus = os.path.join(out_dir, 'nucl_' + k + '.aligned.nexus')
        try:
            AlignIO.convert(aligned_nucl_fasta, 'fasta', aligned_nucl_nexus, 'nexus', molecule_type='DNA')
        except:
            log.warning("Unable to convert alignment of `%s` from FASTA to NEXUS" % k)
            continue
        try:
            alignm_nexus = AlignIO.read(aligned_nucl_nexus, 'nexus')
            hndl = io.StringIO()
            AlignIO.write(alignm_nexus, hndl, 'nexus')
            nexus_string = hndl.getvalue()
            nexus_string = nexus_string.replace('\n' + k + '_', '\nconcatenated_')
            alignm_nexus = Nexus.Nexus.Nexus(nexus_string)
            successList.append((k, alignm_nexus))
        except:
            log.warning("Unable to add alignment of `%s` to concatenation" % k)
            pass

    return successList

def concatenate_alignments(successList, out_dir):
    action = "concatenate alignments (in no particular order)"
    log.info("%s" % action)
    out_fn_nucl_concatenated_fasta = os.path.join(out_dir, 'nucl_' + str(len(successList)) + 'concatenated.aligned.fasta')
    out_fn_nucl_concatenated_nexus = os.path.join(out_dir, 'nucl_' + str(len(successList)) + 'concatenated.aligned.nexus')

    try:
        alignm_concatenated = Nexus.Nexus.combine(successList)
    except:
        log.critical("Unable to concatenate alignments")
        raise Exception()
    
    alignm_concatenated.write_nexus_data(filename=open(out_fn_nucl_concatenated_nexus, 'w'))
    AlignIO.convert(out_fn_nucl_concatenated_nexus, 'nexus', out_fn_nucl_concatenated_fasta, 'fasta')

# ------------------------------------------------------------------------------#
## MAIN
# ------------------------------------------------------------------------------#
def main(args):
    in_dir, out_dir, fileext, exclude_list, len_cutoff, tax_cutoff, select_mode, verbose = unpack_input_parameters(args)
    log = setup_logger(verbose)
    
    main_d_nucl, main_d_prot = extract_and_collect_annotations(in_dir, fileext, select_mode, len_cutoff, log)

    remove_duplicates(main_d_nucl)
    remove_duplicates(main_d_nucl)
    
    if select_mode == 'cds':
        remove_duplicates(main_d_prot)
        remove_duplicates(main_d_prot)

    remove_annotations_below_taxa_cutoff(main_d_nucl, main_d_prot, tax_cutoff)
    remove_orfs(main_d_nucl, main_d_prot)
    remove_user_defined_genes(main_d_nucl, main_d_prot, exclude_list, select_mode, log)
    
    extract_and_combine_nucleotide_sequences(main_d_nucl, out_dir)
    
    if not select_mode == 'cds':
        multiple_sequence_alignment_nucleotide(main_d_nucl, out_dir)
    
    if select_mode == 'cds':
        conduct_protein_alignment_and_back_translation(main_d_prot, out_dir, log)
    
    successList = convert_fasta_alignment_to_nexus(main_d_nucl, out_dir)
    concatenate_alignments(successList, out_dir)
    
    action = "end of script\n"
    log.info("%s" % action)
    quit()


# ------------------------------------------------------------------------------#
# ARGPARSE
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Author|Version: ' + __version__)
    # Required
    parser.add_argument("--inpd", "-i", type=str, required=True,
                        help="path to input directory (which contains the GenBank files)",
                        default="./input")
    # Optional
    parser.add_argument("--outd", "-o", type=str, required=False,
                        help="(Optional) Path to output directory",
                        default="./output")
    parser.add_argument("--selectmode", "-s", type=str, required=False,
                        help="(Optional) Type of regions to be extracted (i.e. `cds`, `int`, or `igs`)",
                        default="cds")
    parser.add_argument("--fileext", "-f", type=str, required=False,
                        help="(Optional) File extension of input files",
                        default=".gb")
    parser.add_argument("--excllist", "-e", type=list, required=False,
                        default=['rps12'],
                        help="(Optional) List of genes to be excluded")
    parser.add_argument("--lencutoff", "-l", type=int, required=False,
                        help="(Optional) Sequence length below which regions will not be extracted",
                        default=1)
    parser.add_argument("--taxcutoff", "-t", type=int, required=False,
                        help="(Optional) Minimum number of taxa in which a region must be present to be extracted",
                        default=1)
    parser.add_argument("--verbose", "-v", action="version", version="%(prog)s " + __version__,
                        help="(Optional) Enable verbose logging", default=True)
    args = parser.parse_args()
    main(args)
# ------------------------------------------------------------------------------#
# EOF
# ------------------------------------------------------------------------------#
