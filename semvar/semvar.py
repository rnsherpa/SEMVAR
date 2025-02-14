import os
from multiprocessing import Pool
from itertools import repeat
from argparse import ArgumentParser
from pyfaidx import Fasta

from semvar.constants import CHR_REFSEQ_DICT
from semvar.io import load_sems, load_baselines, get_spdi
from semvar.utils import ref_mismatch, contains_invalid_chars, no_valid_kmers
from semvar.seq import get_sequences
from semvar.scoring import get_best_sem_score

def annotate_variant(ref_score, alt_score, baseline):
    '''Compare the SEM sum of a variant to the baseline of that SEM to annotate
    Params
        ref_score (float): SEM sum score of reference allele
        alt_score (float): SEM sum score of alternate allele
        baseline (float): SEM sum score of random background for a given motif
    Returns
        annot (str): Annotation label of variant effect on motif
        annot_score (float): Difference of SEM scores between alt and ref. Since SEM scores are log2 transformed, this can be interpreted as a log2 fold change of binding affinity due to the variant.
    '''
    if alt_score < baseline < ref_score: 
        annot = "binding_ablated" 
    elif baseline < alt_score < ref_score:
        annot = "binding_decreased"
    elif alt_score > baseline > ref_score:
        annot = "binding_created"
    elif alt_score > ref_score > baseline:
        annot = "binding_increased"
    elif alt_score == ref_score:
        annot = "binding_unchanged"
    elif (ref_score < baseline) and (alt_score < baseline):
        annot = "no_binding"
    else: 
        raise ValueError(f"{alt_score},{ref_score},{baseline} Scores do not fit to any of the annotation labels, check score calculation/assignment")
    
    annot_score = alt_score-ref_score

    return annot, annot_score

def run_annotation(tf_name, sem_dict, variants_file, output_dir, baselines, assembly, only_report_effects):
    '''Create annotation file of variants for a given SEM
    Params
        tf_name (str): Name of TF the SEM is modeling
        sem_dict (dict): Dict of the SEM data (matrices, baselines, filenames) corresponding to tf_name
        variants_file (str): Path to variants file
        output_dir (str): Path to output directory
        assembly (str): Path to genome assembly fasta (must have corrsponding .fai index file in the same directory)
        only_report_effects (bool): If true, exclude variants with "no_binding" or "binding_unchanged" from final output
    '''         
    mat = sem_dict['mat']
    baseline = sem_dict['baseline']
    sem_filename = sem_dict['filename']
    fasta = Fasta(assembly)
    
    with open(variants_file) as f, open(os.path.join(output_dir,f'{tf_name}_SEM_predictions.tsv'), 'w+') as output:
        output.writelines([f'#SEM_file={sem_filename}\n',
                           f'#TF={tf_name}\n'])
        colnames = '\t'.join(['chrom', 'start', 'end', 'spdi', 'ref', 'alt', 'kmer_coord', 'ref_score', 'alt_score', 'log2fc', 'effect_on_binding'])
        output.write(f'{colnames}\n')

        for line in f:
            variant_info = line.strip().split('\t')
            chrom = variant_info[0]
            start = int(variant_info[1])
            end = int(variant_info[2])
            ref = variant_info[3]
            alt = variant_info[4]
            if contains_invalid_chars(ref) or contains_invalid_chars(alt):
                print(f'Warning: Invalid characters for ref and/or alt alleles. Skipping this variant. Check variant {chrom}:{end}:{ref}:{alt}')
                continue

            if ref_mismatch(chrom, end, ref, assembly):
                continue

            sem_len = len(mat)
            ref_seq, alt_seq = get_sequences(sem_len, fasta, chrom, end, alt)
            if no_valid_kmers(ref_seq, sem_len) or no_valid_kmers(alt_seq, sem_len):
                print(f'Warning: No valid kmers. Skipping this variant. Check variant {chrom}:{end}:{ref}:{alt}')
                continue

            annot = 0
            annot_score = 0
            
            ref_score = get_best_sem_score(ref_seq, mat)
            alt_score = get_best_sem_score(alt_seq, mat)
            annot, annot_score = annotate_variant(ref_score, alt_score, baseline)
            ref_score = 2**ref_score # unlog
            alt_score = 2**alt_score
            
            if (only_report_effects == True) and (annot in ["no_binding", "binding_unchanged"]):
                continue
            
            # Write to output file
            spdi = get_spdi(CHR_REFSEQ_DICT, chrom, start, ref, alt)
            variant_output = '\t'.join([chrom, str(start), str(end), spdi, ref, alt, ref_score, alt_score, annot_score, annot])     
            output.write(f'{variant_output}\n')

if __name__ == "__main__":

    parser = ArgumentParser(prog='SEM_variant_annotation')
    parser.add_argument('--file', '-f', help='Input variant list', required=True)
    parser.add_argument('--semsdir', '-d', help='Path to directory containing SEMs', required=True)
    parser.add_argument('--sems', '-s', nargs='*', default=None, help='List of sems. If not specified, all SEMs from the semsdir will be used')
    parser.add_argument('--assembly', '-a', help='Path to indexed reference genome', required=True)
    parser.add_argument('--n_processes', '-n', type=int, default=1, help='Number of processes for multiprocessing')
    parser.add_argument('--only-report-effects', '-e', action='store_true', help='Only include variants with annotated effect in final output')
    parser.add_argument('--outdir', '-o', help='Output path of variant annotations')
    args = parser.parse_args()

    variants_file = args.file
    sems_dir = args.semsdir
    sems_list = args.sems
    assembly = args.assembly
    n_processes = args.n_processes
    output_dir = args.outdir

    sems_dicts = load_sems(sems_dir, sems_list)
    tf_name_list = list(sems_dicts.keys())
    sems_dict_list = [sems_dicts[tf_name] for tf_name in tf_name_list]

    # Run annotation with multiprocessing
    with Pool(n_processes) as pool:
        pool.starmap(run_annotation, zip(tf_name_list, sems_dict_list, repeat(variants_file), repeat(output_dir), repeat(assembly), repeat(args.only_report_effects)))
        