import os
from multiprocessing import Pool
from itertools import repeat
from argparse import ArgumentParser
from pyfaidx import Fasta

from semvar.constants import CHR_REFSEQ_DICT
from semvar.io import load_sems
from semvar.utils import ref_mismatch, contains_invalid_chars, no_valid_kmers, get_spdi
from semvar.seq import get_sequences, get_seq_context
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

def run_annotation(tf_name, sem_dict, variants_file, output_dir, assembly, only_report_effects):
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
    
    with open(variants_file) as f, open(os.path.join(output_dir, f'{tf_name}_SEM_predictions.tsv'), 'w+') as output:
        output.writelines([
                            f"#Description: Predictions of variant effects on transcription factor binding\n",
                            f"#TFName: {tf_name}\n",
                            f"#BiosampleOntologyTermName: N/A\n",
                            f"#BiosampleOntologyTermID: N/A\n",
                            f"#AssayContext: ChIP-seq\n",
                            f"#Model: SEMpl_v1.0.0\n"
                            ])
        colnames = '\t'.join(['chr', 'pos', 'spdi', 'ref', 'alt', 'ref_seq_context', 'alt_seq_context', 'ref_score', 'alt_score', 'variant_effect_score', 'pvalue', 'SEMpl.annotation', 'SEMpl.baseline', 'SEMpl.SEM_filename'])
        output.write(f'{colnames}\n')

        for line in f:
            variant_info = line.strip().split('\t')
            chrom = variant_info[0]
            pos = int(variant_info[1])
            var_id = variant_info[2]
            ref = variant_info[3]
            alt = variant_info[4]
            if contains_invalid_chars(ref) or contains_invalid_chars(alt):
                print(f'Warning: Invalid characters for ref and/or alt alleles. Skipping this variant. Check variant {chrom}:{pos}:{ref}:{alt}')
                continue

            if ref_mismatch(chrom, pos, ref, fasta):
                print(f'Warning: Reference allele mismatch at {chrom}:{pos}:{ref}:{alt}. Skipping this variant.')
                continue

            sem_len = len(mat)
            ref_seq, alt_seq = get_sequences(sem_len, fasta, chrom, pos, ref, alt)
            if no_valid_kmers(ref_seq, sem_len) or no_valid_kmers(alt_seq, sem_len):
                print(f'Warning: No valid kmers. Skipping this variant. Check variant {chrom}:{pos}:{ref}:{alt}')
                continue

            annot = 0
            annot_score = 0
            
            ref_score, ref_score_idx = get_best_sem_score(ref_seq, mat)
            alt_score, alt_score_idx = get_best_sem_score(alt_seq, mat)
            annot, annot_score = annotate_variant(ref_score, alt_score, baseline)
            ref_score = 2**ref_score # unlog
            alt_score = 2**alt_score
            ref_seq_context = get_seq_context(ref_score_idx, chrom, pos, sem_len)
            alt_seq_context = get_seq_context(alt_score_idx, chrom, pos, sem_len)
            
            if (only_report_effects == True) and (annot in ["no_binding", "binding_unchanged"]):
                continue
            
            # Write to output file
            # spdi = get_spdi(CHR_REFSEQ_DICT, chrom, pos-1, ref, alt) SPDI given in input VCF
            variant_output = '\t'.join([chrom, str(pos), var_id, ref, alt, ref_seq_context, alt_seq_context, str(ref_score), str(alt_score), str(annot_score), 'N/A', annot, str(baseline), sem_filename])     
            output.write(f'{variant_output}\n')

if __name__ == "__main__":

    parser = ArgumentParser(prog='SEM_variant_annotator')
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
        