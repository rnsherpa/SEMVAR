import os
from multiprocessing import Pool
from itertools import repeat
from argparse import ArgumentParser
from pyfaidx import Fasta

from constants import CHR_REFSEQ_DICT
from io import load_sems, load_baselines, get_spdi
from utils import ref_mismatch
from scoring import get_best_scores

def get_kmers(sem_len, fasta, chr, pos, alt):
    '''Get all possible reference and variant k-mers overlapping the variant position where k is the length of the SEM
    Params
        sem_len (int): length of SEM
        fasta (pyfaidx.Fasta): Pyfaidx Fasta object of the reference genome assembly
        chr (str): Chromosome of SNV in 'chrN' format
        pos (int): 1-based position of SNV
        alt (str): Alternate allele
    Returns
        kmers (dict of lists): Lists of k-mers and their coordinates
    '''
    kmers = {
        'ref': [],
        'alt': [],
        'coords': [],
    }

    for idx, kmer_start in enumerate(range(pos-sem_len, pos)): # begin at kmer where variant is at end of motif, last kmer has variant at beginning of motif
        if kmer_start>=0: # positions cannot be negative values
            kmer_end = kmer_start+sem_len
            ref_kmer = fasta[chr][kmer_start:kmer_end].seq.upper()
            var_pos_kmer = (sem_len - 1) - idx # position in kmer where variant is located
            alt_kmer = ref_kmer[:var_pos_kmer] + alt + ref_kmer[var_pos_kmer+1:]
            if 'N' in ref_kmer: # don't include kmers that contain N
                continue
            else:
                kmers['ref'].append(ref_kmer) 
                kmers['alt'].append(alt_kmer)
                kmers['coords'].append(f'{chr}:{kmer_start}-{kmer_end}')

    return kmers

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
        annot = ""
        raise ValueError(f"{alt_score},{ref_score},{baseline} Scores do not fit to any of the annotation labels, check score calculation/assignment")
    
    annot_score = alt_score-ref_score

    return annot, annot_score

def run_annotation(sem, sem_filename, variants_file, output_dir, baselines, assembly, only_report_effects):
    '''Create annotation file of variants for a given SEM
    Params
        sem (str): Name of TF the SEM was created for
        sem_filename (str): Filename of SEM used to annotate variants
        variants_file (str): Path to variants file
        output_dir (str): Path to output directory
        baselines (dict): Keys are motif names and values are baseline values for the corresponding SEM
        assembly (str): Path to genome assembly fasta (must have corrsponding .fai index file in the same directory)
        only_report_effects (bool): If true, exclude variants with "no_binding" or "binding_unchanged" from final output
    '''         
    mat = sems[sem]
    baseline = baselines[sem]
    fasta = Fasta(assembly)
    
    with open(variants_file) as f, open(os.path.join(output_dir,f'{sem}_annotations.tsv'), 'w+') as output:
        output.writelines([f'#SEM_file={sem_filename}\n',
                           f'#TF={sem}\n'
                           f'#Baseline={baseline}\n'])
        colnames = '\t'.join(['chrom', 'start', 'end', 'spdi', 'ref', 'alt', 'kmer_coord', 'ref_score', 'alt_score', 'log2fc', 'effect_on_binding'])
        output.write(f'{colnames}\n')

        for line in f:
            variant_info = line.strip().split('\t')
            chrom = variant_info[0]
            start = int(variant_info[1])
            end = int(variant_info[2])
            ref = variant_info[3]
            alt = variant_info[4]
            spdi = get_spdi(CHR_REFSEQ_DICT, chrom, start, ref, alt)
            variant_output = '\t'.join([chrom, str(start), str(end), spdi, ref, alt])

            if ref_mismatch(chrom, end, ref, assembly):
                continue

            kmers = get_kmers(len(mat), fasta, chrom, end, alt)

            annot = 0
            annot_score = 0
            
            if kmers['ref']: # check if kmer list is not empty (due to exclusion of kmers containing N from get_kmer_pairs())
                ref_score, alt_score, kmer_coord = get_best_scores(mat, baseline, kmers)
                annot, annot_score = annotate_variant(ref_score, alt_score, baseline)
                ref_score = 2**ref_score # unlog
                alt_score = 2**alt_score
            
            if (only_report_effects == True) and (annot in ["no_binding", "binding_unchanged"]):
                continue

            variant_output = f'{variant_output}\t{kmer_coord}\t{ref_score}\t{alt_score}\t{annot_score}\t{annot}'        
            output.write(f'{variant_output}\n')

if __name__ == "__main__":

    parser = ArgumentParser(prog='SEM_variant_annotation')
    parser.add_argument('--file', '-f', help='Input variant list', required=True)
    parser.add_argument('--semsdir', '-d', help='Path to directory containing SEMs', required=True)
    parser.add_argument('--sems', '-s', nargs='*', default=None, help='List of sems. If not specified, all SEMs from the semsdir will be used')
    parser.add_argument('--assembly', '-a', help='Path to indexed reference genome', required=True)
    parser.add_argument('--baselines', '-b', help='File containing all the baseline values', required=True)
    parser.add_argument('--n_processes', '-n', type=int, default=1, help='Number of processes for multiprocessing')
    parser.add_argument('--only-report-effects', '-e', action='store_true', help='Only include variants with annotated effect in final output')
    parser.add_argument('--outdir', '-o', help='Output path of variant annotations')
    args = parser.parse_args()

    variants_file = args.file
    sems_dir = args.semsdir
    sems_list = args.sems
    assembly = args.assembly
    baselines_file = args.baselines
    n_processes = args.n_processes
    output_dir = args.outdir

    sems, sem_filenames = load_sems(sems_dir, sems_list)
    baselines = load_baselines(baselines_file, sems_dir)

    # Run annotation with multiprocessing
    with Pool(n_processes) as pool:
        pool.starmap(run_annotation, zip(list(sems.keys()), sem_filenames, repeat(variants_file), repeat(output_dir), repeat(baselines), repeat(assembly), repeat(args.only_report_effects)))
        