import os
from multiprocessing import Pool
from itertools import repeat
from argparse import ArgumentParser
from pyfaidx import Fasta

def read_motif_file(motif_file):
    mat = []
    with open(motif_file) as file:
        header = file.readline().strip()
        for line in file:
            wt = line.strip().split("\t")
            mat.append([float(i) for i in wt[1:]])
    return header, mat

def load_sems(sems_dir, sems_list = None):
    '''Load all SEMs in a directory into a dictionary, ensuring only one SEM is used per TF 
    Params:
        sems_dir (str): Path to the directory containing SEMs
        sems_list (list(str) | None): List of SEMs to generate annotations with. SEMs identified by TF name.
    Returns:
        sems (dict): Keys are TF names and values are the corresponding matrices
        sem_filenames (list): List of SEM filenames 
    '''
    sems = {}
    sem_filenames = []
    for filename in sorted(os.listdir(sems_dir)):
        header, mat = read_motif_file(os.path.join(sems_dir, filename))
        tf_name = header.split('\t')[0] # Use first identifier
        if sems_list is not None:
            if tf_name not in sems.keys():
                if tf_name in sems_list:
                    sems[tf_name] = mat
                    sem_filenames.append(filename)
            else:
                ValueError(f'{tf_name} appears to have multiple SEMs in the directory. Only using first instance.')
        else:
            if tf_name not in sems.keys():
                sems[tf_name] = mat
                sem_filenames.append(filename)
            else:
                ValueError(f'{tf_name} appears to have multiple SEMs in the directory. Only using first instance.')
    return sems, sem_filenames

def load_baselines(baselines_file, sems_dir):
    '''Load all baselines from file to a dictionary
    Params:
        baselines_file (str): Path to file containing baselines
        sems_dir (str): Path to the directory containing SEMs
    Returns:
        baselines (dict): Key are motifs and values are the corresponding baselines
    '''
    motif_name_map = {}
    for filename in sorted(os.listdir(sems_dir)):
        header, _ = read_motif_file(os.path.join(sems_dir, filename))
        tf_name = header.split('\t')[0]
        filename_no_ext = os.path.splitext(filename)[0]
        motif_name_map[filename_no_ext] = tf_name

    baselines = {}
    with open(baselines_file) as f:
            for line in f:
                line = line.strip()
                motif_name = line.split('\t')[0]
                if motif_name in motif_name_map.keys():
                    tf_name = motif_name_map[motif_name] # map from filename to TF name
                else:
                    continue
                baseline = line.split('\t')[1]
                if tf_name not in baselines.keys():
                    baselines[tf_name] = float(baseline)
                else:
                    ValueError(f'{tf_name} appears to have multiple baseline values in the file. Only using first instance.')
    return baselines

def char_to_num(sequence):
    char_list = ['A', 'C', 'G', 'T']
    numeric_sequence = []
    for char in sequence:
        if char in char_list:
            numeric_sequence.append(char_list.index(char))
        else:
            raise ValueError(f"Unknown character! {sequence}")
    return numeric_sequence

def revdnacomp(dna):
    revcomp = dna[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
    return revcomp

def get_sem_sum_score(seq, mat):
    NUM_SEQ = char_to_num(seq)

    max_score = -10000000

    for k in range(2):
        for i in range(len(NUM_SEQ) - len(mat) + 1):
            bit_score = sum(mat[j][NUM_SEQ[i + j]] for j in range(len(mat)))
            if bit_score > max_score:
                max_score = bit_score

        NUM_SEQ = char_to_num(revdnacomp(seq))

    return max_score

def ref_mismatch(chr, pos, list_ref, ref_fasta):
    '''Check for Ns in ref allele column in variant list
    Params
        chr (str): Chromosome of SNV in 'chrN' format
        pos (int): 1-based position of SNV
        list_ref (str): Reference allele at position from variant list
        ref_fasta (path): Path to the reference genome fasta (must have corrsponding .fai index file in the same directory)
    Returns
        (bool): True if variant list ref does not match fasta ref, False if it does
    '''
    fasta = Fasta(ref_fasta)
    if list_ref != fasta[chr][pos-1:pos].seq.upper():
        return True
    else: return False

def get_kmer_pairs(sem_len, fasta, chr, pos, alt):
    '''Get all possible reference and variant k-mers pairs overlapping the variant position where k is the length of the SEM
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

def get_best_scores(mat, baseline, kmers):
    '''Go through each possible kmer pair and return the ref and alt scores with the largest absolute difference
    Params
        mat (list of lists): SEM
        baseline (float): SEM sum score of randomly scrambled motif 
        kmers (dict of lists): Dictionary of ref and alt k-mers
    Returns
        chosen_ref_score (float): Ref score of chosen k-mer pair
        chosen_alt_score (float): Alt score of chosen k-mer pair
        chosen_kmer_coord (str): Coordinate of chosen k-mer pair
    '''
    max_score_diff = 0
    chosen_ref_score = -10000000
    chosen_alt_score = -10000000
    chosen_kmer_coord = '.'
    for ref, alt, coord in zip(kmers['ref'], kmers['alt'], kmers['coords']):
        ref_score = get_sem_sum_score(ref, mat)
        alt_score = get_sem_sum_score(alt, mat)
        if (ref_score > baseline) or (alt_score > baseline): # Now checking exactly which kmer pairs have binding
            abs_score_diff = abs(alt_score-ref_score)
            if abs_score_diff > max_score_diff:
                max_score_diff = abs_score_diff
                chosen_ref_score = ref_score
                chosen_alt_score = alt_score
                chosen_kmer_coord = coord

    return chosen_ref_score, chosen_alt_score, chosen_kmer_coord

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
    chr_refseq_dict = get_chr_refseq_dict('data/refseq_chr_map.txt')
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
            spdi = get_spdi(chr_refseq_dict, chrom, start, ref, alt)
            variant_output = '\t'.join([chrom, str(start), str(end), spdi, ref, alt])

            if ref_mismatch(chrom, end, ref, assembly):
                continue

            kmers = get_kmer_pairs(len(mat), fasta, chrom, end, alt)

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

def get_chr_refseq_dict(mapping_file):
    '''
    Keys are chromosome in chrN format (eg. chr1), values are refseq equivalent (eg. NC_000001.11)
    '''
    chr_refseq_dict = {}
    with open(mapping_file) as f:
        for line in f:
            line = line.strip().split('\t')
            chr_refseq_dict[line[1]] = line[2]
    return chr_refseq_dict

def get_spdi(chr_refseq_dict, chrom, start, ref, alt):
    '''
    Outputs SNV identification in SPDI format: https://academic.oup.com/bioinformatics/article/36/6/1902/5628222
    '''
    spdi = f'{chr_refseq_dict[chrom]}:{start}:{ref}:{alt}'
    return spdi

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
        