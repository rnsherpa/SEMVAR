import numpy as np

def char_to_num(sequence):
    char_list = ['A', 'C', 'G', 'T']
    numeric_sequence = []
    for char in sequence:
        if char in char_list:
            numeric_sequence.append(char_list.index(char))
        else:
            raise ValueError(f"Unknown character! {sequence}")
    return np.array(numeric_sequence, dtype=np.int32)

def revdnacomp(dna):
    revcomp = dna[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
    return revcomp

def get_spdi(chr_refseq_dict, chrom, start, ref, alt):
    '''
    Outputs SNV identification in SPDI format: https://academic.oup.com/bioinformatics/article/36/6/1902/5628222
    '''
    spdi = f'{chr_refseq_dict[chrom]}:{start}:{ref}:{alt}'
    return spdi

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
    if list_ref != ref_fasta[chr][pos-1:pos].seq.upper():
        return True
    else: return False