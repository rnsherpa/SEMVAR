import numpy as np

def char_to_num(sequence):
    char_list = ['A', 'C', 'G', 'T']
    invalid_char_num = 7
    numeric_sequence = []
    for char in sequence:
        if char in char_list:
            numeric_sequence.append(char_list.index(char))
        else:
            numeric_sequence.append(invalid_char_num) # Marking invalid characters such as 'N' to be 7
    return np.array(numeric_sequence, dtype=np.int32)

def contains_invalid_chars(sequence):
    valid_bases = {'A', 'C', 'T', 'G'}
    if not set(sequence).issubset(valid_bases):
        return True
    else:
        return False
    
def no_valid_kmers(sequence, sem_len):
    valid_bases = {'A', 'C', 'T', 'G'}
    current_length = 0
    
    for char in sequence:
        if char in valid_bases:
            current_length += 1
            if current_length >= sem_len:
                return False
        else:
            current_length = 0
    
    return True

def revdnacomp(dna):
    revcomp = dna[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
    return revcomp

def get_spdi(chr_refseq_dict, chrom, start, ref, alt):
    '''
    Outputs SNV identification in SPDI format: https://academic.oup.com/bioinformatics/article/36/6/1902/5628222
    '''
    spdi = f'{chr_refseq_dict[chrom]}:{start}:{ref}:{alt}'
    return spdi

def ref_mismatch(chr, pos, ref, ref_fasta):
    '''Check for mismatches of ref allele in variant list to reference genome
    Params
        chr (str): Chromosome of SNV in 'chrN' format
        pos (int): 1-based position of SNV
        ref (str): Reference allele at position from variant list
        ref_fasta (pyfaidx.Fasta): Pyfaidx Fasta object of the reference genome assembly
    Returns
        (bool): True if variant list ref does not match fasta ref, False if it does
    '''
    if ref != ref_fasta[chr][pos-1:pos].seq.upper():
        return True
    else: return False