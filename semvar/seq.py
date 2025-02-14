def get_sequences(sem_len, fasta, chr, pos, ref, alt):
    '''Get the sequence windows around the variant to be scanned by the SEM
    Params
        sem_len (int): length of SEM
        fasta (pyfaidx.Fasta): Pyfaidx Fasta object of the reference genome assembly
        chr (str): Chromosome of SNV in 'chrN' format
        pos (int): 1-based position of SNV
        ref (str): Reference allele # To be used in indel implementation
        alt (str): Alternate allele
    Returns
        ref_seq, alt_seq (str, str): Ref and alt sequences to be scanned by the SEM
    '''
    seq_start = pos-sem_len # pos is 1-based
    seq_end = pos+sem_len-1
    ref_seq = fasta[chr][seq_start:seq_end].seq.upper()
    middle_index = sem_len - 1
    alt_seq = ref_seq[:middle_index] + alt + ref_seq[middle_index + 1:]

    return ref_seq, alt_seq