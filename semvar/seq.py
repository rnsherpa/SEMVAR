def get_sequences(sem_len, fasta, chr, pos, ref, alt):
    '''Get all possible reference and variant k-mers overlapping the variant position where k is the length of the SEM
    Params
        sem_len (int): length of SEM
        fasta (pyfaidx.Fasta): Pyfaidx Fasta object of the reference genome assembly
        chr (str): Chromosome of SNV in 'chrN' format
        pos (int): 1-based position of SNV
        ref (str): Reference allele
        alt (str): Alternate allele
    Returns
        kmers (dict of lists): Lists of k-mers and their coordinates
    '''
    seq_start = pos-sem_len # currently pos is 1-based
    seq_end = pos+sem_len-1
    ref_seq = fasta[chr][seq_start:seq_end].seq.upper()
    middle_index = sem_len - 1
    alt_seq = ref_seq[:middle_index] + alt + ref_seq[middle_index + 1:]

    return ref_seq, alt_seq