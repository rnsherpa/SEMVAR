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