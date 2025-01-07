import numpy as np

from utils import revdnacomp, char_to_num

def _get_max_score(num_seq, mat, len_mat, max_score):
    windows = np.lib.stride_tricks.sliding_window_view(num_seq, len_mat)
    bit_scores = mat[np.arange(len_mat)[:, None], windows.T].sum(axis=0)
    current_max = bit_scores.max()
    if current_max > max_score:
        return current_max
    else:
        return max_score

def _get_sem_sum_score(seq, mat):
    num_seq = char_to_num(seq)
    max_score = -10000000

    for _ in range(2):
        max_score = _get_max_score(num_seq, mat, len(mat), max_score)
        num_seq = char_to_num(revdnacomp(seq))

    return max_score

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
        ref_score = _get_sem_sum_score(ref, mat)
        alt_score = _get_sem_sum_score(alt, mat)
        if (ref_score > baseline) or (alt_score > baseline): # Now checking exactly which kmer pairs have binding
            abs_score_diff = abs(alt_score-ref_score)
            if abs_score_diff > max_score_diff:
                max_score_diff = abs_score_diff
                chosen_ref_score = ref_score
                chosen_alt_score = alt_score
                chosen_kmer_coord = coord

    return chosen_ref_score, chosen_alt_score, chosen_kmer_coord