import numpy as np

from semvar.utils import revdnacomp, char_to_num

def _get_max_score(num_seq, mat, len_mat):
    windows = np.lib.stride_tricks.sliding_window_view(num_seq, len_mat)
    windows = windows[np.all(windows != 7, axis=1)]
    bit_scores = mat[np.arange(len_mat)[:, None], windows.T].sum(axis=0)
    max_score = bit_scores.max()
    max_score_idx = np.argmax(bit_scores)
    return max_score, max_score_idx, len(windows)

def get_best_sem_score(seq, mat):
    num_seq = char_to_num(seq)
    fwd_score, fwd_score_idx, _ = _get_max_score(num_seq, mat, len(mat))

    num_seq_revcomp = char_to_num(revdnacomp(seq))
    rvs_score, rvs_score_idx, len_windows = _get_max_score(num_seq_revcomp, mat, len(mat))

    best_score = max(fwd_score, rvs_score)
    if best_score == fwd_score:
        return fwd_score, fwd_score_idx
    else:
        fwd_score_idx = len_windows - 1 - rvs_score_idx # Want to retrive the index in terms of the fwd sequence in order to calculate coordinates on ref genome
        return rvs_score, fwd_score_idx