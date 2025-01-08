import numpy as np

from semvar.utils import revdnacomp, char_to_num

def _get_max_score(num_seq, mat, len_mat, max_score):
    windows = np.lib.stride_tricks.sliding_window_view(num_seq, len_mat)
    windows = windows[np.all(windows != 7, axis=1)]
    bit_scores = mat[np.arange(len_mat)[:, None], windows.T].sum(axis=0)
    current_max = bit_scores.max()
    if current_max > max_score:
        return current_max
    else:
        return max_score

def get_best_sem_score(seq, mat):
    max_score = -10000000

    num_seq = char_to_num(seq)
    max_score = _get_max_score(num_seq, mat, len(mat), max_score)

    num_seq_revcomp = char_to_num(revdnacomp(seq))
    max_score = _get_max_score(num_seq_revcomp, mat, len(mat), max_score)

    return max_score