import numpy as np

from semvar.io import load_sems
from semvar.scoring import get_best_sem_score

def test_get_best_sem_score():
    sems_dict = load_sems('tests/test_data/test_SEMs')
    usf1_sem = sems_dict['USF1']['mat']
    score, idx = get_best_sem_score('GCCACGTGATAT', usf1_sem)
    np.testing.assert_almost_equal(score, -22.692536, decimal=7)
    assert idx == 0

    # Reverse complement has best score. Seq is same length as SEM so idx should be 0
    dummy_sem = sems_dict['dummy']['mat']
    score, idx = get_best_sem_score('AACT', dummy_sem) 
    assert score == float(18)
    assert idx == 0

    # Fwd has best score at the second window
    score, idx = get_best_sem_score('AGCTT', dummy_sem)
    assert score == float(19)
    assert idx == 1