import numpy as np

from semvar.io import load_sems
from semvar.scoring import get_best_sem_score

def test_get_best_sem_score():
    sems, _ = load_sems('tests/test_data/test_SEMs')
    mat = sems['BHLHE40']
    np.testing.assert_almost_equal(get_best_sem_score('GCCACGTGAT', mat), -1.9147285, decimal=7)
    np.testing.assert_almost_equal(get_best_sem_score('ACCACGTGAT', mat), -2.8756463, decimal=7)