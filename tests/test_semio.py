import numpy as np

from semvar.semio import _read_motif_file, load_sems, load_baselines

def test_read_motif_file():
    header, mat = _read_motif_file('tests/test_data/test_SEMs/BHLHB2_GM12878.sem')
    assert header == 'BHLHE40\tA\tC\tG\tT'
    true_array = np.array([
        [-0.916001, -0.456257, 0.0449168, -0.830239],
        [-0.952639, -1.48185, 0.00568567, 0.0552664],
        [-1.81496, 0.0460032, -1.98482, -2.0076],
        [0.0959402, -1.80794, -1.28597, -1.83786],
        [-3.23656, 0.0467237, -2.52416, -2.64658],
        [-2.6633, -2.52651, 0.0467237, -3.20208],
        [-1.82936, -1.29516, -1.78077, 0.1013],
        [-2.03469, -1.99447, 0.0458512, -1.79719],
        [0.0461967, -0.0127797, -1.50707, -0.981267],
        [-0.874196, 0.039353, -0.464275, -0.906534]
        ])
    np.testing.assert_array_equal(mat, true_array)
    
def test_load_sems():
    sems, sem_filenames = load_sems('tests/test_data/test_SEMs')
    pass # TODO: Implement

def test_load_baselines():
    baselines = load_baselines('tests/test_data/baselines.txt', 'tests/test_data/test_SEMs')
    assert baselines == {
        'BHLHE40': -2.473448,
        'CREB3L1': -1.009556,
    }