import numpy as np

from semvar.io import load_sems
    
def test_load_sems():
    sems_dict = load_sems('tests/test_data/test_SEMs')

    true_dummy_sem = np.array([
        [1, 2, 3, 4],
        [2, 3, 4, 5],
        [3, 4, 5, 6],
        [4, 5, 6, 7]
        ])
    np.testing.assert_array_equal(sems_dict['dummy']['mat'], true_dummy_sem)

    true_usf1_sem = np.array([
        [-0.644176,	-0.301248,	0.088734,	-0.942887],
        [-0.724766,	-1.162370,	0.091925,	-1.881332],
        [-1.664509,	-0.805443,	-1.402396,	0.127123],
        [-2.931446,	0.046641,	-2.213511,	-3.054659],
        [0.078802,	-2.417371,	-1.752331,	-2.696168],
        [-3.281289,	0.096484,	-2.583739,	-1.974596],
        [-1.963534,	-2.570744,	0.119178,	-3.271426],
        [-2.707923,	-1.755313,	-2.434221,	0.078391],
        [-3.093634,	-2.237835,	0.045490,	-2.936693],
        [0.165581,	-1.350345,	-0.755858,	-1.618773],
        [-1.853719,	0.119772,	-1.124698,	-0.702419],
        [-0.980905,	0.073014,	-0.339407,	-0.663971]
        ])
    np.testing.assert_array_equal(sems_dict['USF1']['mat'], true_usf1_sem)

    assert sems_dict['USF1']['baseline'] == -2.6771622537803843

    assert sems_dict['dummy']['filename'] == 'dummy.sem'
    assert sems_dict['CTCF']['filename'] == 'CTCF_HUMAN.HepG2.sem'
    assert sems_dict['USF1']['filename'] == 'USF1_HUMAN.HepG2.sem'
