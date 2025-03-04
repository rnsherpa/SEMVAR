import os
from collections import defaultdict

import numpy as np

def _read_motif_file(motif_file):
    mat = []
    with open(motif_file) as file:
        baseline_header = file.readline().strip()
        baseline = float(baseline_header.split(':')[1])

        tf_header = file.readline().strip()
        tf_name = tf_header.split('\t')[0]

        for line in file:
            wt = line.strip().split("\t")
            mat.append([float(i) for i in wt[1:]])
    return tf_name, np.array(mat), baseline

def load_sems(sems_dir, sems_list = None):
    '''Load SEM matrices, baselines, and filenames from a directory into a dictionary
    Params:
        sems_dir (str): Path to the directory containing SEMs
        sems_list (list(str) | None): List of SEMs to generate annotations with. SEMs identified by TF name.
    Returns:
        sems_dicts (defaultdict(dict)): Keys are TF names and values are dicts containing SEM matrices, baselines, and filenames
    '''
    sems_dicts = defaultdict(dict)
    for filename in sorted(os.listdir(sems_dir)):
        tf_name, mat, baseline = _read_motif_file(os.path.join(sems_dir, filename))

        if sems_list is not None:
            if tf_name in sems_list:
                    sems_dicts[tf_name]['mat'] = mat
                    sems_dicts[tf_name]['filename'] = filename
                    sems_dicts[tf_name]['baseline'] = baseline
        else:
            sems_dicts[tf_name]['mat'] = mat
            sems_dicts[tf_name]['filename'] = filename
            sems_dicts[tf_name]['baseline'] = baseline

    return sems_dicts