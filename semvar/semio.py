import os
import numpy as np

def _read_motif_file(motif_file):
    mat = []
    with open(motif_file) as file:
        header = file.readline().strip()
        for line in file:
            wt = line.strip().split("\t")
            mat.append([float(i) for i in wt[1:]])
    return header, np.array(mat)

def load_sems(sems_dir, sems_list = None):
    '''Load all SEMs in a directory into a dictionary, ensuring only one SEM is used per TF 
    Params:
        sems_dir (str): Path to the directory containing SEMs
        sems_list (list(str) | None): List of SEMs to generate annotations with. SEMs identified by TF name.
    Returns:
        sems (dict): Keys are TF names and values are the corresponding matrices
        sem_filenames (list): List of SEM filenames 
    '''
    sems = {}
    sem_filenames = []
    for filename in sorted(os.listdir(sems_dir)):
        header, mat = _read_motif_file(os.path.join(sems_dir, filename))
        tf_name = header.split('\t')[0] # Use first identifier
        if sems_list is not None:
            if tf_name not in sems.keys():
                if tf_name in sems_list:
                    sems[tf_name] = mat
                    sem_filenames.append(filename)
            else:
                print(f'{tf_name} appears to have multiple SEMs in the directory. Only using first instance.')
        else:
            if tf_name not in sems.keys():
                sems[tf_name] = mat
                sem_filenames.append(filename)
            else:
                print(f'{tf_name} appears to have multiple SEMs in the directory. Only using first instance.')
    return sems, sem_filenames

def load_baselines(baselines_file, sems_dir):
    '''Load all baselines from file to a dictionary
    Params:
        baselines_file (str): Path to file containing baselines
        sems_dir (str): Path to the directory containing SEMs
    Returns:
        baselines (dict): Key are motifs and values are the corresponding baselines
    '''
    motif_name_map = {}
    for filename in sorted(os.listdir(sems_dir)):
        header, _ = _read_motif_file(os.path.join(sems_dir, filename))
        tf_name = header.split('\t')[0]
        filename_no_ext = os.path.splitext(filename)[0]
        motif_name_map[filename_no_ext] = tf_name

    baselines = {}
    with open(baselines_file) as f:
            for line in f:
                line = line.strip()
                motif_name = line.split('\t')[0]
                if motif_name in motif_name_map.keys():
                    tf_name = motif_name_map[motif_name] # map from filename to TF name
                else:
                    continue
                baseline = line.split('\t')[1]
                if tf_name not in baselines.keys():
                    baselines[tf_name] = float(baseline)
                else:
                    print(f'{tf_name} appears to have multiple baseline values in the file. Only using first instance.')
    return baselines