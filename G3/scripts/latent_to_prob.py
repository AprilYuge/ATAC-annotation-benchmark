import pandas as pd
import numpy as np
import os
from os import sys, path
from sklearn.neighbors import KNeighborsClassifier

def read_txt_np(filename):
    with open(filename) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
        return np.array(lines)

def rename_duplicates(strings):
    renamed = []
    counts = {}
    for string in strings:
        if string in counts:
            count = counts[string]
            new_name = f"{string}_{count}"
            counts[string] += 1
        else:
            new_name = string
            counts[string] = 2
        renamed.append(new_name)
    return renamed

tissue = str(sys.argv[1])
method = str(sys.argv[2])
result_folder = str(sys.argv[3])
if len(sys.argv) >= 5:
    rna_cell_subset = str(sys.argv[4])
    if rna_cell_subset == '!':
        rna_cell_subset = None
else:
    rna_cell_subset = None

dic_map_rna = {'PBMC': 'RNA', 'BMMC': 'RNA_20k', 'PBMC-COVID': 'RNA_20k', 'BCC': 'RNA_20k', 'Endo': 'RNA', 
               'FetalAtlas': 'RNA_20k', 'Heart': 'RNA_20k', 'HSPC': 'RNA', 'Kidney': 'RNA', 'kidney': 'RNA_20k', 
               'MouseBrain_10x': 'RNA_20k', 'MouseBrain_snap': 'RNA_20k', 'MOp': 'RNA_20k', 'MouseEmbryo': 'RNA_20k', 
               'MouseAtlas-brain': 'FACS', 'MouseAtlas_facs': 'FACS_20k', 'MouseAtlas_droplet': 'droplet_20k', 
               'MouseSpleen': 'RNA', 'PBMC_paired': '10x/RNA', 'BMMC_paired': 'multiome/RNA', 'HSPC_paired': '10x/RNA', 
               'TDBM_paired': '10x/RNA', 'kidney_paired': 'scicar_annotated/RNA', 'MouseRetina_paired': '10x/RNA', 
               'MouseSkin_paired': 'share/RNA_20k', 'MouseEmbryo_paired': '10x/RNA_20k', 'MOp_paired': 'snare/RNA'}
method_map = {'StabMap': ['latent_to_rna.csv', 'latent_to_atac.csv'], 
              'StabMap_semi': ['latent_to_rna.csv', 'latent_to_rna_LD.csv', 'latent_to_atac.csv'], 
              'LIGER': ['latent.csv', 'latent_useRNA.csv'], 
              'UINMF': ['latent.csv'], 'scMC': ['latent.csv'], 'bindSC': ['latent.csv'],
              'StabMap_multi': ['latent.csv'],
              'StabMap_multi_semi': ['latent.csv']}
    
load_dir = '/gpfs/gibbs/pi/zhao/yw599/Multiome/data_pp'    
    
rna_path = '%s/%s/%s' % (load_dir, tissue.split('_')[0], dic_map_rna[tissue])

rna_cell_names = read_txt_np(path.join(rna_path, 'cells.txt'))
if (method != 'StabMap_multi') and (method != 'StabMap_multi_semi'):
    rna_cell_names = [i+'-rna' for i in rna_cell_names]
rna_label = read_txt_np(path.join(rna_path, 'annotations.txt'))
rna_label = pd.DataFrame(rna_label, index=rna_cell_names)

if rna_cell_subset is not None:
    rna_label = rna_label.loc[rna_cell_subset]
    rna_label.index = rename_duplicates(rna_label.index)

for file in method_map[method]:
    # Get suffix
    suffix = file[6:-4]
    
    path = '%s/%s' % (result_folder, file)
    latent = pd.read_csv(path, index_col=0)

    # KNN classifier
    X_train = latent.loc[latent.index.isin(rna_label.index)]
    y_train = rna_label.loc[X_train.index].values.ravel()
    # X_train = X_train.values
    X_test = latent.loc[~latent.index.isin(rna_label.index)]

    knn = KNeighborsClassifier(n_neighbors=30).fit(X_train, y_train)
    tmp = np.unique(y_train).tolist()
    tmp.sort()

    # prob matrix
    prob = pd.DataFrame(knn.predict_proba(X_test), index=X_test.index, columns=tmp)
    prob.to_csv(result_folder+'/prob%s.csv' % suffix)
    # predicted label
    pred = prob.idxmax(axis=1)
    pred.to_csv(result_folder+'/pred%s.csv' % suffix)
