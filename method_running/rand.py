import numpy as np
import pandas as pd
import sys
import os
from os import sys, path
from scipy.sparse import csr_matrix
from scipy.io import mmread

print(len(sys.argv))
assert len(sys.argv) == 4 or len(sys.argv) == 6 or len(sys.argv) == 7, "parameters needed: rna path, result folder"
rna_path = str(sys.argv[1])
atac_path = str(sys.argv[2])
result_folder = str(sys.argv[3])
if len(sys.argv) >= 6:
    subset_rna = str(sys.argv[4])
    if subset_rna == '!':
        subset_rna = None
    subset_atac = str(sys.argv[5])
    if subset_atac == '!':
        subset_atac = None
else:
    subset_rna = None
    subset_atac = None
if len(sys.argv) >= 7:
    rna_new_annot = sys.argv[6]
    if rna_new_annot == '!':
        rna_new_annot = None
else:
    rna_new_annot = None

def read_txt_np(filename):
    with open(filename) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
        return np.array(lines)    

def run_random(rna_path, atac_path, result_folder, subset_rna, subset_atac, rna_new_annot):
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    if result_folder[-1] != '/':
        result_folder += '/'

    rna_cell_names = read_txt_np(path.join(rna_path, 'cells.txt'))
    atac_cell_names = read_txt_np(path.join(atac_path, 'cells.txt'))
    if rna_new_annot is not None:
        rna_label = read_txt_np(rna_new_annot)
    else:
        rna_label = read_txt_np(path.join(rna_path, 'annotations.txt'))
    
    ## subset
    if subset_rna is not None:
        subset_rna_barcodes = read_txt_np(subset_rna)
        ids_series = pd.Series(np.arange(len(rna_cell_names)), index=rna_cell_names)
        idx_bc_rna = ids_series[subset_rna_barcodes]
        rna_label = rna_label[idx_bc_rna]
        rna_cell_names = subset_rna_barcodes
    if subset_atac is not None:
        subset_atac_barcodes = read_txt_np(subset_atac)
        ids_series = pd.Series(np.arange(len(atac_cell_names)), index=atac_cell_names)
        idx_bc_atac = ids_series[subset_atac_barcodes]
        atac_cell_names = subset_atac_barcodes
    
    dic = {}
    for ct in set(rna_label):
        dic[ct] = np.sum(rna_label == ct)/len(rna_label)
    prob = np.array(list(dic.values())).repeat(len(atac_cell_names)).reshape(-1,len(atac_cell_names)).transpose()
    prob = pd.DataFrame(prob, index=atac_cell_names, columns=list(dic.keys()))
    prob.to_csv(result_folder+'prob.csv')
    # predicted label
    pred = prob.idxmax(axis=1)
    pred.to_csv(result_folder+'pred.csv')

run_random(rna_path, atac_path, result_folder, subset_rna, subset_atac, rna_new_annot)