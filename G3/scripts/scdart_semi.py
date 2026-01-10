from memory_profiler import memory_usage
save = True
binarize = False

import time
from scDART import scDART_semi
from itertools import chain
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import os
from os import sys, path
from scipy.sparse import csc_matrix, csr_matrix
from scipy.io import mmwrite, mmread
from sklearn.neighbors import KNeighborsClassifier
import re
import torch

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
            counts[string] = 1
        renamed.append(new_name)
    return renamed

rna_path = str(sys.argv[1])
atac_path = str(sys.argv[2])
result_folder = str(sys.argv[3])
if not os.path.exists(result_folder):
    os.makedirs(result_folder)
reg_path = '%s/region2gene.csv' % result_folder

load_dir = '/gpfs/gibbs/pi/zhao/yw599/Multiome/data_pp'    

def main():  
    
    if len(sys.argv) >= 6:
        rna_cell_subset = str(sys.argv[4])
        if rna_cell_subset == '!':
            rna_cell_subset = None
        atac_cell_subset = str(sys.argv[5])
        if atac_cell_subset == '!':
            atac_cell_subset = None
    else:
        rna_cell_subset = None
        atac_cell_subset = None
    if len(sys.argv) >= 7:
        rna_new_annot = sys.argv[6]
        if rna_new_annot == '!':
            rna_new_annot = None
    else:
        rna_new_annot = None
    
    start = time.time()

    rna_counts = csr_matrix(mmread(path.join(rna_path, 'counts.mtx')))
    rna_gene_names = read_txt_np(path.join(rna_path, 'genes.txt'))
    rna_cell_names = read_txt_np(path.join(rna_path, 'cells.txt'))
    rna_label = read_txt_np(path.join(rna_path, 'annotations.txt'))
    if binarize:
        rna_counts = (rna_counts > 0).astype('int')

    rna = sc.AnnData(rna_counts)
    # rna.var_names = [i.upper() for i in rna_gene_names]
    rna.var_names = rna_gene_names
    rna.obs_names = rna_cell_names
    rna.var_names.name = 'Gene'
    rna.obs_names.name = 'CellID'
    rna.obs = pd.DataFrame({'celltype': rna_label}, index=rna_cell_names)
    del rna_counts, rna_gene_names, rna_cell_names
    if rna_cell_subset is not None:
        rna_cell_subset = read_txt_np(rna_cell_subset)
        rna = rna[rna_cell_subset]
        rna.obs_names = rename_duplicates(rna.obs_names)

    atac_counts = csr_matrix(mmread(path.join(atac_path, 'counts.mtx')))
    atac_gene_names = read_txt_np(path.join(atac_path, 'genes.txt'))
    atac_cell_names = read_txt_np(path.join(atac_path, 'cells.txt'))
    atac_label = read_txt_np(path.join(atac_path, 'annotations.txt'))
    if binarize:
        atac_counts = (atac_counts > 0).astype('int')

    atac = sc.AnnData(atac_counts)
    # atac.var_names = [i.upper() for i in atac_gene_names]
    atac.var_names = atac_gene_names
    atac.obs_names = atac_cell_names
    atac.var_names.name = 'Gene'
    atac.obs_names.name = 'CellID'
    atac.obs = pd.DataFrame({'celltype': atac_label}, index=atac_cell_names)
    del atac_counts, atac_gene_names, atac_cell_names
    if atac_cell_subset is not None:
        atac_cell_subset = read_txt_np(atac_cell_subset)
        atac = atac[atac_cell_subset]
        atac.obs_names = rename_duplicates(atac.obs_names)

    # Data preprocessing
    coarse_reg = pd.read_csv(reg_path, sep = ",", index_col = 0)

    sc.pp.filter_cells(rna, min_genes=200)
    sc.pp.filter_genes(rna, min_cells=3)
    if 'Spleen' not in rna_path: 
        sc.pp.normalize_total(rna)
        sc.pp.log1p(rna)
    coarse_reg = coarse_reg.loc[:, coarse_reg.columns.isin(rna.var_names)]
    rna = rna[:, coarse_reg.columns].copy()

    sc.pp.filter_cells(atac, min_genes=200)
    sc.pp.filter_genes(atac, min_cells=3)
    coarse_reg = coarse_reg.loc[coarse_reg.index.isin(atac.var_names)]
    atac = atac[:, coarse_reg.index].copy()
    atac.X = (atac.X > 0).astype(np.float)

    # Process RNA labels
    rna_label = rna.obs['celltype']
    label_idx_mapping = {}
    unique_labels = np.unique(rna_label)
    for i, name in enumerate(unique_labels):
        label_idx_mapping[name] = i
    print(label_idx_mapping)
    rna.obs['label'] = rna.obs['celltype'].replace(label_idx_mapping)
    
    # Run scDART
    scDART_op = scDART_semi(n_epochs = 500, use_anchor = False, k = 10, batch_size=128,
                            device = torch.device('cuda' if torch.cuda.is_available() else 'cpu'))
    scDART_op = scDART_op.fit(rna_count = rna.X.todense(), rna_label = rna.obs['label'].values, 
                              atac_count = atac.X.todense(), reg = coarse_reg.values,
                              rna_anchor = None, atac_anchor = None)
    z_rna, z_atac = scDART_op.transform(rna_count = rna.X.todense(), atac_count = atac.X.todense(), rna_anchor = None, atac_anchor = None)

    # KNN classifier
    X_train = z_rna
    y_train = rna.obs['celltype']
    X_test = z_atac

    knn = KNeighborsClassifier(n_neighbors=30).fit(X_train, y_train)
    tmp = y_train.unique().tolist()
    tmp.sort()
    
    if save:
        # prob matrix
        prob = pd.DataFrame(knn.predict_proba(X_test), index=atac.obs_names, columns=tmp)
        prob.to_csv(result_folder+'/prob.csv')
        # predicted label
        pred = prob.idxmax(axis=1)
        pred.to_csv(result_folder+'/pred.csv')
        # latent emb
        latent_rna = pd.DataFrame(z_rna, index=rna.obs_names)
        latent_atac = pd.DataFrame(z_atac, index=atac.obs_names)
        latent = pd.concat([latent_rna, latent_atac])
        latent.to_csv(result_folder+'/latent.csv')

    end = time.time()
    print('Running time: %.2f sec' % (end-start))
    with open(path.join(result_folder, 'time.txt'), "w") as file:
        file.write(str(end-start))
    
peak_mem_usage = memory_usage(main, max_iterations=1, max_usage=True)
print('Peak memory usage: %.2f MB' % peak_mem_usage)
with open(path.join(result_folder, 'memory.txt'), "w") as file:
    file.write(str(peak_mem_usage))
