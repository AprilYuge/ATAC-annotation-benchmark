from memory_profiler import memory_usage
save = True
binarize = False

import time
import sys
sys.path.append('/gpfs/gibbs/pi/zhao/yw599/Multiome/SCOT/src')
from os import sys, path
version = str(sys.argv[1])
if version == 'v1':
    from scotv1 import *
else:
    from scotv2 import *

import scglue
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
from scipy.sparse import csc_matrix, csr_matrix
from scipy.io import mmwrite, mmread
from sklearn.neighbors import KNeighborsClassifier
import re

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

print(len(sys.argv))
assert len(sys.argv) == 5 or len(sys.argv) == 7 or len(sys.argv) == 8, "parameters needed: SCOT version, rna path, atac path, result folder"
rna_path = str(sys.argv[2])
atac_path = str(sys.argv[3])
result_folder = str(sys.argv[4])
if not os.path.exists(result_folder):
    os.makedirs(result_folder)

def main():     
    
    if len(sys.argv) >= 7:
        rna_cell_subset = str(sys.argv[5])
        if rna_cell_subset == '!':
            rna_cell_subset = None
        atac_cell_subset = str(sys.argv[6])
        if atac_cell_subset == '!':
            atac_cell_subset = None
    else:
        rna_cell_subset = None
        atac_cell_subset = None
    if len(sys.argv) >= 8:
        rna_new_annot = sys.argv[7]
        if rna_new_annot == '!':
            rna_new_annot = None
    else:
        rna_new_annot = None
        
    load_dir = '/gpfs/gibbs/pi/zhao/yw599/Multiome/data_pp'  
    
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
    sc.pp.filter_cells(rna, min_genes=200)
    sc.pp.filter_genes(rna, min_cells=3)
    if 'Spleen' not in rna_path: 
        sc.pp.normalize_total(rna)
        sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=2000)
    rna = rna[:, rna.var.highly_variable]
    sc.pp.scale(rna)
    sc.pp.pca(rna, n_comps=50)

    sc.pp.filter_cells(atac, min_genes=200)
    sc.pp.filter_genes(atac, min_cells=3)
    scglue.data.lsi(atac, n_components=50)

    if version == 'v1':
        scot= SCOT(rna.obsm['X_pca'], atac.obsm['X_lsi'])
        X_train, X_test = scot.align(k=50, XontoY=False)
    else:
        scot= SCOTv2([rna.obsm['X_pca'], atac.obsm['X_lsi']])
        X_train, X_test = scot.align(k=50)

    # KNN classifier
    y_train = rna.obs['celltype']

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
        latent_rna = pd.DataFrame(X_train, index=rna.obs_names)
        latent_atac = pd.DataFrame(X_test, index=atac.obs_names)
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
