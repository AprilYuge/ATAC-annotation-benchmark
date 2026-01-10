from memory_profiler import memory_usage
save = True
binarize = False

import time
import uniport as up
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
from os import sys, path
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

rna_path = str(sys.argv[1])
atac_path = str(sys.argv[2])
result_folder = str(sys.argv[3])
if not os.path.exists(result_folder):
    os.makedirs(result_folder)

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
    sc.pp.filter_cells(rna, min_genes=200)
    sc.pp.filter_genes(rna, min_cells=3)
    rna.obs['domain_id'] = 1
    rna.obs['domain_id'] = rna.obs['domain_id'].astype('category')
    rna.obs['source'] = 'RNA'

    sc.pp.filter_cells(atac, min_genes=200)
    sc.pp.filter_genes(atac, min_cells=3)
    atac.obs['domain_id'] = 0
    atac.obs['domain_id'] = atac.obs['domain_id'].astype('category')
    atac.obs['source'] = 'ATAC'
    
    # Process RNA and ATAC labels
    rna_label = rna.obs['celltype']
    label_idx_mapping = {}
    unique_labels = np.unique(rna_label)
    for i, name in enumerate(unique_labels):
        label_idx_mapping[name] = i
    print(label_idx_mapping)
    rna_label_int = rna.obs['celltype'].replace(label_idx_mapping)
    atac_label_int = atac.obs['celltype'].replace(label_idx_mapping)
    atac_label_int[~atac_label_int.isin(list(label_idx_mapping.values()))] = -1
    rna.obs['label'] = rna_label_int
    atac.obs['label'] = atac_label_int

    adata_cm = atac.concatenate(rna, join='inner', batch_key='domain_id')
    if 'Spleen' not in rna_path: 
        sc.pp.normalize_total(adata_cm)
        sc.pp.log1p(adata_cm)
    sc.pp.highly_variable_genes(adata_cm, n_top_genes=2000, inplace=False, subset=True)
    up.batch_scale(adata_cm)

    if 'Spleen' not in rna_path: 
        sc.pp.normalize_total(rna)
        sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, inplace=False, subset=True)
    up.batch_scale(rna)

    if 'Spleen' not in rna_path: 
        sc.pp.normalize_total(atac)
        sc.pp.log1p(atac)
    sc.pp.highly_variable_genes(atac, n_top_genes=2000, inplace=False, subset=True)
    up.batch_scale(atac)

    # Run uniPort
    # adata = up.Run(adatas=[atac, rna], adata_cm=adata_cm, lambda_s=1.0, rep_celltype='celltype', outdir=result_folder)
    adata = up.Run_semi(adatas=[atac, rna], adata_cm=adata_cm, lambda_s=1.0, outdir=result_folder, label_name='label')

    # KNN classifier
    X_train = adata[adata.obs.source=='RNA'].obsm['latent']
    y_train = adata[adata.obs.source=='RNA'].obs['celltype']
    X_test = adata[adata.obs.source=='ATAC'].obsm['latent']

    knn = KNeighborsClassifier(n_neighbors=30).fit(X_train, y_train)
    tmp = y_train.unique().tolist()
    tmp.sort()
    
    if save:
        # prob matrix
        prob = pd.DataFrame(knn.predict_proba(X_test), index=adata[adata.obs.source=='ATAC'].obs_names, columns=tmp)
        prob.to_csv(result_folder+'/prob.csv')
        # predicted label
        pred = prob.idxmax(axis=1)
        pred.to_csv(result_folder+'/pred.csv')
        # latent emb
        latent = pd.DataFrame(adata.obsm['latent'], index=adata.obs_names)
        latent.to_csv(result_folder+'/latent.csv')

    end = time.time()
    print('Running time: %.2f sec' % (end-start))
    with open(path.join(result_folder, 'time.txt'), "w") as file:
        file.write(str(end-start))
    
peak_mem_usage = memory_usage(main, max_iterations=1, max_usage=True)
print('Peak memory usage: %.2f MB' % peak_mem_usage)
with open(path.join(result_folder, 'memory.txt'), "w") as file:
    file.write(str(peak_mem_usage))
