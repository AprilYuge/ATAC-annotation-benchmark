from memory_profiler import memory_usage
save = True
binarize = False

import time
from jamie import JAMIE
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

    # Read in RNA data
    rna_counts = csr_matrix(mmread(path.join(rna_path, 'counts.mtx')))
    if binarize:
        rna_counts = (rna_counts > 0).astype('int')
    rna_gene_names = read_txt_np(path.join(rna_path, 'genes.txt'))
    rna_cell_names = read_txt_np(path.join(rna_path, 'cells.txt'))
    rna_label = read_txt_np(path.join(rna_path, 'annotations.txt'))

    rna = sc.AnnData(rna_counts)
    rna.var_names = rna_gene_names
    rna.obs_names = rna_cell_names
    rna.var_names.name = 'Gene'
    rna.obs_names.name = 'CellID'
    rna.obs = pd.DataFrame({'celltype': rna_label}, index=rna_cell_names)
    rna.obs['dataset'] = 'RNA'
    del rna_counts, rna_gene_names, rna_cell_names
    
    # Read in ATAC data
    atac_counts = csr_matrix(mmread(path.join(atac_path, 'counts.mtx')))
    if binarize:
        atac_counts = (atac_counts > 0).astype('int')
    atac_gene_names = read_txt_np(path.join(atac_path, 'genes.txt'))
    atac_cell_names = read_txt_np(path.join(atac_path, 'cells.txt'))
    atac_label = read_txt_np(path.join(atac_path, 'annotations.txt'))

    atac = sc.AnnData(atac_counts)
    atac.var_names = atac_gene_names
    atac.obs_names = atac_cell_names
    atac.var_names.name = 'Gene'
    atac.obs_names.name = 'CellID'
    atac.obs = pd.DataFrame({'celltype': atac_label}, index=atac_cell_names)
    # atac.obs['label'] = -1
    atac.obs['dataset'] = 'ATAC'
    del atac_counts, atac_gene_names, atac_cell_names
        
    # Data preprocessing
    sc.pp.filter_cells(rna, min_genes=200)
    sc.pp.filter_genes(rna, min_cells=3)
    sc.pp.filter_cells(atac, min_genes=200)
    sc.pp.filter_genes(atac, min_cells=3)
    sc.pp.scale(rna)
    sc.pp.scale(atac)
    dataset = [rna.X, atac.X]
    
    # Priors
    priors = np.zeros((len(dataset[0]), len(dataset[1])))
    
    # Parameters
    reduced_dim = 32
    kwargs = {
        'output_dim': reduced_dim,
        'epoch_DNN': 10000,
        'min_epochs': 2500,
        'log_DNN': 500,
        'use_early_stop': True,
        'batch_size': 512,
        'pca_dim': 2*[512],
        'dist_method': 'euclidean',
        'loss_weights': [1,1,1,1],
        # 'loss_weights': [1,1,1,0],
        # 'use_f_tilde': False,
        'dropout': 0,
        'enable_memory_logging': True,
    }
    
    # Train JAMIE
    jm = JAMIE(**kwargs, debug=True)
    jm_data = jm.fit_transform(dataset=dataset, P=priors)

    # Get embeddings
    X_train = jm_data[0]
    X_test = jm_data[1]

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
