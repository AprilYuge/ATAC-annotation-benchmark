from memory_profiler import memory_usage
save = True
binarize = False

import time
import numpy as np
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix
from scipy.io import mmread
import os
import pandas as pd
import anndata as ad
from sklearn.neighbors import KNeighborsClassifier
from os import sys, path

rna_path = str(sys.argv[1])
atac_path = str(sys.argv[2])
result_folder = str(sys.argv[3])
if not os.path.exists(result_folder):
    os.makedirs(result_folder)
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

def run_scvi(rna_path, atac_path, result_folder, rna_cell_subset, atac_cell_subset, rna_new_annot):

    # Read in RNA data
    rna_counts = csr_matrix(mmread(path.join(rna_path, 'counts.mtx')))
    rna_gene_names = read_txt_np(path.join(rna_path, 'genes.txt'))
    rna_cell_names = read_txt_np(path.join(rna_path, 'cells.txt'))
    rna_label = read_txt_np(path.join(rna_path, 'annotations.txt'))
    if binarize:
        rna_counts = (rna_counts > 0).astype('int')

    rna = sc.AnnData(rna_counts)
    rna.var_names = rna_gene_names
    rna.obs_names = rna_cell_names
    rna.var_names.name = 'Gene'
    rna.obs_names.name = 'CellID'
    rna.obs = pd.DataFrame({'celltype': rna_label}, index=rna_cell_names)
    rna.var.modality = 'Gene Expression'
    del rna_counts, rna_gene_names, rna_cell_names
    if rna_cell_subset is not None:
        rna_cell_subset = read_txt_np(rna_cell_subset)
        rna = rna[rna_cell_subset]
        rna.obs_names = rename_duplicates(rna.obs_names)
    sc.pp.filter_cells(rna, min_genes=200)
    sc.pp.filter_genes(rna, min_cells=3)
    
    label_idx_mapping = {}
    unique_labels = np.unique(rna.obs['celltype'])
    for i, name in enumerate(unique_labels):
        label_idx_mapping[name] = i
    print(label_idx_mapping)
    rna_label_int = rna.obs['celltype'].replace(label_idx_mapping)
    rna.obs['label'] = rna_label_int

    # Read in ATAC data
    atac_counts = csr_matrix(mmread(path.join(atac_path, 'counts.mtx')))
    atac_gene_names = read_txt_np(path.join(atac_path, 'genes.txt'))
    atac_cell_names = read_txt_np(path.join(atac_path, 'cells.txt'))
    atac_label = read_txt_np(path.join(atac_path, 'annotations.txt'))
    if binarize:
        atac_counts = (atac_counts > 0).astype('int')

    atac = sc.AnnData(atac_counts)
    atac.var_names = atac_gene_names
    atac.obs_names = atac_cell_names
    atac.var_names.name = 'Gene'
    atac.obs_names.name = 'CellID'
    atac.obs = pd.DataFrame({'celltype': atac_label}, index=atac_cell_names)
    atac.obs['label'] = -1
    atac.var.modality = 'Peaks'
    del atac_counts, atac_gene_names, atac_cell_names
    if atac_cell_subset is not None:
        atac_cell_subset = read_txt_np(atac_cell_subset)
        atac = atac[atac_cell_subset]
        atac.obs_names = rename_duplicates(atac.obs_names)
    sc.pp.filter_cells(atac, min_genes=200)
    sc.pp.filter_genes(atac, min_cells=3)

    # Concatenate RNA and ATAC
    adata = ad.concat([rna, atac], label = 'batch', keys = ['RNA', 'ATAC'])
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch', subset=True)
    
    # Run scVI
    scvi.settings.seed = 0
    scvi.model.SCVI_semi.setup_anndata(adata, batch_key='batch', layer="counts", labels_key='label')
    model = scvi.model.SCVI_semi(adata, n_latent=20)
    model.train()

    # Calculate the corrected latent variables.
    latent = model.get_latent_representation()
    latent = pd.DataFrame(latent, index=adata.obs_names, columns=['D%d' % i for i in range(1,latent.shape[1]+1)])

    # KNN classifier
    X_train = latent[adata.obs['batch'] == 'RNA']
    y_train = adata[adata.obs['batch'] == 'RNA'].obs['celltype']
    X_test = latent[adata.obs['batch'] == 'ATAC']

    knn = KNeighborsClassifier(n_neighbors=30).fit(X_train, y_train)
    tmp = y_train.unique().tolist()
    tmp.sort()

    if save:
        # prob matrix
        prob = pd.DataFrame(knn.predict_proba(X_test), index=X_test.index, columns=tmp)
        prob.to_csv(result_folder+'/prob.csv')
        # predicted label
        pred = prob.idxmax(axis=1)
        pred.to_csv(result_folder+'/pred.csv')
        # latent emb
        latent.to_csv(result_folder+'/latent.csv')

def main(): 
    
    start = time.time()

    run_scvi(rna_path, atac_path, result_folder, rna_cell_subset, atac_cell_subset, rna_new_annot)

    end = time.time()
    print('Running time: %.2f sec' % (end-start))
    with open(path.join(result_folder, 'time.txt'), "w") as file:
        file.write(str(end-start))
    
peak_mem_usage = memory_usage(main, max_iterations=1, max_usage=True)
print('Peak memory usage: %.2f MB' % peak_mem_usage)
with open(path.join(result_folder, 'memory.txt'), "w") as file:
    file.write(str(peak_mem_usage))
