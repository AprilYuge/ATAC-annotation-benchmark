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
omics_rna_path = str(sys.argv[3])
omics_atac_path = str(sys.argv[4])
result_folder = str(sys.argv[5])
if not os.path.exists(result_folder):
    os.makedirs(result_folder)
if len(sys.argv) >= 8:
    rna_cell_subset = str(sys.argv[6])
    if rna_cell_subset == '!':
        rna_cell_subset = None
    atac_cell_subset = str(sys.argv[7])
    if atac_cell_subset == '!':
        atac_cell_subset = None
else:
    rna_cell_subset = None
    atac_cell_subset = None
if len(sys.argv) >= 9:
    rna_new_annot = sys.argv[8]
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

def run_multivi(rna_path, atac_path, omics_rna_path, omics_atac_path, result_folder, rna_cell_subset, atac_cell_subset, rna_new_annot):

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

    # Read in RNA data in multiomics
    rna_counts = csr_matrix(mmread(path.join(omics_rna_path, 'counts.mtx')))
    rna_gene_names = read_txt_np(path.join(omics_rna_path, 'genes.txt'))
    rna_cell_names = read_txt_np(path.join(omics_rna_path, 'cells.txt'))
#     rna_label = read_txt_np(path.join(omics_rna_path, 'annotations.txt'))
    if binarize:
        rna_counts = (rna_counts > 0).astype('int')

    multirna = sc.AnnData(rna_counts)
    multirna.var_names = rna_gene_names
    multirna.obs_names = rna_cell_names
    multirna.var_names.name = 'Gene'
    multirna.obs_names.name = 'CellID'
    multirna.obs = pd.DataFrame({'label': -1}, index=rna_cell_names)
#     multirna.obs = pd.DataFrame({'celltype': rna_label}, index=rna_cell_names)
    del rna_counts, rna_gene_names, rna_cell_names

    # Read in ATAC data in multiomics
    atac_counts = csr_matrix(mmread(path.join(omics_atac_path, 'counts.mtx')))
    atac_gene_names = read_txt_np(path.join(omics_atac_path, 'genes.txt'))
    atac_cell_names = read_txt_np(path.join(omics_atac_path, 'cells.txt'))
#     atac_label = read_txt_np(path.join(omics_atac_path, 'annotations.txt'))
    if binarize:
        atac_counts = (atac_counts > 0).astype('int')

    multiatac = sc.AnnData(atac_counts)
    multiatac.var_names = atac_gene_names
    multiatac.obs_names = atac_cell_names
    multiatac.var_names.name = 'Gene'
    multiatac.obs_names.name = 'CellID'
#     multiatac.obs = pd.DataFrame({'celltype': atac_label}, index=atac_cell_names)
    del atac_counts, atac_gene_names, atac_cell_names

    # Concatenate RNA and ATAC in multiomics together
    multi = ad.concat([multirna, multiatac], axis=1, label = 'modality', keys = ['Gene Expression', 'Peaks'])
    multi.obs = multirna.obs
    del multirna, multiatac

    # We can now use the organizing method from scvi to concatenate these anndata
    adata_mvi = scvi.data.organize_multiome_anndatas(multi, rna, atac)

    # MultiVI requires the features to be ordered, such that genes appear before genomic regions
    adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()

    # Filter features to remove those that appear in fewer than 1% of the cells
    sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))

    # Setup and Training MultiVI
    scvi.model.MULTIVI_semi.setup_anndata(adata_mvi, batch_key="modality", labels_key='label')
    mvi = scvi.model.MULTIVI_semi(
            adata_mvi,
            n_genes=(adata_mvi.var["modality"] == "Gene Expression").sum(),
            n_regions=(adata_mvi.var["modality"] == "Peaks").sum(),
            n_latent=20
    )
    mvi.train()


    # Calculate the corrected latent variables.
    latent = mvi.get_latent_representation()
    latent = pd.DataFrame(latent, index=adata_mvi.obs_names, columns=['D%d' % i for i in range(1,latent.shape[1]+1)])
    latent = latent.loc[list(adata_mvi[adata_mvi.obs.modality == 'expression'].obs_names) + 
                        list(adata_mvi[adata_mvi.obs.modality == 'accessibility'].obs_names)]
    latent.index = ['_'.join(i.split('_')[:-1]) for i in latent.index]

    # KNN classifier
    X_train = latent.loc[latent.index.isin(rna.obs_names)]
    y_train = rna[X_train.index].obs['celltype']
    X_test = latent.loc[latent.index.isin(atac.obs_names)]

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

    run_multivi(rna_path, atac_path, omics_rna_path, omics_atac_path, result_folder, rna_cell_subset, atac_cell_subset, rna_new_annot)

    end = time.time()
    print('Running time: %.2f sec' % (end-start))
    with open(path.join(result_folder, 'time.txt'), "w") as file:
        file.write(str(end-start))
    
peak_mem_usage = memory_usage(main, max_iterations=1, max_usage=True)
print('Peak memory usage: %.2f MB' % peak_mem_usage)
with open(path.join(result_folder, 'memory.txt'), "w") as file:
    file.write(str(peak_mem_usage))
