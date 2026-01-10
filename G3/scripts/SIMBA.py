from memory_profiler import memory_usage
save = True
binarize = False

import time
import simba as si
import anndata as ad
import pandas as pd
import numpy as np
import os
from os import sys, path
from scipy.sparse import csc_matrix, csr_matrix
from scipy.io import mmwrite, mmread
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
            counts[string] = 1
        renamed.append(new_name)
    return renamed

rna_path = str(sys.argv[1])
atac_path = str(sys.argv[2])
ga_path = str(sys.argv[3])
result_folder = str(sys.argv[4])
if not os.path.exists(result_folder):
    os.makedirs(result_folder)
    
si.settings.set_workdir(result_folder)

load_dir = '/gpfs/gibbs/pi/zhao/yw599/Multiome/data_pp' 

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
    
    start = time.time()

    # Read in RNA data
    rna_counts = csr_matrix(mmread(path.join(rna_path, 'counts.mtx'))).astype('float')
    if binarize:
        rna_counts = (rna_counts > 0).astype('float')
    rna_gene_names = read_txt_np(path.join(rna_path, 'genes.txt'))
    rna_cell_names = read_txt_np(path.join(rna_path, 'cells.txt'))
    rna_label = read_txt_np(path.join(rna_path, 'annotations.txt'))

    rna = ad.AnnData(rna_counts)
    rna.var_names = rna_gene_names
    rna.obs_names = [b+'_rna' for b in rna_cell_names]
    rna.var_names.name = 'Gene'
    rna.obs_names.name = 'CellID'
    rna.obs = pd.DataFrame({'celltype': rna_label}, index=[b+'_rna' for b in rna_cell_names])
    rna.obs['dataset'] = 'RNA'
    del rna_counts, rna_gene_names, rna_cell_names
    
    # Read in ATAC data
    atac_counts = csr_matrix(mmread(path.join(atac_path, 'counts.mtx'))).astype('float')
    if binarize:
        atac_counts = (atac_counts > 0).astype('float')
    atac_gene_names = read_txt_np(path.join(atac_path, 'genes.txt'))
    atac_cell_names = read_txt_np(path.join(atac_path, 'cells.txt'))
    atac_label = read_txt_np(path.join(atac_path, 'annotations.txt'))

    atac = ad.AnnData(atac_counts)
    atac.var_names = atac_gene_names
    atac.obs_names = atac_cell_names
    atac.var_names.name = 'Gene'
    atac.obs_names.name = 'CellID'
    atac.obs = pd.DataFrame({'celltype': atac_label}, index=atac_cell_names)
    # atac.obs['label'] = -1
    atac.obs['dataset'] = 'ATAC'
    del atac_counts, atac_gene_names, atac_cell_names
    
    # Read in GA data
    ga_counts = csr_matrix(mmread(path.join(ga_path, 'counts.mtx'))).astype('float')
    if binarize:
        ga_counts = (ga_counts > 0).astype('float')
    ga_gene_names = read_txt_np(path.join(ga_path, 'genes.txt'))
    ga_cell_names = read_txt_np(path.join(ga_path, 'cells.txt'))
    ga_label = read_txt_np(path.join(ga_path, 'annotations.txt'))

    ga = ad.AnnData(ga_counts)
    ga.var_names = ga_gene_names
    ga.obs_names = ga_cell_names
    ga.var_names.name = 'Gene'
    ga.obs_names.name = 'CellID'
    ga.obs = pd.DataFrame({'celltype': ga_label}, index=ga_cell_names)
    del ga_counts, ga_gene_names, ga_cell_names
        
    # Data preprocessing
    si.pp.filter_genes(rna, min_n_cells=3)
    si.pp.filter_cells_rna(rna, min_n_genes=200)
    si.pp.normalize(rna, method='lib_size')
    si.pp.log_transform(rna)
    si.pp.select_variable_genes(rna, n_top_genes=4000)
    si.tl.discretize(rna, n_bins=5)
    
    si.pp.filter_peaks(atac, min_n_cells=3)
    si.pp.filter_cells_atac(atac, min_n_peaks=200)
    si.pp.pca(atac, n_components=50)
    si.pp.select_pcs_features(atac)

    si.pp.filter_genes(ga, min_n_cells=3)
    cells_use = set(atac.obs_names).intersection(ga.obs_names)
    ga = ga[list(cells_use)]
    atac = atac[list(cells_use)]
    si.pp.normalize(ga, method='lib_size')
    si.pp.log_transform(ga)
    CrnaCatac = si.tl.infer_edges(rna, ga, n_components=15, k=15)
    
    si.tl.trim_edges(CrnaCatac, cutoff=0.5)
    
    si.tl.gen_graph(list_CP=[atac],
                    list_CG=[rna],
                    list_CC=[CrnaCatac],
                    copy=False,
                    use_highly_variable=True,
                    use_top_pcs=True,
                    dirname='graph0')
    # Training
    # modify parameters
    dict_config = si.settings.pbg_params.copy()
    dict_config['workers'] = 1
    ## start training
    si.tl.pbg_train(pbg_params = dict_config, auto_wd=True, save_wd=True, output='model')
    #si.tl.pbg_train(auto_wd=True, save_wd=True, output='model')

    # Get embeddings
    dict_adata = si.read_embedding()
    atac_emb = dict_adata['C']  # embeddings for ATACseq cells
    rna_emb = dict_adata['C2']  # embeddings for RNAseq cells
    adata_all = si.tl.embed(adata_ref=rna_emb,
                            list_adata_query=[atac_emb],
                            use_precomputed=False)
    X_train = adata_all[adata_all.obs.id_dataset == 'ref'].X
    X_test = adata_all[adata_all.obs.id_dataset != 'ref'].X

    # KNN classifier
    y_train = rna.obs.loc[adata_all.obs_names[adata_all.obs.id_dataset == 'ref']]['celltype']

    knn = KNeighborsClassifier(n_neighbors=30).fit(X_train, y_train)
    tmp = y_train.unique().tolist()
    tmp.sort()
    
    if save:
        # prob matrix
        prob = pd.DataFrame(knn.predict_proba(X_test), index=adata_all.obs_names[adata_all.obs.id_dataset != 'ref'], columns=tmp)
        prob.to_csv(result_folder+'/prob.csv')
        # predicted label
        pred = prob.idxmax(axis=1)
        pred.to_csv(result_folder+'/pred.csv')
        # latent emb
        latent_rna = pd.DataFrame(X_train, index=[b[:-4] for b in adata_all.obs_names[adata_all.obs.id_dataset == 'ref']])
        latent_atac = pd.DataFrame(X_test, index=adata_all.obs_names[adata_all.obs.id_dataset != 'ref'])
        latent = pd.concat([latent_rna, latent_atac])
        latent.to_csv(result_folder+'/latent.csv')

    end = time.time()
    print('Running time: %.2f sec' % (end-start))
    with open(path.join(result_folder, 'time.txt'), "w") as file:
        file.write(str(end-start))
    
if __name__ == '__main__':
    peak_mem_usage = memory_usage(main, max_iterations=1, max_usage=True)
    print('Peak memory usage: %.2f MB' % peak_mem_usage)
    with open(path.join(result_folder, 'memory.txt'), "w") as file:
        file.write(str(peak_mem_usage))
