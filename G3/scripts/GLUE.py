from memory_profiler import memory_usage
save = True
binarize = False

import time
import scglue
import networkx as nx
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

gtf_dic = {'mm10': '/gpfs/gibbs/pi/zhao/yw599/Multiome/data/GTF/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz', 
           'mm9': '/gpfs/gibbs/pi/zhao/yw599/Multiome/data/GTF/gencode.vM1.annotation.gtf.gz', 
           'hg19': '/gpfs/gibbs/pi/zhao/yw599/Multiome/data/GTF/gencode.v43lift37.annotation.gtf.gz', 
           'hg38': '/gpfs/gibbs/pi/zhao/yw599/Multiome/data/GTF/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz'}

rna_path = str(sys.argv[1])
atac_path = str(sys.argv[2])
ref_genome = str(sys.argv[3])
result_folder = str(sys.argv[4])
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

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
    
    rna_counts = csr_matrix(mmread(path.join(rna_path, 'counts.mtx'))).astype(int)
    if binarize:
        rna_counts = (rna_counts > 0).astype('int')
    rna_gene_names = read_txt_np(path.join(rna_path, 'genes.txt'))
    rna_cell_names = read_txt_np(path.join(rna_path, 'cells.txt'))
    rna_label = read_txt_np(path.join(rna_path, 'annotations.txt'))

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

    atac_counts = csr_matrix(mmread(path.join(atac_path, 'counts.mtx'))).astype(int)
    if binarize:
        atac_counts = (atac_counts > 0).astype('int')
    atac_gene_names = read_txt_np(path.join(atac_path, 'genes.txt'))
    atac_cell_names = read_txt_np(path.join(atac_path, 'cells.txt'))
    atac_label = read_txt_np(path.join(atac_path, 'annotations.txt'))

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

    # Get gene annotations
    scglue.data.get_gene_annotation(
        rna, gtf=gtf_dic[ref_genome],
        gtf_by="gene_name"
    )
    rna = rna[:,~rna.var.chrom.isna()]
    print('Remained genes number is :', rna.shape[1])
    rna.var["chromStart"] = rna.var["chromStart"].astype(int)
    rna.var["chromEnd"] = rna.var["chromEnd"].astype(int)

    split = atac.var_names.str.split(r"[:-]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)

    # Data preprocessing
    sc.pp.filter_cells(rna, min_genes=200)
    sc.pp.filter_genes(rna, min_cells=3)
    rna.layers["counts"] = rna.X.copy()
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
    if 'Spleen' not in rna_path: 
        sc.pp.normalize_total(rna)
        sc.pp.log1p(rna)
    sc.pp.scale(rna)
    sc.tl.pca(rna, n_comps=100, svd_solver="auto")
    sc.pp.filter_cells(atac, min_genes=200)
    sc.pp.filter_genes(atac, min_cells=3)
    scglue.data.lsi(atac, n_components=100, n_iter=15)

    # Graph construction
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)

    # Configure data
    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca"
    )
    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=True,
        use_rep="X_lsi"
    )

    guidance_hvf = guidance.subgraph(chain(
        rna.var.query("highly_variable").index,
        atac.var.query("highly_variable").index
    )).copy()

    # Run GLUE
    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance_hvf,
        fit_kws={"directory": result_folder}
    )
    glue.save("%s/glue.dill" % result_folder)

    # Check integration consistency
    dx = scglue.models.integration_consistency(
        glue, {"rna": rna, "atac": atac}, guidance_hvf
    )
    print(dx)

    # KNN classifier
    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)
    X_train = rna.obsm['X_glue']
    y_train = rna.obs['celltype']
    X_test = atac.obsm['X_glue']

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
        latent_rna = pd.DataFrame(rna.obsm['X_glue'], index=rna.obs_names)
        latent_atac = pd.DataFrame(atac.obsm['X_glue'], index=atac.obs_names)
        latent = pd.concat([latent_rna, latent_atac])
        latent.to_csv(result_folder+'/latent.csv')
        # graph
        nx.write_graphml(guidance_hvf, "%s/guidance-hvf.graphml.gz" % result_folder)

    end = time.time()
    print('Running time: %.2f sec' % (end-start))
    with open(path.join(result_folder, 'time.txt'), "w") as file:
        file.write(str(end-start))
    
peak_mem_usage = memory_usage(main, max_iterations=1, max_usage=True)
print('Peak memory usage: %.2f MB' % peak_mem_usage)
with open(path.join(result_folder, 'memory.txt'), "w") as file:
    file.write(str(peak_mem_usage))
    


