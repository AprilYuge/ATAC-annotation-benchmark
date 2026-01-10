from memory_profiler import memory_usage
save = True
binarize = False

import time
from cobolt.utils import SingleData, MultiomicDataset
from cobolt.model import Cobolt
import os
from scipy.io import mmread
import pandas as pd
import numpy as np
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

def run_cobolt(rna_path, atac_path, omics_rna_path, omics_atac_path, result_folder, rna_cell_subset, atac_cell_subset, rna_new_annot):

    # Read in RNA data
    count = mmread(os.path.join(rna_path, "counts.mtx")).tocsr().astype(float)
    if binarize:
        count = (count > 0).astype('int')
    feature = pd.read_csv(os.path.join(rna_path, "genes.txt"), header=None, usecols=[0])[0].values.astype('str')
    rna_barcode = pd.read_csv(os.path.join(rna_path, "cells.txt"), header=None, usecols=[0])[0].values.astype('str')
    rna_label = read_txt_np(path.join(rna_path, 'annotations.txt'))
    rna_label = pd.DataFrame({'celltype': rna_label}, index=rna_barcode)
    rna = SingleData("GeneExpr", "RNA", feature, count, rna_barcode)
    if rna_cell_subset is not None:
        rna_cell_subset = read_txt_np(rna_cell_subset)
        rna_cell_subset = np.array(['RNA~%s' % i for i in rna_cell_subset])
        rna.filter_barcode(rna_cell_subset)
    rna.filter_features(min_count=0, min_cell=3, upper_quantile=0.99, lower_quantile=0.7)
    rna.filter_cells(min_count=0, min_feature=200)
    
    # Read in ATAC data
    count = mmread(os.path.join(atac_path, "counts.mtx")).tocsr().astype(float)
    if binarize:
        count = (count > 0).astype('int')
    feature = pd.read_csv(os.path.join(atac_path, "genes.txt"), header=None, usecols=[0])[0].values.astype('str')
    atac_barcode = pd.read_csv(os.path.join(atac_path, "cells.txt"), header=None, usecols=[0])[0].values.astype('str')
    atac_label = read_txt_np(path.join(atac_path, 'annotations.txt'))
    atac_label = pd.DataFrame({'celltype': atac_label}, index=atac_barcode)
    atac = SingleData("ChromAccess", "ATAC", feature, count, atac_barcode)
    if atac_cell_subset is not None:
        atac_cell_subset = read_txt_np(atac_cell_subset)
        atac_cell_subset = np.array(['ATAC~%s' % i for i in atac_cell_subset])
        atac.filter_barcode(atac_cell_subset)
    atac.filter_features(min_count=0, min_cell=3, upper_quantile=0.99, lower_quantile=0.7)
    atac.filter_cells(min_count=0, min_feature=200)

    # Read in RNA data in multiomics
    count = mmread(os.path.join(omics_rna_path, "counts.mtx")).tocsr().astype(float)
    if binarize:
        count = (count > 0).astype('int')
    feature = pd.read_csv(os.path.join(omics_rna_path, "genes.txt"), header=None, usecols=[0])[0].values.astype('str')
    barcode = pd.read_csv(os.path.join(omics_rna_path, "cells.txt"), header=None, usecols=[0])[0].values.astype('str')
    multirna = SingleData("GeneExpr", "multiomics", feature, count, barcode)
    multirna.filter_features(min_count=0, min_cell=0, upper_quantile=0.99, lower_quantile=0.7)

    # Read in ATAC data in multiomics
    count = mmread(os.path.join(omics_atac_path, "counts.mtx")).tocsr().astype(float)
    if binarize:
        count = (count > 0).astype('int')
    feature = pd.read_csv(os.path.join(omics_atac_path, "genes.txt"), header=None, usecols=[0])[0].values.astype('str')
    barcode = pd.read_csv(os.path.join(omics_atac_path, "cells.txt"), header=None, usecols=[0])[0].values.astype('str')
    multiatac = SingleData("ChromAccess", "multiomics", feature, count, barcode)
    multiatac.filter_features(min_count=0, min_cell=0, upper_quantile=0.99, lower_quantile=0.7)

    # Merge the above datasets
    multi_dt = MultiomicDataset.from_singledata(rna, atac, multirna, multiatac)

    # Training
    if 'Embryo' in result_folder:
        model = Cobolt(dataset=multi_dt, lr=0.0002, n_latent=10)
    elif 'results_binarize/MouseBrain_snap' in result_folder:
        model = Cobolt(dataset=multi_dt, lr=0.0002, n_latent=10)
    elif 'results_balanced_ct/kidney' in result_folder:
        model = Cobolt(dataset=multi_dt, lr=0.0002, n_latent=10)
    else:
        model = Cobolt(dataset=multi_dt, lr=0.001, n_latent=10)
    model.train()

    # Calculate the corrected latent variables.
    model.calc_all_latent()
    latent = model.get_all_latent()
    latent = pd.DataFrame(latent[0], index=latent[1], columns=['D%d' % i for i in range(1,11)])
    latent = latent.loc[list(rna.barcode) + list(atac.barcode)]
    latent.index = [i.split('~')[1] for i in latent.index]

    # KNN classifier
    X_train = latent.loc[latent.index.isin(rna_label.index)]
    y_train = rna_label.loc[X_train.index]['celltype']
    X_test = latent.loc[latent.index.isin(atac_label.index)]

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

    run_cobolt(rna_path, atac_path, omics_rna_path, omics_atac_path, result_folder, rna_cell_subset, atac_cell_subset, rna_new_annot)

    end = time.time()
    print('Running time: %.2f sec' % (end-start))
    with open(path.join(result_folder, 'time.txt'), "w") as file:
        file.write(str(end-start))
    
peak_mem_usage = memory_usage(main, max_iterations=1, max_usage=True)
print('Peak memory usage: %.2f MB' % peak_mem_usage)
with open(path.join(result_folder, 'memory.txt'), "w") as file:
    file.write(str(peak_mem_usage))
