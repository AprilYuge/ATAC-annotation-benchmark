import numpy as np
import pandas as pd
import random
from os import path
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, precision_score, recall_score

def read_txt_np(filename):
    with open(filename) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
        return np.array(lines)

def cal_entropy_enrichment(prob, RNA_path, RNA_subset = None):
    """Calculate scaled entropy and enrichment for ATAC-specific cell types
    
    Parameters:
    -------------
    prob: `pandas.core.frame.DataFrame`
        The predicted probability matrix of cells from cell types unique to scATAC-seq data. 
        Rows are cells with cell barcodes as Index and columns are cell types that are 
        contained in the reference (scRNA-seq) data.
    RNA_path: `string`
        The path to the annotations.txt file for scRNA-seq. This is to calculated the
        background composition of cell types.
    
    """
    from collections import Counter
    
    if RNA_path.endswith('h5ad'):
        import scanpy as sc
        adata = sc.read_h5ad(RNA_path)
        anno = adata.obs.celltype.tolist()
    elif RNA_path.endswith('txt'):
        anno = read_txt_np(RNA_path)
    else:
        anno = read_txt_np(path.join(RNA_path, 'annotations.txt'))
    if RNA_subset is not None:
        cells = np.loadtxt(path.join(RNA_path, 'cells.txt'), dtype=object)
        anno = pd.DataFrame(anno, index = cells, columns = ['annotation'])
        anno = anno.loc[RNA_subset]
        anno = list(anno.annotation)
    N = len(anno)
    tmp = Counter(anno)
    dic = {}
    for ct in prob.columns:
        dic[ct] = tmp[ct]
    N = sum(dic.values())
    dic = {k:dic[k]/N for k in dic.keys()}
    prob_rna = pd.DataFrame.from_dict(dic, orient='index')
    S = prob.divide(prob_rna.transpose().iloc[0], axis=1)
    S = S.divide(S.sum(1), axis=0).values
    entropy = np.nansum(-S*np.log2(S), axis=1).mean()/np.log2(S.shape[1])
    enrichment = np.max(S, axis=1).mean()
    return entropy, enrichment

def get_metrics(prob, true_label, sim, method, RNA_path, atac_unique, RNA_subset = None):
    """Evaluate the accuracy of label transfer
    
    Parameters:
    -------------
    prob: `pandas.core.frame.DataFrame`
        The predicted probability matrix. Rows are cells with cell barcodes as Index and 
        columns are cell types that are contained in the reference (scRNA-seq) data.
    true_label: `pandas.core.frame.DataFrame`
        The true labels of cells in scATAC-seq data. Indexes are cell barcodes and
        the only column is called 'annotation'.
    sim: `pandas.core.frame.DataFrame`
        The similarity matrix showing the relationship between each pair of cell types.
        Rows are cell types that are contained in scATAC-seq and columns are cell types
        that are contained in scRNA-seq. Values are between 0 and 1. For the pair that is
        composed of the same cell type, the value should be 1.
    method: `string`
        Name of method.
    RNA_path: `string`
        The path to the annotations.txt file for scRNA-seq. This is to calculated the
        background composition of cell types.
    atac_unique: `boolean`
        Whether or not to calculate metrics designed for unique cell types in ATAC.
    
    """
    
    # Check if cells are contained in the index of true_label
    assert set(prob.index) <= set(true_label.index), 'not all cells in scATAC-seq are in the row of true_label'
      
    # Make the order of cells same in prob and true_label
    true_label = true_label.loc[prob.index]
    
    # Check if cell types in true_label are contained in the index of sim
    assert set(true_label.annotation) <= set(sim.index), 'not all cell types in scATAC-seq are in the row of sim'
    
    # Check if cell types in column of prob are contained in the column of sim
    assert set(prob.columns) <= set(sim.columns), 'set of cell types in column of prob and sim are not the same'
    
    # Get commons ct
    common_ct = list(set(sim.index).intersection(sim.columns))
    true_label_common = true_label.loc[true_label.annotation.isin(common_ct)]
    prob_common = prob.loc[true_label_common.index]
    
    # Make the order of column names same in prob and sim
    sim = sim[prob.columns]
    sim_common = sim[prob_common.columns]
    
    # Calculate weighted acc using sim
    sim_extend = sim.loc[true_label.annotation]
    weighted_acc = (sim_extend.values * prob.values).sum(axis=1).mean()
    sim_extend_common = sim_common.loc[true_label_common.annotation]
    weighted_acc_common = (sim_extend_common.values * prob_common.values).sum(axis=1).mean()
    
    if atac_unique:
        # Calculate weighted acc using sim for unique cell types in atac only
        unique_ct = list(set(true_label.annotation).difference(prob.columns))
        mask = true_label.annotation.isin(unique_ct)
        if sum(mask) == 0:
            weighted_acc_unique = entropy_unique = enrichment_unique = f1_unique = np.nan
        else:
            # weighted acc
            weighted_acc_unique = (sim_extend.values[mask,:] * prob.values[mask,:]).sum(axis=1).mean()    
            # Scaled entropy & enrichment
            entropy_unique, enrichment_unique = cal_entropy_enrichment(prob = prob[mask], RNA_path = RNA_path, RNA_subset = RNA_subset)
            # F1 of entropy and enrichment (unknown -> higher entropy & lower enrichment)
            f1_unique = 2*entropy_unique*(1-enrichment_unique)/(entropy_unique+(1-enrichment_unique))
    
    # Calculate other metrics
    y_true = true_label.annotation.to_list()
    y_true_common = true_label_common.annotation.to_list()
    if method == 'random':
        y_pred = random.choices(prob.columns, prob.iloc[0], k=prob.shape[0])
        y_pred_common = random.choices(prob_common.columns, prob_common.iloc[0], k=prob_common.shape[0])
    else:
        y_pred = prob.idxmax(axis=1).to_list()
        y_pred_common = prob_common.idxmax(axis=1).to_list()
    
    if atac_unique:
        met = [accuracy_score(y_true, y_pred),
               weighted_acc,
               f1_score(y_true, y_pred, average='macro'),
               precision_score(y_true, y_pred, average='macro'),
               recall_score(y_true, y_pred, average='macro'),
               accuracy_score(y_true_common, y_pred_common),
               weighted_acc_common,
               f1_score(y_true_common, y_pred_common, average='macro'),
               precision_score(y_true_common, y_pred_common, average='macro'),
               recall_score(y_true_common, y_pred_common, average='macro'),
               weighted_acc_unique, 
               entropy_unique, enrichment_unique, f1_unique]
    else:
        met = [accuracy_score(y_true, y_pred),
               weighted_acc,
               f1_score(y_true, y_pred, average='macro'),
               precision_score(y_true, y_pred, average='macro'),
               recall_score(y_true, y_pred, average='macro'),
               accuracy_score(y_true_common, y_pred_common),
               weighted_acc_common,
               f1_score(y_true_common, y_pred_common, average='macro'),
               precision_score(y_true_common, y_pred_common, average='macro'),
               recall_score(y_true_common, y_pred_common, average='macro'),]
    return met

def aggregate_res_xp(params, atac_unique=True):
    """Aggregating evaluation metrics for a list of dataset*method
    
    Paramter
    --------------
    params: This is a list of lists. Each list contains probability matrix, true labels,
        similarity matrix, data name, method name, and path to RNA data. 
        (e.g. [prob, true_label, sim, 'brain-share', 'bridge', RNA_path])
    
    """
    
    results = []
    if atac_unique:
        names = ['overall_accuracy', 'weighted_accuracy',
                 'weighted_accuracy_for_unique_ct_in_atac',
                 'entropy_for_unique_ct_in_atac',
                 'enrichment_for_unique_ct_in_atac',
                 'f1_for_unique_ct_in_atac',
                 'f1 macro', 'precision macro', 'recall macro',
                 'data', 'seed', 'method']
    else:
        names = ['overall_accuracy', 'weighted_accuracy',
                 'f1 macro', 'precision macro', 'recall macro',
                 'data', 'seed', 'method']
    
    for prob, true_label, sim, data, seed, method, RNA_path in params:
        metrics = get_metrics(prob, true_label, sim, method, RNA_path, atac_unique)
        results.append(metrics + [data, seed, method])
        res = pd.DataFrame(results, columns=names)
        
    return res    

def aggregate_res(params, atac_unique=True):
    """Aggregating evaluation metrics for a list of dataset*method
    
    Paramter
    --------------
    params: This is a list of lists. Each list contains probability matrix, true labels,
        similarity matrix, data name, method name, and path to RNA data. 
        (e.g. [prob, true_label, sim, 'brain-share', 'bridge', RNA_path])
    
    """
    
    results = []
    if atac_unique:
        names = ['overall_accuracy', 'weighted_accuracy',
                 'f1_macro', 'precision_macro', 'recall_macro',
                 'overall_accuracy_common', 'weighted_accuracy_common',
                 'f1_macro_common', 'precision_macro_common', 'recall_macro_common',
                 'weighted_accuracy_for_unique_ct_in_atac',
                 'entropy_for_unique_ct_in_atac',
                 'enrichment_for_unique_ct_in_atac',
                 'f1_for_unique_ct_in_atac',
                 'data', 'method']
    else:
        names = ['overall_accuracy', 'weighted_accuracy',
                 'f1_macro', 'precision_macro', 'recall_macro',
                 'overall_accuracy_common', 'weighted_accuracy_common',
                 'f1_macro_common', 'precision_macro_common', 'recall_macro_common',
                 'data', 'method']
    
    for prob, true_label, sim, data, method, RNA_path, RNA_subset in params:
        print(data)
        metrics = get_metrics(prob, true_label, sim, method, RNA_path, atac_unique, RNA_subset)
        results.append(metrics + [data, method])
        res = pd.DataFrame(results, columns=names)
        
    return res    
    
    
