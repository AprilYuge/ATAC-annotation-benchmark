import TOSICA
import os
import sys
import time
import scanpy as sc
import numpy as np
import pandas as pd
import argparse
import torch
print(torch.__version__)
print(torch.cuda.get_device_capability(device=None),  torch.cuda.get_device_name(device=None))

parser = argparse.ArgumentParser()
parser.add_argument("--gene_num", type=int, default=2000, help='Number of highly variable genes.')
parser.add_argument("--epoch", type=int, default=5, help='Number of epochs.')
parser.add_argument("--seed", type=int, default=2021, help='Random seed.')
parser.add_argument("--label_name", type=str, default='celltype', help='Label name')
parser.add_argument("--novel_type", type=bool, default=False, help='Novel cell tpye exists or not.')
parser.add_argument("--unassign_thres", type=float, default=0.5, help='The confidence score threshold for novel cell type annotation.')
parser.add_argument("--model_path", type=str, default=None, help='Path of finetuned model.')
parser.add_argument("--load_dir", type=str, default='.', help='Directory to load data.')
parser.add_argument("--save_dir", type=str, default='.', help='Directory to save data.')
parser.add_argument("--gmt", type=str, default='human_gobp', help='GMT file to use.')

args = parser.parse_args()

SEED = args.seed
GENE_NUM = args.gene_num
EPOCHS = args.epoch
LABEL_NAME = args.label_name
LOAD_DIR = args.load_dir
GMT = args.gmt
SAVE_DIR = '%s/%s' % (args.save_dir, GMT)

if os.path.exists(SAVE_DIR) is False:
    os.makedirs(SAVE_DIR)
os.chdir(SAVE_DIR)

ref_adata = sc.read('%s/preprocessed_rna_data.h5ad' % LOAD_DIR)
# ref_adata = ref_adata[:,ref_adata.var_names]
sc.pp.highly_variable_genes(ref_adata, n_top_genes=GENE_NUM)
ref_adata = ref_adata[:, ref_adata.var.highly_variable]
print(ref_adata)
print(ref_adata.obs[LABEL_NAME].value_counts())

query_adata = sc.read('%s/preprocessed_ga_data.h5ad' % LOAD_DIR)
query_adata = query_adata[:,ref_adata.var_names]
print(query_adata)
print(query_adata.obs[LABEL_NAME].value_counts())

tissue = args.save_dir.split('/')[-1]
if tissue not in ['PBMC', 'BMMC']:
    ref_adata.var_names = [a[0]+a.lower()[1:] for a in ref_adata.var_names]
    query_adata.var_names = [a[0]+a.lower()[1:] for a in query_adata.var_names]

TOSICA.train(ref_adata, gmt_path=GMT, label_name=LABEL_NAME, epochs=EPOCHS)

today = time.strftime('%Y%m%d', time.localtime(time.time()))
model_weight_path = './weights%s/model-%s.pth' % (today, EPOCHS-1)
new_adata, prob = TOSICA.pre(query_adata, model_weight_path = model_weight_path, 
                             mask_path = '%s/mask.npy' % SAVE_DIR)
new_adata_ref, prob_ref = TOSICA.pre(ref_adata, model_weight_path = model_weight_path, 
                             mask_path = '%s/mask.npy' % SAVE_DIR)

# Save results for test
# adata
new_adata.obs = new_adata.obs.drop(columns=['n_genes'])
new_adata.write('adata_attn.h5ad')
# prob matrix
prob.to_csv('prob.csv')
# predicted label
pred = prob.idxmax(axis=1)
pred.to_csv('pred.csv')

# Save results for train
# adata
new_adata_ref.obs = new_adata_ref.obs.drop(columns=['n_genes'])
new_adata_ref.write('adata_attn_train.h5ad')
# prob matrix
prob_ref.to_csv('prob_train.csv')
# predicted label
pred_ref = prob_ref.idxmax(axis=1)
pred_ref.to_csv('pred_train.csv')