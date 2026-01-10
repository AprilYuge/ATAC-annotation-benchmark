from memory_profiler import memory_usage
save = True
binarize = False

import time
import scglue
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
from sklearn.preprocessing import normalize

import math
# import sys
import torch
from torch.utils import data
import torch.optim
from torch.autograd import Variable
import torch.nn as nn
import torch.cuda

from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module
import torch.nn.functional as F

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

    # Get kernel inputs
    tmp = normalize(rna.obsm['X_pca'])
    K1 = tmp @ tmp.T
    tmp = normalize(atac.obsm['X_lsi'])
    K2 = tmp @ tmp.T

    # Run MMD-MA
    if torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device ="cpu"

    print("Running on",device)

    def compute_pairwise_distances(x, y): #function to calculate the pairwise distances
        if not len(x.size()) == len(y.size()) == 2:
            raise ValueError('Both inputs should be matrices.')

        if list(x.size())[1] != list(y.size())[1]:
            raise ValueError('The number of features should be the same.')

        diff =  (x.unsqueeze(2) - y.t())
        diff = (diff ** 2).sum(1)
        return diff.t()

    def gaussian_kernel_matrix(x, y, sigmas): #function to calculate Gaussian kernel
        beta = 1.0 / (2.0 * (sigmas.unsqueeze(1)))
        dist = compute_pairwise_distances(x, y)
        s = beta * (dist.contiguous()).view(1,-1)
        result =  ((-s).exp()).sum(0)
        return (result.contiguous()).view(dist.size())


    def stream_maximum_mean_discrepancy(x, y,  sigmas, kernel=gaussian_kernel_matrix): #This function has been implemented  to caculate MMD value for large number of samples (N>5,000)
        n_x = x.shape[0]
        n_y = y.shape[0]

        n_small = np.minimum(n_x, n_y)
        n = (n_small // 2) * 2

        cost = (kernel(x[:n:2], x[1:n:2], sigmas)  + kernel(y[:n:2], y[1:n:2], sigmas)
              - kernel(x[:n:2], y[1:n:2], sigmas) - kernel(x[1:n:2], y[:n:2], sigmas)).mean()
        if cost.data.item()<0:
            cost = torch.FloatTensor([0.0]).to(device)
        return cost

    def maximum_mean_discrepancy(x, y, sigmas, kernel=gaussian_kernel_matrix): #Function to calculate MMD value

        cost = (kernel(x, x, sigmas)).mean()
        cost += (kernel(y, y, sigmas)).mean()
        cost -= 2.0 * (kernel(x, y, sigmas)).mean()

        if cost.data.item()<0:
            cost = torch.FloatTensor([0.0]).to(device)

        return cost

    def calc_sigma(x1,x2): #Automatic sigma calculation 
        const = 8
        mat = np.concatenate((x1,x2))
        dist = []
        nsamp = mat.shape[0]
        for i in range(nsamp):
            euc_dist = np.sqrt(np.sum(np.square(np.subtract(mat[i,:], mat)), axis=1))
            dist.append(sorted(euc_dist)[1])

        sigma = np.square(const*np.median(dist))
        print("Calculated sigma:",sigma)
        return sigma

    class manifold_alignment(nn.Module): #MMD objective function

        def __init__(self, nfeat, num_k1, num_k2, seed=0):
            super(manifold_alignment, self).__init__()
            #Initializing alpha and beta 
            if torch.cuda.is_available():
                torch.cuda.manual_seed(seed)
            else:
                torch.manual_seed(seed)
            self.alpha = Parameter(torch.FloatTensor(num_k1, nfeat).uniform_(0.0,0.1).to(device))
            self.beta = Parameter(torch.FloatTensor(num_k2, nfeat).uniform_(0.0,0.1).to(device))


        def forward(self, k1, k2, ip, sigmas, lambda1, lambda2):


            if sigmas == 0: #If the user does not specify sigma values for the kernel calculation, they will be caclulated automatically
                x1 = (torch.matmul(k1,self.alpha)).detach().cpu().numpy()
                x2 = (torch.matmul(k2,self.beta)).detach().cpu().numpy()

                sigma = calc_sigma(x1,x2)
                sigmas = torch.FloatTensor([sigma]).to(device)

            mmd = maximum_mean_discrepancy(torch.matmul(k1,self.alpha),torch.matmul(k2,self.beta), sigmas)
            #mmd = stream_maximum_mean_discrepancy(torch.matmul(k1,self.alpha),torch.matmul(k2,self.beta), sigmas) #remove comment and comment the previous line if number of samples are large (N>5,000)

            penalty = lambda1 * ((((torch.matmul(self.alpha.t(),torch.matmul(k1,self.alpha))) - ip).norm(2))
                      + (((torch.matmul(self.beta.t(),torch.matmul(k2,self.beta))) - ip).norm(2)))

            distortion = lambda2 * ((((torch.matmul((torch.matmul(k1,self.alpha)),(torch.matmul(self.alpha.t(),k1.t()))))-k1).norm(2))
                      + (((torch.matmul((torch.matmul(k2,self.beta)),(torch.matmul(self.beta.t(),k2.t()))))-k2).norm(2)))

            return mmd, penalty, distortion, sigmas


    # #Functions to plot function values
    # def plot_data(filename,k,i,obj,mmd,pen,dist,nfeat,sigma,lambda1,lambda2):
    #     plt.xlabel('Iteration')
    #     plt.ylabel('log(Function value)')
    #     plt.title('nfeat:'+str(nfeat)+',seed:'+str(k)+', sigma:'+str(sigma)+', lambda1:'+str(lambda1)+', lambda2:'+str(lambda2))

    #     plt.plot(obj, 'k--', label='Objective')
    #     plt.plot(mmd, 'r--', label='MMD')
    #     plt.plot(pen, 'b--', label='Penalty')
    #     plt.plot(dist, 'g--', label='Distortion')
    #     if i == 10000:
    #         plt.legend(loc='upper right')
    #     plt.savefig(filename)
    #     plt.close()

    print("Loading data...")

    seed = 0

    k1_matrix = K1.astype(np.float32)
    k2_matrix = K2.astype(np.float32)
    del K1, K2

    print(k1_matrix.shape)
    print(k2_matrix.shape)

    print("Size of matrices...")
    num_k1 = k1_matrix.shape[0]                    #number of samples in dataset 1
    num_k2 = k2_matrix.shape[0]                    #number of samples in dataset 2

    nfeat = 5      # number features in joint embedding
    print("Number of dimensions of latent space...",nfeat)
    sigma = 0.0
    sigmas = torch.FloatTensor([sigma]).to(device)
    lambda_1 = 1e-5  # Lambda1 coefficient for penalty term
    lambda_2 = 1e-5  # Lambda2 coefficient for distortion term

    Ip = np.identity(nfeat).astype(np.float32)          #identity matrix of size nfeatxnfeat

    K1 = torch.from_numpy(k1_matrix).to(device)
    K2 = torch.from_numpy(k2_matrix).to(device)
    I_p = torch.from_numpy(Ip).to(device)

    obj_val=[]
    mmd_val=[]
    pen_val=[]
    dist_val=[]

    model = manifold_alignment(nfeat, num_k1, num_k2, seed)
    model = model.to(device)

    optimizer = torch.optim.Adam(model.parameters(),lr=0.0005,amsgrad=True)

    model.train()

    for i in range(10001): #Training takes place for 10,000 iterations

        optimizer.zero_grad()

        mmd, penalty, distortion, sigmas = model(K1, K2, I_p, sigmas, lambda_1, lambda_2)
        obj = mmd + penalty + distortion

        obj.backward()

        optimizer.step()

        obj_value = obj.data.item()
        mmd_value = mmd.data.item()
        pen_value = penalty.data.item()
        dist_value = distortion.data.item()

        if mmd_value > 0 : 
            obj_val.append(math.log(obj_value))
            mmd_val.append(math.log(mmd_value))
            pen_val.append(math.log(pen_value))
            dist_val.append(math.log(dist_value))

        if (i%200 == 0): #the weights can be saved every 200 iterations
            weights=[]

            for p in model.parameters():
                if p.requires_grad:
                    weights.append(p.data)

    #         plot_data(result_folder+"/Functions_"+str(seed)+".png",
    #                   seed,i,obj_val,mmd_val,pen_val,dist_val,nfeat,sigma,lambda_1,lambda_2)

    #         if (i==0 or i==10000): #This saves the weights at the beginning and end of the training
    #             np.savetxt(result_folder+"/alpha_hat_"+str(seed)+"_"+str(i)+".txt", weights[0].cpu().numpy())
    #             np.savetxt(result_folder+"/beta_hat_"+str(seed)+"_"+str(i)+".txt", weights[1].cpu().numpy())  

    # KNN classifier
    X_train = k1_matrix @ weights[0].cpu().numpy()
    X_test = k2_matrix @ weights[1].cpu().numpy()
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
