# ATAC-annotation-benchmark (metric calculation)

This folder contains notebooks for metrics calculations in the manuscript “*Benchmarking automated cell type annotation tools for single-cell ATAC-seq data*”. 

In our study, we considered two set of metrics, which are overall metrics that assess using all ATAC cells and metrics that only assess the performance of different methods on ATAC-specific cell types. The first set of metrics contains overall accuracy, weighted accuracy and F1 (macro) of precision and recall. The second set of metrics contains weighted accuracy on ATAC-specific cell types and F1 of enrichment and entropy. Details of each metric can be found in our [manuscript](https://www.biorxiv.org/content/10.1101/2022.10.05.511014v1). 

In this folder, `metrics.py` contains functions used to calculate all metrics written in Python. The two notebooks contain code for metric calculation and visualizations on all tissues (`tissue_metrics.ipynb`) and experiments designed based on BMMC (`BMMC_metrics.ipynb`) separately. 
