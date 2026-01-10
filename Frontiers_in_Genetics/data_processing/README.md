# ATAC-annotation-benchmark (data processing)

This folder contains notebooks for data processing in the manuscript “*Benchmarking automated cell type annotation tools for single-cell ATAC-seq data*”. 

We collected data from five different tissues, including mouse lung, mouse brain, mouse kidney, human peripheral blood mononuclear cell (PBMC) and human bone marrow mononuclear cells (BMMC) to benchmark differnet methods for automated scATAC-seq label annotation. For mouse lung, scRNA-seq data from both 10x Chromium (droplet-based) and Smart-seq2 (FACS-based) were collected. Among all the methods, only Bridge integration required multimodal data where scATAC-seq and scRNA-seq were simultaneously measured. Therefore, we collected multimodal data for each tissue except for mouse lung (the SHARE-seq data for mouse lung were sequenced too shallowly to be used). For the mouse brain, both SHARE-seq and SNARE-seq data were used as the multimodal data to benchmark Bridge integration separately.

Scripts for data processing are aggregated by each tissue type and separated into three parts, including label unification between scATAC-seq and scRNA-seq datasets, other processing steps, and cell type similarity matrix calculation. For BMMC, other processing steps also incude code used for generating data under each experimental scenario. 

All the single-cell data used in this manuscript are publicly available. Detailed information of each data and their downloadable links can be found in `STable1.xlsx`.