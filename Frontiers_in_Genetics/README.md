# Benchmarking automated cell type annotation tools for single-cell ATAC-seq data

This reporsity contains scripts and notebooks for reproducing results in the manuscript “*Benchmarking automated cell type annotation tools for single-cell ATAC-seq data*” published in *Frontiers in Genetics*. Scripts are divided into three parts. The first part contains notebooks used for data clearning, formatting, peak set mapping, gene activity matrix calculation and similarity matrix calculation in the [data_processing](https://github.com/AprilYuge/ATAC-annotation-benchmark/tree/main/data_processing) folder. The second part contains all scripts used for running the five methods we benchmarked (Bridge integration [1], scJoint [2], scGCN [3], Seurat v3 [4] and Conos [5]) and can be found in folder [method_running](https://github.com/AprilYuge/ATAC-annotation-benchmark/tree/main/method_running). The third part includes scripts and notebooks used for metric calculation in the [metric_calculation](https://github.com/AprilYuge/ATAC-annotation-benchmark/tree/main/metric_calculation) folder.

More details about our benchmaking study can be found in our [manuscript](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.1063233/full).

### Paper references

[1] Hao, Y., Stuart, T., Kowalski, M., Choudhary, S., Hoffman, P., Hartman, A., et al. (2022). Dictionary learning for integrative, multimodal, and scalable single-cell analysis. bioRxiv.

[2] Lin, Y., Wu, T.-Y., Wan, S., Yang, J.Y., Wong, W.H., and Wang, Y. (2022). scJoint integrates atlas-scale single-cell RNA-seq and ATAC-seq data with transfer learning. Nature Biotechnology 40(5), 703-710.

[3] Song, Q., Su, J., and Zhang, W. (2021). scGCN is a graph convolutional networks algorithm for knowledge transfer in single cell omics. Nature communications 12(1), 1-11.

[4] Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W.M., et al. (2019). Comprehensive integration of single-cell data. Cell 177(7), 1888-1902. e1821.

[5] Barkas, N., Petukhov, V., Nikolaeva, D., Lozinsky, Y., Demharter, S., Khodosevich, K., et al. (2019). Joint analysis of heterogeneous single-cell RNA-seq dataset collections. Nature methods 16(8), 695-698.
