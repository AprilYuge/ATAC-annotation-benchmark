# ATAC-annotation-benchmark (method running)

This reporsity contains scripts for reproducing results in the manuscript “*Benchmarking automated cell type annotation tools for single-cell ATAC-seq data*”. Scripts are divided into three parts. The first part contains scripts used for data clearning, formatting, peak set mapping, gene activity matrix calculation and similarity matrix calculation in the [data_processing](https://github.com/AprilYuge/ATAC-annotation-benchmark/tree/main/data_processing) folder. The second part contains all scripts used for running the five methods we benchmarked (Bridge integration, scGCN, Seurat v3, Conos and scJoint) and can be found in folder [method_running](https://github.com/AprilYuge/ATAC-annotation-benchmark/tree/main/method_running). The third part includes scripts and notebooks used for metric calculation in the [metric_calculation](https://github.com/AprilYuge/ATAC-annotation-benchmark/tree/main/metric_calculation) folder.

More details about our benchmaking study can be found in our [manuscript](https://www.biorxiv.org/content/10.1101/2021.11.08.467781v3.full).

### Package requirement

ResPAN is implemented in Python and based on the framework of PyTorch. Before downloading and installing ResPAN, some packages need to be installed first. These required packages along with their versions used in our manuscript are listed below.

| Package    | Version      |
|------------|--------------|
| numpy      | 1.18.1       |
| pandas     | 1.3.5        |
| scipy      | 1.8.0        |
| scanpy     | 1.8.2        |
| pytorch    | 1.10.2+cu113 |

### Code references

For the implementation of ResPAN, we referred to [WGAN-GP](https://github.com/Zeleni9/pytorch-wgan) for the structure of Generative Adversarial Network and [iMAP](https://github.com/Svvord/iMAP) for the random walk mutual nearest neighbor method. Many thanks to their open-source treasure.