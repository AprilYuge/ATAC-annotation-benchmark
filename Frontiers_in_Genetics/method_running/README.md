# ATAC-annotation-benchmark (method running)

This reporsity contains all scripts for running each automated ATAC annotation method benchamrked in the manuscript “*Benchmarking automated cell type annotation tools for single-cell ATAC-seq data*”. 

### Package requirement

The versions of packages we used are shown below.

|        | Package       | Version   |
|--------|---------------|-----------|
| Python | numpy         | 1.18.1    |
|        | pandas        | 1.3.5     |
|        | tensorflow    | 2.3.0     |
|        | scanpy        | 1.8.2     |
|        | scipy         | 1.8.0     |
|        | scikit-learn  | 1.0.2     |
|        | scJoint       | see below |
|        | scGCN         | see below |
|        | seaborn       | 0.11.2    |
|        | matplotlib    | 3.5.1     |
| R      | Seurat        | 4.1.0     |
|        | Conos         | 1.4.6     |
|        | Signac        | 1.6.0     |
|        | Pagoda2       | 1.0.10    | 

### Code references

To make the implementation of [scGCN](https://github.com/QSong-github/scGCN) and [scJoint](https://github.com/SydneyBioX/scJoint) easier, we modified their source code by adding some helper functions. Moreover for sGCN, since the original code rely on older versions of tensorflow and numpy, we also modified their code to make it compatible with more updated tensorflow and numpy. The modified versions of scGCN and scJoint are under the folders named as `scGCN` and `scJoint` in this directory, respectively.