args = commandArgs(trailingOnly=TRUE)
if (!(length(args) %in% c(3,5,6))) {
      stop("One argument must be supplied (input file).n", call.=FALSE)
}
source("./utils.r")
scGCN_path <- "./scGCN"
rna.path <- args[1]
atac.path <- args[2]
tmp_path <- args[3]
if (length(args) > 3){
    subset.rna <- args[4]
    if (args[4] == '!') { subset.rna <- NULL }
    subset.atac <- args[5]
    if (args[5] == '!') { subset.atac <- NULL }
} else {
    subset.rna <- NULL
    subset.atac <- NULL
}
if (length(args) == 6) {
    annot.new.rna <- args[6]
} else {
    annot.new.rna <- NULL
}
#data.path <- "/gpfs/gibbs/pi/zhao/xs272/Multiomics/sc_benchmark/data/BMMC_processed_s3d7"
# obj.rna <- load.rna(data.path)
# atac.ga <- load.atac.gene.activity.count(data.path)
# rna_mat <- obj.rna[["RNA"]][, ]
# atac_mat <- atac.ga
rna.data <- load.data.folder(rna.path, subset.rna, annot.new.rna)
rna_mat <- rna.data$counts
atac.data <- load.data.folder(atac.path, subset.atac)
atac_mat <- atac.data$counts
common.genes <- intersect(rownames(rna_mat), rownames(atac_mat))
rna_mat <- rna_mat[common.genes, ]
atac_mat <- atac_mat[common.genes, ]
rna_mat <- as(rna_mat, "matrix")
atac_mat <- as(atac_mat, "matrix")
rna.label <- data.frame(rna.data$annotations)
atac.label <- data.frame(atac.data$annotations)
atac.label[,1] <- rna.label[1,1] # this is used to solve th issue that scGCN can not deal with ATAC-specific cell types
count.list <- list(rna_mat, atac_mat)
label.list <- list(rna.label, atac.label)
dir.create(file.path(tmp_path, 'input'), recursive = TRUE)
source(file.path(scGCN_path, 'data_preprocess_utility.R'))
# stopifnot(dir.exists(file.path(tmp_path, 'input', 'Data1.csv')))
# stopifnot(dir.exists(file.path(tmp_path, 'input', 'Data2.csv')))
# stopifnot(dir.exists(file.path(tmp_path, 'input', 'Label1.csv')))
# stopifnot(dir.exists(file.path(tmp_path, 'input', 'Label2.csv')))

save_processed_data(count.list,label.list,folder=tmp_path)
# writeLines(dim(rna_mat), tmp_path+'/rna.dim')
