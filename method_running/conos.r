args = commandArgs(trailingOnly=TRUE)
if (!(length(args) %in% c(3,5,6))) {
      stop("three arguments must be supplied (rna path, gene activtity path, result.folder)", call.=FALSE)
}
rna.path <- args[1]
# atac.path <- args[2]
activity.path <- args[2]
result.folder <- args[3]
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
source("/gpfs/gibbs/pi/zhao/yw599/Multiome/utils.r")
library(conos)
library(pagoda2)
library(Seurat)
library(parallel)
run.conos <- function(rna.path, activity.path, result.folder, subset.rna, subset.atac, annot.new.rna) {
    if (!dir.exists(result.folder)){
        dir.create(result.folder, recursive = TRUE)
    }
    # obj.rna <- load.rna(data.path)
    # obj.atac <- load.atac(data.path)
    # obj.rna <- RenameCells(obj.rna, new.names = paste0(Cells(obj.rna), "-rna"))
    # obj.atac <- RenameCells(obj.atac, new.names = paste0(Cells(obj.atac), "-atac"))
    rna.data <- load.data.folder(rna.path, subset.rna, annot.new.rna)
    rna.counts <- rna.data$counts
    gene.activities.data <- load.data.folder(activity.path, subset.atac)
    gene.activities <- gene.activities.data$counts
    colnames(rna.counts) <- paste0(colnames(rna.counts), '-rna')
    panel <- list(rna = rna.counts, atac = gene.activities)
    if (endsWith(result.folder, 'pagoda')){
        p2l <- lapply(panel,basicP2proc,n.odgenes=3e3,min.cells.per.gene=-1,nPcs=30,make.geneknn=FALSE,n.cores=1,get.tsne=FALSE)
    }else{
        p2l <- lapply(panel, basicSeuratProc, tsne=FALSE)
    }
    ## instantiate Conos object
    con <- Conos$new(p2l, n.cores=1)
    ## build joint graph
    con$buildGraph(k=15, k.self=5, k.self.weigh=0.01, ncomps=30, n.odgenes=5e3, space='PCA')
    rna.annotation <- rna.data$annotations 
    atac.annotation <- gene.activities.data$annotations
    names(rna.annotation) <- paste0(names(rna.annotation), '-rna')
    cellannot <- rna.annotation
#     new.label.info <- con$propagateLabels(labels = cellannot, method='solver')
    if (endsWith(result.folder, 'PBMC/conos_pagoda')){
        new.label.info <- con$propagateLabels(labels = cellannot, method='solver')
    }else{
        new.label.info <- con$propagateLabels(labels = cellannot, verbose=TRUE)
    }
    # save results
    if (endsWith(result.folder, 'pagoda')){
        cell.names <- rownames(p2l[[2]]$counts)
    }else{
        cell.names <- colnames(p2l[[2]])
        cell.names <- intersect(cell.names, rownames(new.label.info$label.distribution))
    }
    # prob matrix
    prob <- new.label.info$label.distribution[cell.names,]
    prob <- as.matrix(prob)
#     prob.cols <- colnames(prob)
#     for (i in 1:length(prob.cols)){
#         ct <- prob.cols[i]
#         if (startsWith(ct, 'X')){
#             prob.cols[i] <- substr(ct, 2, nchar(ct))
#         }
#     }
#     colnames(prob) <- prob.cols
    write.csv(prob, file.path(result.folder, "prob.csv"))
    # predicted label
    pred <- new.label.info$labels[cell.names]
    write.csv(pred, file.path(result.folder, "pred.csv"))
    # uncertainty
    uncertainty <- new.label.info$uncertainty[cell.names]
    write.csv(uncertainty, file.path(result.folder, 'uncertainty.csv'))
    # writeLines(dim(rna.counts), file.path(result.folder, 'rna.dim'))
}
run.conos(rna.path, activity.path, result.folder, subset.rna, subset.atac, annot.new.rna)
