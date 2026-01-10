args = commandArgs(trailingOnly=TRUE)
if (!(length(args) %in% c(4, 7, 8))) {
      stop("five arguments must be supplied (rna path, atac path, multiome path, result.folder)", call.=FALSE)
}
rna.path <- args[1]
atac.path <- args[2]
multiome.path <- args[3]
result.folder <- args[4]
if (length(args) >= 7){
    subset.rna <- args[5]
    if (args[5] == '!') { subset.rna <- NULL }
    subset.atac <- args[6]
    if (args[6] == '!') { subset.atac <- NULL }
    subset.bridge <- args[7]
    if (args[7] == '!') { subset.bridge <- NULL }
} else {
    subset.rna <- NULL
    subset.atac <- NULL
    subset.bridge <- NULL
}
if (length(args) == 8) {
    annot.new.rna <- args[8]
} else {
    annot.new.rna <- NULL
}
source("/gpfs/gibbs/pi/zhao/yw599/Multiome/utils.r")
run.bridge <- function(rna.path, atac.path, multiome.path, result.folder, subset.rna, subset.atac, subset.bridge, annot.new.rna, multi.subsample=NULL) {
    ##############################################################################
    # multi.subsample: the subsample ratio of cells in the multiome data. Default is
    #                  NULL, representing using the entire data.
    ##############################################################################
    if (!dir.exists(result.folder)){
        dir.create(result.folder, recursive = TRUE)
    }
    obj.rna <- load.rna.folder(rna.path, subset.rna, annot.new.rna)
    obj.atac <- load.atac.folder(atac.path, subset.atac)
    obj.multi <- load.multi.folder(multiome.path, subset.bridge)
    # rename cells to avoid same names
    obj.rna <- RenameCells(obj.rna, new.names = paste0(Cells(obj.rna), "-rna"))
    obj.atac <- RenameCells(obj.atac, new.names = paste0(Cells(obj.atac), "-atac"))
    obj.multi <- RenameCells(obj.multi, new.names = paste0(Cells(obj.multi), "-multi"))
    if (any(duplicated(rownames(obj.multi[['RNA']])))) { stop('duplicate Rr') }
    if (any(duplicated(colnames(obj.multi[['RNA']])))) { stop('duplicate Rc') }
    if (any(duplicated(rownames(obj.multi[['ATAC']])))) { stop('duplicate Ar') }
    if (any(duplicated(colnames(obj.multi[['ATAC']])))) { stop('duplicate Ac') }

    common.genes <- intersect(rownames(obj.rna), rownames(obj.multi[['RNA']]))
    obj.rna <- obj.rna[common.genes, ]
    print(1)
    # obj.multi[["RNA"]] <- obj.multi[["RNA"]][common.genes, ]
    obj.multi[['RNA']] <- subset(obj.multi[['RNA']], features = common.genes)
    print(2)
    common.peaks <- intersect(rownames(obj.atac), rownames(obj.multi[['ATAC']]))
    obj.atac <- obj.atac[common.peaks, ]
    print(3)
    obj.multi[['ATAC']] <- subset(obj.multi[['ATAC']], features = common.peaks)
    print(4)
    # normalize rna
    obj.rna <- SCTransform(object = obj.rna) %>% RunPCA# %>% RunUMAP(dims = 1:50, return.model = TRUE) 
    # normalize multiome RNA
    DefaultAssay(obj.multi) <- "RNA"
    obj.multi <- SCTransform(obj.multi, verbose = FALSE) 
    # normalize multiome ATAC
    DefaultAssay(obj.multi) <- "ATAC"
    obj.multi <- RunTFIDF(obj.multi)
    obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")
    # normalize query
    obj.atac <- RunTFIDF(obj.atac)
    # Drop first dimension for ATAC reduction
    dims.atac <- 2:50
    dims.rna <- 1:50
    DefaultAssay(obj.multi) <-  "RNA"
    DefaultAssay(obj.rna) <- "SCT"
    # Subsample multiome data is needed
    if (!is.null(multi.subsample)){
        obj.multi <- obj.multi[, sample(ncol(obj.multi), as.integer(ncol(obj.multi)*multi.subsample))]
    }
    obj.rna.ext <- PrepareBridgeReference(reference = obj.rna,
                                          bridge = obj.multi, 
    #                                       reference.reduction = "spca",
                                          reference.reduction = "pca",
                                          reference.dims = dims.rna,
                                          normalization.method = "SCT"
    )
    # build anchor
    bridge.anchor <- FindBridgeTransferAnchors(extended.reference = obj.rna.ext, 
                                               query = obj.atac,
                                               reduction = "lsiproject",
                                               dims = dims.atac
    )
    # transfer label
    obj.atac <- MapQuery(anchorset = bridge.anchor, 
                         reference = obj.rna, 
                         query = obj.atac, 
                         refdata = list(
                             l1 = "seurat_annotations"
                                        )#
    #                      reduction.model = "wnn.umap" 
    )
    # save results
    # prob matrix
    prob <- t(obj.atac$prediction.score.l1[,])
    atac.cells <- rownames(prob)
    rownames(prob) <- substr(atac.cells, 1, nchar(atac.cells) - 5)
    cell.types <- names(table(obj.rna$seurat_annotations)) # get original cell type names
    prob.cols <- colnames(prob) # get modified names in prob's columns
    for (ct in cell.types) {
        ct_ <- gsub('_', '-', ct)
        # ct_ <- ct
        prob.cols[prob.cols == ct_] <- ct # recover original name
    }
    colnames(prob) <- prob.cols
    write.csv(prob, file.path(result.folder, "prob.csv"))
    # predicted label
    pred <- obj.atac$predicted.l1
    names(pred) <- substr(atac.cells, 1, nchar(atac.cells) - 5)
    write.csv(pred, file.path(result.folder, "pred.csv"))
    # writeLines(dim(obj.rna["RNA"]), file.path(result.folder, 'rna.dim'))
    # writeLines(dim(obj.multi["RNA"]), file.path(result.folder, 'multi.rna.dim'))
    # writeLines(dim(obj.multi["ATAC"]), file.path(result.folder, 'multi.rna.dim'))
}
run.bridge(rna.path, atac.path, multiome.path, result.folder, subset.rna, subset.atac, subset.bridge, annot.new.rna)

