library(peakRAM)
save <- TRUE
binarize <- FALSE

args = commandArgs(trailingOnly=TRUE)
if (!(length(args) %in% c(5, 7, 8, 9))) {
      stop("five arguments must be supplied (rna path, atac path, multiome path, result.folder)", call.=FALSE)
}
rna.path <- args[1]
atac.path <- args[2]
multiome.rna.path <- args[3]
multiome.atac.path <- args[4]
result.folder <- args[5]
if (length(args) >= 7){
    subset.rna <- args[6]
    if (args[6] == '!') { subset.rna <- NULL }
    subset.atac <- args[7]
    if (args[7] == '!') { subset.atac <- NULL }
} else {
    subset.rna <- NULL
    subset.atac <- NULL
}
if (length(args) >= 8){
    subset.bridge <- args[8]
    if (args[8] == '!') { subset.bridge <- NULL }
} else {
    subset.bridge <- NULL
}
if (length(args) == 9) {
    annot.new.rna <- args[9]
} else {
    annot.new.rna <- NULL
}
source("/gpfs/gibbs/pi/zhao/yw599/Multiome/utils.r")

run.bridge <- function(rna.path, atac.path, multiome.rna.path, multiome.atac.path, result.folder, subset.rna, subset.atac, subset.bridge, annot.new.rna, multi.subsample=NULL) {
    ##############################################################################
    # multi.subsample: the subsample ratio of cells in the multiome data. Default is
    #                  NULL, representing using the entire data.
    ##############################################################################
    if (!dir.exists(result.folder)){
        dir.create(result.folder, recursive = TRUE)
    }
    obj.rna <- load.rna.folder(rna.path, subset.rna, annot.new.rna, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
    obj.atac <- load.atac.folder(atac.path, subset.atac, assay='ATAC', min.cells=3, min.features=200, binarize=binarize)
#     RenameAssays(object = obj.atac, RNA = 'ATAC')
    obj.multi <- load.multi.folder(rna.path = multiome.rna.path, atac.path = multiome.atac.path, subset.file = subset.bridge, binarize=binarize)
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
#     obj.rna <- SCTransform(object = obj.rna) %>% RunPCA# %>% RunUMAP(dims = 1:50, return.model = TRUE) 
    obj.rna <- NormalizeData(obj.rna) %>% FindVariableFeatures %>% ScaleData %>% RunPCA
    # normalize multiome RNA
    DefaultAssay(obj.multi) <- "RNA"
#     obj.multi <- SCTransform(obj.multi, verbose = FALSE) 
    obj.multi <- NormalizeData(obj.multi) %>% FindVariableFeatures %>% ScaleData
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
#     DefaultAssay(obj.rna) <- "SCT"
    DefaultAssay(obj.rna) <- "RNA"
    # Subsample multiome data is needed
    if (!is.null(multi.subsample)){
        obj.multi <- obj.multi[, sample(ncol(obj.multi), as.integer(ncol(obj.multi)*multi.subsample))]
    }
    obj.rna.ext <- PrepareBridgeReference(reference = obj.rna,
                                          bridge = obj.multi, 
    #                                       reference.reduction = "spca",
                                          reference.reduction = "pca",
                                          reference.dims = dims.rna,
#                                           normalization.method = "SCT"
                                          normalization.method = "LogNormalize"
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
    if (save) {
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
}

out <- peakRAM(
    run.bridge(rna.path, atac.path, multiome.rna.path, multiome.atac.path, result.folder, subset.rna, subset.atac, subset.bridge, annot.new.rna)
)

mem <- out$Peak_RAM_Used_MiB
time <- out$Elapsed_Time_sec
print(paste0('Running time: ', time, ' sec'))
print(paste0('Peak memory usage: ', mem, ' MB'))
write(time, file.path(result.folder, "time.txt"))
write(mem, file.path(result.folder, "memory.txt"))
