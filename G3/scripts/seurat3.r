library(peakRAM)
save <- TRUE
binarize <- FALSE

args = commandArgs(trailingOnly=TRUE)
if (!(length(args) %in% c(4, 6, 7))) {
      stop("three arguments must be supplied (rna path, atac path, gene activtity path, result.folder)", call.=FALSE)
}
rna.path <- args[1]
atac.path <- args[2]
activity.path <- args[3]
result.folder <- args[4]
if (length(args) >= 6) {
    subset.rna <- args[5]
    if (args[5] == '!') { subset.rna <- NULL }
    subset.atac <- args[6]
    if (args[6] == '!') { subset.atac <- NULL }
} else {
    subset.rna <- NULL
    subset.atac <- NULL
}
if (length(args) == 7) {
    annot.new.rna <- args[7]
} else {
    annot.new.rna <- NULL
}
print(args)
source("/gpfs/gibbs/pi/zhao/yw599/Multiome/utils.r")

run.seurat3 <- function(rna.path, atac.path, activity.path, result.folder, subset.rna, subset.atac, annot.new.rna) {
    if (!dir.exists(result.folder)){
        dir.create(result.folder, recursive = TRUE)
    }
    obj.rna <- load.rna.folder(rna.path, subset.rna, annot.new.rna, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
    obj.atac <- load.atac.folder(atac.path, subset.atac, assay='RNA', min.cells=3, min.features=200, binarize=binarize)

    # Perform standard analysis of each modality independently RNA analysis
    if (!grepl('Spleen', rna.path)){
        obj.rna <- NormalizeData(obj.rna)
    }
    obj.rna <- FindVariableFeatures(obj.rna)
    obj.rna <- ScaleData(obj.rna)
    obj.rna <- RunPCA(obj.rna)
    #obj.rna <- RunUMAP(obj.rna, dims = 1:30)

    # ATAC analysis add gene annotation information
 #   annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
 #   seqlevelsStyle(annotations) <- "UCSC"
 #   genome(annotations) <- "hg38"
 #   Annotation(obj.atac) <- annotations

    # We exclude the first dimension as this is typically correlated with sequencing depth
    obj.atac <- RunTFIDF(obj.atac)
    obj.atac <- FindTopFeatures(obj.atac, min.cutoff = "q0")
    obj.atac <- RunSVD(obj.atac)
    # obj.atac <- RunUMAP(obj.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

    # quantify gene activity
    #gene.activities <- GeneActivity(obj.atac, features = VariableFeatures(obj.rna))
    gene.activities <- load.atac.gene.activity.count.folder(activity.path, subset.atac, binarize=binarize)    
    # add gene activities as a new assay
    ga.assay <- CreateAssayObject(counts = gene.activities)
    ga.assay.sub <- subset(ga.assay, features=VariableFeatures(obj.rna), cells=colnames(obj.atac))
    obj.atac[["ACTIVITY"]] <- ga.assay.sub

    # normalize gene activities
    DefaultAssay(obj.atac) <- "ACTIVITY"
    if (!grepl('Spleen', rna.path)){
        obj.atac <- NormalizeData(obj.atac)
    }
    obj.atac <- ScaleData(obj.atac, features = rownames(obj.atac))

    # Identify anchors
    transfer.anchors <- FindTransferAnchors(reference = obj.rna, query = obj.atac, features = VariableFeatures(object = obj.rna),
        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

    celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = obj.rna$seurat_annotations,
        weight.reduction = obj.atac[["lsi"]], dims = 2:30)

    obj.atac <- AddMetaData(obj.atac, metadata = celltype.predictions)

    if (save) {
        # save results
        # prob matrix
        prob <- celltype.predictions[,2:(ncol(celltype.predictions)-1)]
        cell.types <- names(table(obj.rna$seurat_annotations)) # get original cell type names
        prob.cols <- colnames(prob) # get modified names in prob's columns
        for (ct in cell.types) {
            ct_ <- paste0("prediction.score.", gsub('\\s|\\+|\\/|-|[()]', '.', ct))
            # ct_ <- ct
            prob.cols[prob.cols == ct_] <- ct # recover original name
        }
        colnames(prob) <- prob.cols
        write.csv(prob, file.path(result.folder, "prob.csv"))
        # predicted label
        pred <- obj.atac$predicted.id
        write.csv(pred, file.path(result.folder, "pred.csv"))
        # writeLines(dim(obj.rna["RNA"]), file.path(result.folder, 'rna.dim'))
    }
}

out <- peakRAM(
    run.seurat3(rna.path, atac.path, activity.path, result.folder, subset.rna, subset.atac, annot.new.rna)
)

mem <- out$Peak_RAM_Used_MiB
time <- out$Elapsed_Time_sec
print(paste0('Running time: ', time, ' sec'))
print(paste0('Peak memory usage: ', mem, ' MB'))
write(time, file.path(result.folder, "time.txt"))
write(mem, file.path(result.folder, "memory.txt"))

