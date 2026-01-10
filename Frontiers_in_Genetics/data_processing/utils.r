library(Seurat)
# library(SeuratDisk)
library(Signac)
# library(EnsDb.Hsapiens.v86)
# library(anndata)
library(dplyr)
library(Matrix)
# library(reticulate)
# np <- import("numpy")
library(genefu)

print.names.match <- function(name1, name2){
    print("Length1, Length2, Length intersection, Identical:")
    print(paste(length(name1), length(name2), length(intersect(name1, name2)), identical(name1, name2)))
}

convert.gene.activtity <- function(counts, frag.path, annotations, genome='hg38'){
    # Change style to UCSC
    # annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
    genome(annotations) <- genome
    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c("-", "-"),
        genome = genome,
        fragments = frag.path,
        min.cells = 10,
        min.features = 200,
        annotation = annotations
    )

    obj.atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
    gene.activities <- GeneActivity(obj.atac)
    return(gene.activities)
}

map.peaks <- function(peaks, frag.path){
    fragments <- CreateFragmentObject(
                  path = frag.path,
                  cells = NULL,
                  validate.fragments = FALSE
                )
    counts <- FeatureMatrix(
                fragments = fragments,
                features = peaks
              )
    return(counts)
}

write.counts.transpose <- function(counts, path, anno=NULL, mod='RNA'){
    if (!dir.exists(path)){
        dir.create(path, recursive = TRUE)
    }
    writeMM(t(counts), file = paste0(path, '/counts.mtx'))
    if (mod == 'ATAC'){
        new.peaks <- c()
        for (p in rownames(counts)){
            new.peaks <- c(new.peaks, gsub(':', '-', p))
        }
        write(new.peaks, paste0(path, '/genes.txt'))
    }else{
        write(rownames(counts), paste0(path, '/genes.txt'))
    }
    write(colnames(counts), paste0(path, '/cells.txt'))
    if (!is.null(anno)){
        write(anno, paste0(path, '/annotations.txt'))
    }       
}

tfidf <- function(mtx, method = 1, scale.factor = 1e4) {
  
  require(Matrix)
  
  npeaks <- colSums(mtx)
  if (any(npeaks == 0)) {
    warning("Some cells contain 0 total counts")
  }
  tf <- tcrossprod(mtx, y = Diagonal(x=1/npeaks))
  rsums <- rowSums(mtx)
  if (any(rsums == 0)) {
    warning("Some features contain 0 total counts")
  }
  idf <- ncol(mtx) / rsums
  if (method == 2) {
    idf <- log(1 + idf)
  } else if (method == 3) {
    tf <- log1p(tf * scale.factor)
    idf <- log(1 + idf)
  }
  mtx.tfidf <- Diagonal(n = length(idf), x = idf) %*% tf
  if (method == 1) {
    mtx.tfidf <- log1p(mtx.tfidf * scale.factor)
  }
  colnames(mtx.tfidf) <- colnames(mtx)
  rownames(mtx.tfidf) <- rownames(mtx)
  # set NA values to 0
  mtx.tfidf[is.na(mtx.tfidf)] <- 0
  return(mtx.tfidf)
}

rename.rna.genes <- function(gene.names){
    gene.names <- str_to_title(gene.names)
    idx <- grep('^.*rik$', gene.names)
    for (i in idx){
        gene <- gene.names[i]
        gene.names[i] <- paste0(toupper(substr(gene, 1, nchar(gene)-3)), 'Rik')
    }
    idx <- grep('^.*rik[0-9]$', gene.names)
    for (i in idx){
        gene <- gene.names[i]
        gene.names[i] <- paste0(toupper(substr(gene, 1, nchar(gene)-4)), 'Rik', substr(gene, nchar(gene), nchar(gene)))
    }
    idx <- grep('^[a-zA-Z]{2}[0-9]{6}$', gene.names)
    for (i in idx){
        gene <- gene.names[i]
        gene.names[i] <- toupper(gene)
    }
    idx <- grep('^[a-zA-Z]{2}[0-9]{6}.[0-9]$', gene.names)
    for (i in idx){
        gene <- gene.names[i]
        gene.names[i] <- toupper(gene)
    }
    return(gene.names)
}

set.names <- function(counts, rowname, colname, allow.duplicate = FALSE) {
    dims <- dim(counts)
    if (dims[1] != length(rowname)) {
        stop('Length of rownames not match!')
    }
    if (dims[2] != length(colname)) {
        stop('Length of colnames not match!')
    }
    if (!allow.duplicate) {
        if (any(duplicated(rowname))) {
            stop('Duplicated rownames!')
        }
        if (any(duplicated(colname))) {
            stop('Duplicated colnames!')
        }
    }else{
        if (any(duplicated(rowname))) {
            print('Duplicated rownames!')
        }
        if (any(duplicated(colname))) {
            print('Duplicated colnames! Renaming...')
            # Use the 'duplicated' function to find duplicates
            duplicates <- duplicated(colname)
            # Create a counter for each unique string
            counter <- rep(1, length(colname))
            # Loop through the duplicates and rename them
            for (i in seq_along(duplicates)) {
              if (duplicates[i]) {
                counter[i] <- counter[i] + 1
                colname[i] <- paste0(colname[i], ".", counter[i])
              }
            }
        }
    }
    rownames(counts) <- rowname
    colnames(counts) <- colname
    counts
}
# DEPRECATED
#load.data <- function(data.path, modality) {
#    data.count <- readMM(file.path(data.path, modality, "counts.mtx"))
#    cell.names <- read.table(file.path(data.path, modality, "cells.txt"), sep = ',')
#    cell.names <- cell.names[,1]
#    # this should be done in the prepare data process.
#    #cell.names <- gsub("-1.*", "-1", cell.names)
#    gene.names <- read.table(file.path(data.path, modality, "genes.txt"), sep = ',')
#    gene.names <- gene.names[,1]
#    if (file.exists(file.path(data.path, modality, "annotations.txt"))) {
#        anno <- read.table(file.path(data.path, modality, "annotations.txt"), sep = ',')
#        anno <- anno[,1]
#        names(anno) <- cell.names 
#    } else {
#        anno <- NULL
#    }
#    # rownames(data.count) <- cell.names
#    # colnames(data.count) <- gene.names
#    data.count <- set.names(data.count, cell.names, gene.names)
#    data.count <- as(t(data.count), "dgCMatrix")
#    data.dat <- list(counts = data.count, annotations = anno)
#    data.dat
#}

load.data.folder <- function(data.path, subset.file = NULL, annotation.file=NULL, allow.duplicate = FALSE) {
    if (endsWith(data.path, 'mtx')){
        data.count <- readMM(data.path)
        data.path <- paste(head(strsplit(data.path, '/')[[1]], -1), collapse='/')
    }else{
        data.count <- readMM(file.path(data.path, "counts.mtx"))
    }
    # cell.names <- read.table(file.path(data.path, "cells.txt"), sep = ',')
    # cell.names <- cell.names[,1]
    cell.names <- readLines(file.path(data.path, "cells.txt"))
    # this should be done in the prepare data process.
    #cell.names <- gsub("-1.*", "-1", cell.names)
    # gene.names <- read.table(file.path(data.path, "genes.txt"), sep = ',')
    # gene.names <- gene.names[,1]
    gene.names <- readLines(file.path(data.path, 'genes.txt'))
    if (!is.null(annotation.file)) {
        anno <- readLines(annotation.file)
        names(anno) <- cell.names 
    } else {
        if (file.exists(file.path(data.path, "annotations.txt"))) {
            # anno <- read.table(file.path(data.path, "annotations.txt"), sep = ',')
            # anno <- anno[,1]
            anno <- readLines(file.path(data.path, "annotations.txt"))
            names(anno) <- cell.names 
        } else {
            anno <- NULL
        }
    }
    # rownames(data.count) <- cell.names
    # colnames(data.count) <- gene.names
    data.count <- set.names(data.count, cell.names, gene.names, allow.duplicate = allow.duplicate)
    data.count <- as(t(data.count), "dgCMatrix")
    data.dat <- list(counts = data.count, annotations = anno)
    if (!is.null(subset.file)) {
        subset.barcodes <- readLines(subset.file)
#         subset.barcodes <- unique(subset.barcodes)  # Yuge
#         subset.barcodes <- c(subset.barcodes, rep(subset.barcodes[1], 5))  # Yuge
        if (length(subset.barcodes) == 0) {
            stop('empty subset!')
        }
        if (!all(subset.barcodes %in% colnames(data.dat$counts))) {
            stop('unseen cell in subset barcodes!')
        }
        data.dat$counts <- data.dat$counts[, subset.barcodes]
        subset.barcodes.dedup <- rename.duplicate(subset.barcodes)$new.x  # Yuge
        colnames(data.dat$counts) <- subset.barcodes.dedup  # Yuge
        if (!is.null(data.dat$annotations)) {
            data.dat$annotations <- data.dat$annotations[subset.barcodes]
            names(data.dat$annotations) <- subset.barcodes.dedup  # Yuge
        }
    }
    data.dat
}

# DEPRECATED
# load.rna <- function(data.path) {
#     data <- load.data(data.path, "RNA")
#     rna.count <- data$counts
#     anno <- data$annotations
#     obj.rna <- CreateSeuratObject(counts = rna.count, assay="RNA")
#     obj.rna$seurat_annotations <- as.factor(anno)
#     obj.rna
# }

load.rna.folder <- function(data.path, subset.file=NULL, annotation.file=NULL, assay='RNA', min.cells=NULL, min.features=NULL) {
    data <- load.data.folder(data.path, subset.file, annotation.file)
#     rna.count <- data$counts
#     obj.rna <- CreateSeuratObject(counts = data$counts, assay=assay, min.cells=min.cells, min.features=min.features)
    obj.rna <- CreateSeuratObject(counts = data$counts, assay=assay)
    if (!is.null(data$annotations)) {
        names(data$annotations) <- colnames(obj.rna)
        obj.rna$seurat_annotations <- as.factor(data$annotations)
    }
    if (!is.null(min.cells)){
        obj.rna <- obj.rna[rowSums(obj.rna[[assay]]@counts > 0) >= min.cells,]
    }
    if (!is.null(min.features)){
        obj.rna <- subset(obj.rna, subset = nFeature_RNA >= min.features)
    }
    obj.rna
}

# DEPRECATED
# load.atac.gene.activity.count <- function(data.path) {
#     data <- load.data(data.path, "ATAC.gene.activities")
#     atac.count <- data$counts
#     atac.count
# }

load.atac.gene.activity.count.folder <- function(data.path, subset.file=NULL, annotation.file=NULL) {
    data <- load.data.folder(data.path, subset.file, annotation.file)
    atac.count <- data$counts
    atac.count
}

# DEPRECATED
# load.atac <- function(data.path) { data <- load.data(data.path, "ATAC")
#     atac.count <- data$counts
#     obj.atac <- CreateSeuratObject(counts = atac.count, assay="ATAC")
#     if (!is.null(data$annotation)) {
#         obj.atac$seurat_annotations <- as.factor(data$annotations)
#     }
#     obj.atac
# }
    
load.atac.folder <- function(data.path, subset.file = NULL, annotation.file=NULL, assay='ATAC', min.cells=NULL, min.features=NULL) {
    data <- load.data.folder(data.path, subset.file, annotation.file)
#     atac.count <- data$counts
#     obj.atac <- CreateSeuratObject(counts = atac.count, assay=assay, min.cells=min.cells, min.features=min.features)
    obj.atac <- CreateSeuratObject(counts = data$counts, assay=assay)
    if (!is.null(data$annotations)) {
        names(data$annotations) <- colnames(obj.atac)
        obj.atac$seurat_annotations <- as.factor(data$annotations)
    }
    if (!is.null(min.cells)){
        obj.atac <- obj.atac[rowSums(obj.atac[[assay]]@counts > 0) >= min.cells,]
    }
    if (!is.null(min.features)){
        obj.atac <- subset(obj.atac, subset = nFeature_RNA >= min.features)
    }
    obj.atac
}
# DEPRECATED
# load.atac.annotated <- function(data.path) {
#     # atac.count <- readMM(file.path(data.path, "ATAC", "count.mtx"))
#     # atac.meta <- np$load(file.path(data.path, "ATAC", "metadata.npz"), allow_pickle = TRUE)
#     # cell.names <- gsub("-1.*", "-1", atac.meta[["cells"]]) # to match that of the frag file.
#     # rownames(atac.count) <- cell.names
#     # colnames(atac.count) <- atac.meta[["features"]]
#     # atac.count <- t(atac.count)
#     atac.data <- load.data(data.path, "ATAC")
#     atac.count <- atac.data$counts
#     # Get gene annotations
#     annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#     # Change style to UCSC
#     annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
#     genome(annotations) <- "hg38"
#     frag.file <- file.path(data.path, "ATAC", "atac_fragments.tsv.gz")
#     chrom_assay <- CreateChromatinAssay(
#       counts = atac.count,
#       sep = c(":", "-"),
#       genome = 'hg38',
#     fragments = frag.file,
#       min.cells = 10,
#       annotation = annotations
#     )

#     obj.atac <- CreateSeuratObject(counts = chrom_assay, assay="ATAC")
#     obj.atac$seurat_annotations <- as.factor(atac.data$annotations)
#     obj.atac
# }

# load.atac.annotated.folder <- function(data.path, subset.file=NULL, annotation.file, min.cells=3, min.features=200) {
#     # atac.count <- readMM(file.path(data.path, "ATAC", "count.mtx"))
#     # atac.meta <- np$load(file.path(data.path, "ATAC", "metadata.npz"), allow_pickle = TRUE)
#     # cell.names <- gsub("-1.*", "-1", atac.meta[["cells"]]) # to match that of the frag file.
#     # rownames(atac.count) <- cell.names
#     # colnames(atac.count) <- atac.meta[["features"]]
#     # atac.count <- t(atac.count)
#     atac.data <- load.data.folder(data.path, subset.file, annotation.file=NULL)
#     atac.count <- atac.data$counts
#     # Get gene annotations
#     annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#     # Change style to UCSC
#     annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
#     genome(annotations) <- "hg38"
#     frag.file <- file.path(data.path, "atac_fragments.tsv.gz")
#     chrom_assay <- CreateChromatinAssay(
#       counts = atac.count,
#       sep = c(":", "-"),
#       genome = 'hg38',
#     fragments = frag.file,
#       min.cells = 10,
#       annotation = annotations
#     )

#     obj.atac <- CreateSeuratObject(counts = chrom_assay, assay="ATAC", min.cells=3, min.features=200)
#     obj.atac$seurat_annotations <- as.factor(atac.data$annotations)
#     obj.atac
# }
# load.multi <- function(data.path) {
#     obj.rna <- load.rna(file.path(data.path, "multiome"))
#     obj.atac <- load.atac(file.path(data.path, "multiome"))
#     obj.multi <- CreateSeuratObject(counts = obj.rna[["RNA"]])
#     # Get % of mitochondrial genes
#     obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")

#     # Add the ATAC assay to the multiome object
#     obj.multi[["ATAC"]] <- obj.atac[["ATAC"]]
#     # NO NEED TO ADD ANNOTATIONS TO THIS OBJECT
#     obj.multi
# }


load.multi.folder <- function(data.path, subset.file = NULL) {
    if (endsWith(data.path, 'mtx')){
        counts.name <- tail(strsplit(data.path, '/')[[1]], 1)
        data.path <- paste(head(strsplit(data.path, '/')[[1]], -1), collapse='/')
        obj.rna <- load.rna.folder(file.path(data.path, "RNA", counts.name), subset.file)
        obj.atac <- load.atac.folder(file.path(data.path, "ATAC", counts.name), subset.file)
    }else{
        obj.rna <- load.rna.folder(file.path(data.path, "RNA"), subset.file)
        obj.atac <- load.atac.folder(file.path(data.path, "ATAC"), subset.file)
    }

    obj.multi <- CreateSeuratObject(counts = obj.rna[["RNA"]])
    # Get % of mitochondrial genes
    obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")

    # Add the ATAC assay to the multiome object
    obj.multi[["ATAC"]] <- obj.atac[["ATAC"]]
    # NO NEED TO ADD ANNOTATIONS TO THIS OBJECT
    obj.multi
}
