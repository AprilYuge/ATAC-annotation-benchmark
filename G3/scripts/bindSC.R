library(peakRAM)
save <- TRUE
binarize <- FALSE

# necessary pacakges
.libPaths(c(.libPaths(), "/gpfs/gibbs/pi/zhao/yw599/conda_envs/myR/lib/R/library"))
# .libPaths(c(.libPaths(), "/gpfs/ysm/apps/software/R/4.1.0-foss-2020b/lib64/R/library")) this other way around doesn't work
library(bindSC)
library(Dict)
library(Seurat)
library(Signac)
source("/gpfs/gibbs/pi/zhao/yw599/Multiome/utils.r")

args = commandArgs(trailingOnly=TRUE)
print(paste0('Number of args is ', length(args)))

rna.path <- args[1]
atac.path <- args[2]
ga.path <- args[3]
result.folder <- args[4]
if (length(args) >= 5){
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
subset.ga <- subset.atac
npcs <- 30

if (!dir.exists(result.folder)){
    dir.create(result.folder, recursive = TRUE)
}

load_dir <- '/gpfs/gibbs/pi/zhao/yw599/Multiome/data_pp'

out <- peakRAM(
    function(){
        rna <- load.rna.folder(rna.path, subset.rna, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
        rna <- RenameCells(rna, new.names = paste0(Cells(rna), "-rna"))
        atac <- load.atac.folder(atac.path, subset.atac, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
        ga <- load.atac.folder(ga.path, subset.ga, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
        atac.cells <- intersect(colnames(atac), colnames(ga))
        atac <- atac[, atac.cells]
        ga <- ga[, atac.cells]

        if (!grepl('Spleen', rna.path)){
            rna <- NormalizeData(rna)
            ga <- NormalizeData(ga)
        }
        rna <- FindVariableFeatures(rna, nfeatures = 10000)
        ga <- FindVariableFeatures(ga, nfeatures = 10000)
        gene.use <- intersect(VariableFeatures(rna), VariableFeatures(ga))

        atac <- RunTFIDF(atac)
        atac <- FindTopFeatures(atac, min.cutoff = 'q0')
        atac <- RunSVD(atac)

        # X <- rna[["RNA"]][gene.use,]
        # Z0 <- ga[["RNA"]][gene.use,]

        combined <- merge(rna, ga)
        combined <- ScaleData(combined, features = gene.use)
        combined <- RunPCA(combined, features = gene.use)
        combined <- combined@reductions$pca@cell.embeddings
        X <- combined[colnames(rna),]
        Z0 <- combined[colnames(ga),]
        Y  <- atac@reductions$lsi@cell.embeddings
        x.clst <- rep("ct", nrow(X))
        y.clst <- rep("ct", nrow(Z0))

        res <- BiCCA( X = t(X) ,
                     Y = t(Y), 
                     Z0 =t(Z0), 
                     X.clst = x.clst,
                     Y.clst = y.clst,
                     alpha = 0.5, 
                     lambda = 0.5,
                     K = 20,
                     temp.path  = result.folder,
                     num.iteration = 50,
                     tolerance = 0.01,
                     save = FALSE,
                     parameter.optimize = FALSE,
                     block.size = 0)

        latent <- rbind(res$u, res$r)
        if (save){
            write.csv(latent, file.path(result.folder, "latent.csv"))
        }
    }
)

mem <- out$Peak_RAM_Used_MiB
time <- out$Elapsed_Time_sec
print(paste0('Running time: ', time, ' sec'))
print(paste0('Peak memory usage: ', mem, ' MB'))
write(time, file.path(result.folder, "time.txt"))
write(mem, file.path(result.folder, "memory.txt"))
