library(peakRAM)
save <- TRUE
binarize <- FALSE

# necessary pacakges
library(rliger)
library(Dict)
library(scMC)
library(Seurat)
library(stringr)
source("/gpfs/gibbs/pi/zhao/yw599/Multiome/utils.r")

args = commandArgs(trailingOnly=TRUE)
print(paste0('Number of args is ', length(args)))
print(args)
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
print(result.folder)
npcs <- 30

if (!dir.exists(result.folder)){
    dir.create(result.folder, recursive = TRUE)
}

load_dir <- '/gpfs/gibbs/pi/zhao/yw599/Multiome/data_pp'

out <- peakRAM(
    function(){
        rna <- load.rna.folder(rna.path, subset.rna, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
        rna <- RenameCells(rna, new.names = paste0(Cells(rna), "-rna"))
        unshared_atac <- load.atac.folder(atac.path, subset.atac, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
        shared_atac <- load.atac.folder(ga.path, subset.ga, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
        cells_atac <- intersect(colnames(shared_atac), colnames(unshared_atac))
        shared_atac <- shared_atac[, cells_atac]
        unshared_atac <- unshared_atac[, cells_atac]

        chunk <- function(x,n){
            vect <- c(1:x)
            num <- ceiling(x/n)
            split(vect,rep(1:num,each=n,len=x))
        }

        liger <- createLiger(list(peaks = unshared_atac[["RNA"]]@counts))
        rm(unshared_atac)
        liger <- rliger::normalize(liger)
        norm <- liger@norm.data$peaks

        # Select the top 2,000 variable features
        se = CreateSeuratObject(norm)
        rm(norm)
        se <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
        top2000 <- head(VariableFeatures(se),2000)
        rm(se)
        # top2000_feats <-  norm[top2000,] 

        # liger <- selectGenes(liger)
        liger@var.genes <- top2000
        liger <- scaleNotCenter(liger)
        unshared_feats = liger@scale.data$peaks

        liger <- createLiger(list(rna = rna[["RNA"]]@counts, atac = shared_atac[["RNA"]]@counts))
        rm(rna)
        rm(shared_atac)
        liger <- rliger::normalize(liger)
        liger <- selectGenes(liger, var.thresh = 0.1, datasets.use = 1, unshared = TRUE, unshared.datasets = list(2), unshared.thresh = 0.2)
        liger <- scaleNotCenter(liger)

        peak_names <- colnames(unshared_feats)
        liger@var.unshared.features[[2]] <- peak_names
        liger@scale.unshared.data[[2]] <- t(unshared_feats)

        liger <- optimizeALS(liger, k=npcs, rand.seed=1)
        liger <- quantile_norm(liger, rand.seed=1)
        latent <- as.data.frame(liger@H.norm)
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
