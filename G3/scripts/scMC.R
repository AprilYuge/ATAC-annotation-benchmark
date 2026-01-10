library(peakRAM)
save <- TRUE
binarize <- FALSE

library(Dict)
library(scMC)
library(Seurat)
source("/gpfs/gibbs/pi/zhao/yw599/Multiome/utils.r")

args = commandArgs(trailingOnly=TRUE)
print(paste0('Number of args is ', length(args)))
rna.path <- args[1]
atac.path <- args[2]
result.folder <- args[3]
if (length(args) >= 4){
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

if (!dir.exists(result.folder)){
    dir.create(result.folder, recursive = TRUE)
}

load_dir <- '/gpfs/gibbs/pi/zhao/yw599/Multiome/data_pp'

out <- peakRAM(
    function(){
        obj.rna <- load.rna.folder(rna.path, subset.rna, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
        obj.rna <- RenameCells(obj.rna, new.names = paste0(Cells(obj.rna), "-rna"))
        obj.atac <- load.atac.folder(atac.path, subset.atac, assay='RNA', min.cells=3, min.features=200, binarize=binarize)

        object.list <- vector("list", 2)
        object.list[[1]] <- obj.rna
        object.list[[2]] <- obj.atac

        object.list <- lapply(X = object.list, FUN = function(x) {
            if (!grepl('Spleen', rna.path)){
                x <- NormalizeData(x, verbose = FALSE)
            }
            x <- FindVariableFeatures(x, verbose = FALSE)
            # perform scaling on the previously identified variable features
            x <- ScaleData(x, verbose = FALSE)
        })

        options(future.rng.onMisuse="ignore")
        combined <- RunscMC(object.list, similarity.cutoff = 0.1)

        latent = combined@reductions$scMC@cell.embeddings
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

