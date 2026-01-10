library(peakRAM)
save <- TRUE
binarize <- FALSE

# necessary pacakges
library(rliger)
library(Dict)
library(Seurat)
source("/gpfs/gibbs/pi/zhao/yw599/Multiome/utils.r")

args = commandArgs(trailingOnly=TRUE)
print(paste0('Number of args is ', length(args)))
rna.path <- args[1]
atac.path <- args[2]
result.folder <- args[3]
datasets.use <- args[4]
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
npcs <- 20

if (!dir.exists(result.folder)){
    dir.create(result.folder, recursive = TRUE)
}

load_dir <- '/gpfs/gibbs/pi/zhao/yw599/Multiome/data_pp'

out <- peakRAM(
    
    function(){
        obj.rna <- load.rna.folder(rna.path, subset.rna, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
        obj.rna <- RenameCells(obj.rna, new.names = paste0(Cells(obj.rna), "-rna"))
        obj.atac <- load.atac.folder(atac.path, subset.atac, assay='RNA', min.cells=3, min.features=200, binarize=binarize)

        chunk <- function(x,n){
            vect <- c(1:x)
            num <- ceiling(x/n)
            split(vect,rep(1:num,each=n,len=x))
        }

        liger_list <- list()
        liger_names <- list()
        liger_list <- append(liger_list, obj.rna[["RNA"]]@counts)
        liger_names <- append(liger_names, 'RNA')
        liger_list <- append(liger_list, obj.atac[["RNA"]]@counts)
        liger_names <- append(liger_names, 'ATAC')
        names(liger_list) <- liger_names
        rm(obj.rna)
        rm(obj.atac)

        data_liger <- createLiger(liger_list)
        data_liger <- rliger::normalize(data_liger)
        if (datasets.use == 'RNA'){
            data_liger <- selectGenes(data_liger, datasets.use = 1)
        }else{
            data_liger <- selectGenes(data_liger)
        }
        if (sum(unlist(lapply(liger_list, ncol))) > 100000){
            data_liger@scale.data <- lapply(1:length(data_liger@norm.data),
                                         function(i) {rliger:::scaleNotCenterFast(t(data_liger@norm.data[[i]][data_liger@var.genes,]))})
            data_liger@scale.data <- lapply(data_liger@scale.data,function(l){
                l2 <- lapply(chunk(nrow(l),2000), function(i){as.matrix(l[i,])})
                res <- do.call(rbind,l2)
                return(res)
            })
            names(data_liger@scale.data) <- names(data_liger@norm.data)
            for (i in 1:length(data_liger@scale.data)) {
              data_liger@scale.data[[i]][is.na(data_liger@scale.data[[i]])] <- 0
              rownames(data_liger@scale.data[[i]]) <- colnames(data_liger@raw.data[[i]])
              colnames(data_liger@scale.data[[i]]) <- data_liger@var.genes
            }
            data_liger <-rliger:::removeMissingObs(data_liger, slot.use = "scale.data", use.cols = F)
            gc()
        } else {data_liger <- scaleNotCenter(data_liger)}
        data_liger <- optimizeALS(data_liger, k=npcs, rand.seed=1)
        data_liger <- quantile_norm(data_liger, rand.seed=1)
        latent <- as.data.frame(data_liger@H.norm)
        if (save) {
            if (datasets.use == 'RNA'){
                write.csv(latent, file.path(result.folder, "latent_useRNA.csv"))
            }else{
                write.csv(latent, file.path(result.folder, "latent.csv"))
            }
        }
    }
)

mem <- out$Peak_RAM_Used_MiB
time <- out$Elapsed_Time_sec
print(paste0('Running time: ', time, ' sec'))
print(paste0('Peak memory usage: ', mem, ' MB'))
if (datasets.use == 'RNA'){
    write(time, file.path(result.folder, "time_useRNA.txt"))
    write(mem, file.path(result.folder, "memory_useRNA.txt"))
}else{
    write(time, file.path(result.folder, "time.txt"))
    write(mem, file.path(result.folder, "memory.txt"))
}
