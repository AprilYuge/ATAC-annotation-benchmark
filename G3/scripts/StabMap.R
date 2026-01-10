library(peakRAM)
save <- TRUE
binarize <- FALSE

.libPaths(c(.libPaths(), "/gpfs/gibbs/pi/zhao/yw599/conda_envs/myR/lib/R/library"))
library(Dict)
library(Seurat)
library(Signac)
source("/gpfs/gibbs/pi/zhao/yw599/Multiome/utils.r")

# library(EnsDb.Hsapiens.v75)
library(StabMap)
# library(SingleCellMultiModal)
library(scran)
library(scuttle)
library(stringr)

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
        rna <- load.rna.folder(rna.path, subset.rna, assay='RNA', min.cells=3, min.features=200, binarize=binarize)
        rna <- RenameCells(rna, new.names = paste0(Cells(rna), "-rna"))
        atac <- load.atac.folder(atac.path, subset.atac, assay='RNA', min.cells=3, min.features=200, binarize=binarize)

        # if (tissue %in% c('brain', 'lung_facs', 'lung_droplet')){
        #     counts <- GetAssayData(rna, assay = "RNA", slot = "counts") 
        #     rownames(counts) <- rename.rna.genes(rownames(counts))
        #     rna <- CreateSeuratObject(counts = counts, meta.data = rna@meta.data, assay = 'RNA')
        # }

        # Convert data to SingleCellExperiment
        rna <- SingleCellExperiment(list(counts = GetAssayData(rna, assay = "RNA", slot = "counts")), colData=rna@meta.data)
        atac <- SingleCellExperiment(list(counts = GetAssayData(atac, assay = "RNA", slot = "counts")), colData=atac@meta.data)

        # Process RNA
        if (!grepl('Spleen', rna.path)){
            rna <- scuttle::logNormCounts(rna)
        }else{
            logcounts(rna) <- counts(rna)
        }
        decomp <- modelGeneVar(rna)
        rna.hvgs <- rownames(decomp)[decomp$mean > 0.01 & decomp$p.value <= 0.05]
        rna <- rna[rna.hvgs,]

        # Process ATAC
        peak_path = paste0(result.folder, '/region2gene.csv')
        peakInfo <- read.csv(peak_path, row.names=1)
        peakInfo <- peakInfo[peakInfo$gene %in% rna.hvgs,]
        peakInfo <- peakInfo[!duplicated(peakInfo$peak),]
        rownames(peakInfo) <- peakInfo$peak
        logcounts(atac) <- tfidf(1*(counts(atac)>0))
        decomp <- modelGeneVar(atac)
        if (grepl('kidney_paired', result.folder)){
            thr_m <- 0.1
        }else{
            thr_m <- 0.25
        }
        atac.hvgs <- rownames(decomp)[decomp$mean > thr_m & decomp$p.value <= 0.05]
        atac.hvgs_all <- union(atac.hvgs, intersect(peakInfo$peak, rownames(decomp)[decomp$mean > thr_m]))
        atac <- atac[atac.hvgs_all,]

        # length(intersect(rownames(atac), rownames(rna)))
        rnames.atac = rownames(atac)
        idx <- rownames(atac) %in% peakInfo$peak
        rnames.atac[idx] <- peakInfo[rnames.atac[idx], 'gene']
        rnames.rna = rownames(rna)
        print(paste0('Number of common features is ', length(intersect(rnames.atac, rnames.rna))))
        print(paste0('Number of unique features in ATAC is ', length(setdiff(rnames.atac, rnames.rna))))
        print(paste0('Number of unique features in RNA is ', length(setdiff(rnames.rna, rnames.atac))))
        # rownames(rna) <- rnames.rna
        rownames(atac) <- rnames.atac
        # nrow(atac)
        # length(unique(rownames(atac)))
        atac <- atac[!duplicated(rownames(atac)),]
        # nrow(rna)
        # length(unique(rownames(rna)))
        # rna <- rna[!duplicated(rownames(rna)),]
        # length(intersect(rownames(atac), rownames(rna)))

        SCE_list = list(atac = atac, rna = rna)
        assayNames = list(atac = "logcounts", rna = "logcounts")
        assay_list = mapply(assay, SCE_list, assayNames)
        latent = stabMap(assay_list,
                         reference_list = c('atac', 'rna'),
        #                  projectAll = TRUE,
                         scale.center = TRUE,
                         scale.scale = TRUE,
                         plot=FALSE,
                         maxFeatures=5000)
        latent <- reWeightEmbedding(latent)
        if (save){
            write.csv(latent[,1:50], file.path(result.folder, "latent_to_atac.csv"))
            write.csv(latent[,51:100], file.path(result.folder, "latent_to_rna.csv"))
        }
    }
)

mem <- out$Peak_RAM_Used_MiB
time <- out$Elapsed_Time_sec
print(paste0('Running time: ', time, ' sec'))
print(paste0('Peak memory usage: ', mem, ' MB'))
write(time, file.path(result.folder, "time.txt"))
write(mem, file.path(result.folder, "memory.txt"))
