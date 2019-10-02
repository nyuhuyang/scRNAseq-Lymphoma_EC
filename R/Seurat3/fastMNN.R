library(Seurat)
library(scater)
library(magrittr)
library(batchelor)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/20190814_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$group %in% c("3119","EC","RO2"))
df_samples <- df_samples[sample_n,]
print(df_samples)
(samples = df_samples$sample)

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_10_20190922.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.Seurat)

for(i in 1:length(samples)){
        object_list[[i]]$orig.ident <- df_samples$sample[i]
        object_list[[i]]$conditions <- df_samples$conditions[i]
        object_list[[i]]$group <- df_samples$group[i]
        object_list[[i]]$project <- df_samples$project[i]
        object_list[[i]]$tests <- df_samples$tests[i]
        object_list[[i]]$tissue <- df_samples$tissue[i]
        object_list[[i]]$date <- gsub(" UTC","",df_samples$date[i])
        Idents(object_list[[i]]) <- df_samples$sample[i]
}

#======1.4 mnncorrect =========================
# https://support.bioconductor.org/p/106010/
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
remove(sce_list,object_list);GC()
sce <- as.SingleCellExperiment(object)
system.time(sce <- fastMNN(sce, batch = object$orig.ident))
save(sce, file = paste0("data/","sce_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))
