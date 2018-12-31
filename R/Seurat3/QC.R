########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/181227_Single_cell_TALL_sample_list.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("test",paste0("test1")))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample.id[sample_n]
conditions <- df_samples$conditions[sample_n]
projects <- df_samples$project[sample_n]
tissues <- df_samples$tissue[sample_n]

# check missing data
(current <- list.files("data")[!grepl(".Rda|RData",list.files("data"))])
(missing_data <- sample.id[!(sample.id %in% current)])

if(length(missing_data)>0){
        # Move files from Download to ./data and rename them
        species <- "hg19"
        for(missing_dat in missing_data){
                old.pth  <- paste("~/Downloads", missing_dat,"outs",
                                  "filtered_gene_bc_matrices",species,
                                  sep = "/")
                list.of.files <- list.files(old.pth)
                new.folder <- paste("./data", missing_dat,"outs",
                                    "filtered_gene_bc_matrices",
                                    species,sep = "/")
                if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
                # copy the files to the new folder
                file.copy(paste(old.pth, list.of.files, sep = "/"), new.folder)
                list.files(new.folder)
        }
}

## Load the dataset
object_raw <- list()
object_Seurat <- list()
for(i in 1:length(samples)){
        object_raw[[i]] <- Read10X(data.dir = paste0("data/",sample.id[i],
                                   "/outs/filtered_gene_bc_matrices/hg19"))
        colnames(object_raw[[i]]) = paste0(samples[i],"_",colnames(object_raw[[i]]))
        rownames(object_raw[[i]]) = gsub("hg19_","",rownames(object_raw[[i]]))
        object_Seurat[[i]] <- CreateSeuratObject(object_raw[[i]],
                                                 project = projects[i],
                                                 min.cells = 0,
                                                 min.features = 0)
        object_Seurat[[i]]$conditions <- conditions[i]
        object_Seurat[[i]]$tissues <- tissues[i]
}

#======1.1.2 QC before merge =========================
cell.number <- sapply(object_Seurat, function(x) length(colnames(x)))
QC_list <- lapply(object_Seurat, function(x) as.matrix(GetAssayData(x, slot = "counts")))
median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))

min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))

QC.list <- cbind(df_samples[sample_n,],cell.number, median.nUMI, median.nGene, 
                 min.nUMI,min.nGene, row.names = samples)
write.csv(QC.list,paste0(path,"QC_list.csv"))
QC.list %>% kable() %>% kable_styling()
remove(QC_list,median.nUMI,median.nGene,min.nUMI,min.nGene,QC.list);GC()

#========1.1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_Seurat)
remove(object_raw,object_Seurat);GC()

# read and select mitochondial genes
#all.mito.genes <- read.csv("./doc/Mitochondrial.csv",row.names = 1) %>% rownames %>%
#        tolower %>% Hmisc::capitalize()
#(mito.genes <-  rownames(x = object@data) %in% all.mito.genes %>%
#        rownames(x = object@data)[.])
(mito.features <- grep(pattern = "^MT-", x = rownames(object), value = TRUE))
percent.mito <- Matrix::colSums(GetAssayData(object, slot = 'counts')[mito.features, ])/
        Matrix::colSums(GetAssayData(object, slot = 'counts'))

object[["percent.mito"]] <- percent.mito
Idents(object) = factor(Idents(object),levels = samples)
g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mito"), function(features){
        VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)
        })
save(g1,file= paste0(path,"g1_11_20181230.Rda"))
