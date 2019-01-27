library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "./data/Lymphoma_EC_12_20190125.Rda"))
(load(file = "../SingleR/data/Blueprint_encode.RData"))
object@scale.data = NULL
GC();GC();GC();GC();GC();GC();GC();GC();GC();
singler = CreateSinglerObject(object@data, annot = NULL, project.name="T_EC",
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              ref.list = list(Blueprint_encode), normalize.gene.length = F, 
                              variable.genes = "de",fine.tune = F, do.signatures = F, clusters = NULL)
# singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object@data)
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = object@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file = paste0(path,"singler_Lymphoma_EC_12_20190125.Rda"))
