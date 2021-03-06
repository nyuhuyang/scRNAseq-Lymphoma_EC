########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
library(Matrix)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa.

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)#args=4,6,13
print(paste0("slurm_arrayid=",args))

# 3.1.1 load data
(load(file="data/Lymphoma_EC_10_20190922.Rda"))
Idents(object) =  "Doublets"
object %<>% subset(idents = "Singlet")

Idents(object) <-  "integrated_snn_res.0.6"
object %<>% RenameIdents("0" = "T cells",
                         "1" = "Endothelial cells",
                         "2" = "T cells",
                         "3" = "T cells",
                         "4" = "Endothelial cells",
                         "5" = "T cells",
                         "6" = "T cells",
                         "7" = "T cells",
                         "8" = "Endothelial cells",
                         "9" = "T cells",
                         "10" = "T cells")
object@meta.data$cell.type = as.character(Idents(object))
# cell.type_samples
object@meta.data$cell.type_samples = paste0(object@meta.data$cell.type,"_",
                                         object@meta.data$orig.ident)
(cell.type_samples = unique(object@meta.data$cell.type_samples))
ident.1.name = c("RO2-plus-EC","RO2-plus-EC-2","EC2")

Idents(object) = "cell.type_samples"
EC_markers <- FindMarkers.UMI(object, assay = "SCT",slot = "data",
                              ident.1 = cell.type_samples[args], #args=4,6,13
                              ident.2 = "Endothelial cells_EC",
                              test.use = "MAST",
                              min.pct = -Inf,
                              min.cells.feature = -Inf, 
                              min.cells.group = -Inf,
                              only.pos = F, # don't change it!
                              logfc.threshold = -Inf)
write.csv(EC_markers,paste0(path,"EC_markers_",cell.type_samples[args],"_vs_EC_SCT_data_.csv"))
