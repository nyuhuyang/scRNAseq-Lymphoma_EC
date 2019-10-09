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
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# 3.1.1 load data
(load(file="data/Lymphoma_EC_10_20190922.Rda"))
(load(file = "data/sce_10_20191002.Rda"))
reconstructed = as.matrix(sce@assays$data$reconstructed)
rownames(reconstructed) = rownames(sce)
RNA_data = as(reconstructed, "sparseMatrix")*log(2)
RNA_data = RNA_data[rownames(object),colnames(object)]
object@assays$RNA@data = RNA_data
Idents(object) <-  "Doublets"
object %<>% subset(idents = "Singlet")

# subset date
Idents(object) = "date"
(dates <- unique(object$date))
(date <- dates[args])
object %<>% subset(idents = date)

# Rename ident
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
# subset cell types
object@meta.data$cell.type_conditions = paste0(object@meta.data$cell.type,"_",
                                               object@meta.data$conditions)
Idents(object) = "cell.type_conditions"
object %<>% subset(idents = c("Endothelial cells_EC","Endothelial cells_EC+T-ALL"))
object %<>% RenameIdents("Endothelial cells_EC" = "EC",
                         "Endothelial cells_EC+T-ALL" = "TEC")
object@meta.data$cell.type_conditions = as.character(Idents(object))
Idents(object) = "cell.type_conditions"
object %<>% sortIdent()
table(Idents(object))
assay = c("RNA","SCT","SCT"); slot = c("data","data","scale.data")
for(i in 1){
        EC_markers <- FindMarkers.UMI(object, assay = assay[i],slot = slot[i],
                                      ident.1 = "TEC", ident.2 = "EC",
                                      test.use = "MAST",
                                      min.pct = -Inf,
                                      min.cells.feature = -Inf, 
                                      min.cells.group = -Inf,
                                      only.pos = F, # don't change it!
                                      logfc.threshold = -Inf)
        write.csv(EC_markers,paste0(path,"EC_markers_",date,"_",assay[i],"_",slot[i],".csv"))
        svMisc::progress(i/length(assay)*100)
}
