library(Seurat)
library(SingleR)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggthemes)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file = "data/Lymphoma_EC_10_20190922.Rda"))
(load(file="output/singler_Lymphoma_EC_10T_20190922.Rda"))

singler@nrows == ncol(object)
# if singler didn't find all cell labels
if(singler@nrows < ncol(object)) object = subset(object, cells = singler@rownames)
table(singler@rownames == colnames(object))
##############################
# add singleR label to Seurat
###############################

singlerDF = data.frame("labels" = singler$labels,
                       row.names = singler@rownames)
head(singlerDF)

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
sce_object <- as.SingleCellExperiment(object)
jpeg(paste0(path,"DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
plotScoreHeatmap(singler, show.labels = TRUE,
                 annotation_col=data.frame(cell.type=sce_object$cell.type,
                                           row.names=rownames(singler)))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

##############################
# process color scheme
##############################
object <- AddMetaData(object = object,metadata = singlerDF)
object <- AddMetaColor(object = object, label= "labels", colors = Singler.colors)
Idents(object) <- "labels"

TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T, label.repel = T,pt.size = 1,
           label.size = 4, repel = T,no.legend = T,do.print = T,
           title = "Cell types identified by Blueprint + Encode")
save(object,file="data/Lymphoma_EC_10_20190922.Rda")