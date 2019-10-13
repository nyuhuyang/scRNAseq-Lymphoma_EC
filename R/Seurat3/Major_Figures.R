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

# load data
# 3.1.1 load data
(load(file="data/Lymphoma_EC_10_20190922.Rda"))

# Fig 4F
Idents(object) = "cell.type"
TSNEPlot.1(object, cols = ExtractMetaColor(object), 
           title.size =15,
           no.legend = T, width = 5.4, height=4,do.print = T)

# Fig 4G
Idents(object) =  "Doublets"
object %<>% subset(idents = "Singlet")
object@meta.data$cell.type_samples = paste0(object@meta.data$cell.type,"_",
                                         object@meta.data$orig.ident)
(cell.type_samples = unique(object@meta.data$cell.type_samples))
Idents(object) = "cell.type_samples"
RO2_EC <- subset(object, idents = c("Endothelial cells_EC",
                                    "Endothelial cells_RO2-plus-EC",
                                    "Endothelial cells_RO2-plus-EC-2"))


Idents(RO2_EC) = "conditions"
DefaultAssay(RO2_EC) = "SCT"
RO2_EC_markers <- FindAllMarkers.UMI(RO2_EC,assay = "SCT",slot = "data",
                                     test.use = "MAST",
                                     logfc.threshold = 0.1,
                                     only.pos = T)
Top_n =50
top = RO2_EC_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
RO2_EC %<>% ScaleData(features=top$gene)
table(rownames(RO2_EC@assays$SCT@scale.data) %in% top$gene)
DoHeatmap.1(RO2_EC, features = top$gene, Top_n = Top_n,
            do.print=T, angle = 0, group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=F, cex.row=4, width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"differentially expressed genes in Endothelial cells"))
remove(RO2_EC);GC()
# Fig 5F
RO2_T <- subset(object, idents = c("T cells_RO2",
                                    "T cells_RO2-plus-EC",
                                    "T cells_RO2-plus-EC-2"))


Idents(RO2_T) = "conditions"
DefaultAssay(RO2_T) = "SCT"
RO2_T_markers <- FindAllMarkers.UMI(RO2_T,assay = "SCT",slot = "data",
                                     test.use = "MAST",
                                     logfc.threshold = 0.1,
                                     only.pos = T)
Top_n =50
top = RO2_T_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
RO2_T %<>% ScaleData(features=top$gene)
table(rownames(RO2_T@assays$SCT@scale.data) %in% top$gene)
DoHeatmap.1(RO2_T, features = top$gene, Top_n = Top_n,
            do.print=T, angle = 0, group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=F, cex.row=4, width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"differentially expressed genes in T cells"))
remove(RO2_T);GC()

# Fig 5G
T_RO2_3119 <- subset(object, idents = c("T cells_RO2",
                                        "T cells_RO2-plus-EC",
                                        "T cells_RO2-plus-EC-2",
                                        "T cells_3119-plus-EC",
                                        "T cells_3119-plus-EC2",
                                        "T cells_3119"))
T_RO2_3119 %<>% RenameIdents("T cells_3119-plus-EC2" = "3119 T-ALL + EC rep 2",
                             "T cells_3119-plus-EC" = "3119 T-ALL + EC rep 1",
                             "T cells_3119" = "3119 T-ALL",
                             "T cells_RO2-plus-EC-2" = "RO2 T-ALL + EC rep 2",
                             "T cells_RO2-plus-EC" = "RO2 T-ALL + EC rep 1",
                             "T cells_RO2" = "RO2 T-ALL")
DefaultAssay(T_RO2_3119) = "SCT"
table(Idents(T_RO2_3119))
jpeg(paste0(path,"Fig_5G_Dotplot.jpeg"),units="in", width=8, height=5,res=600)
DotPlot(T_RO2_3119, features = c("DTX1","NOTCH3","HES1","HES4","IGF1R")) + 
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
