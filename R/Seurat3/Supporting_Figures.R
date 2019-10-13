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

#S4C
jpeg(paste0(path,"Fig_S4C1_T.jpeg"),units="in", width=6, height=2.4,res=600)
FeaturePlot(object, features = c("CD3G","CD3D"),reduction = "tsne")
dev.off()

jpeg(paste0(path,"Fig_S4C1_EC.jpeg"),units="in", width=6, height=2.4,res=600)
FeaturePlot(object, features = c("VWF","CDH5"),reduction = "tsne")
dev.off()

#S4D
(load(file="output/singler_Lymphoma_EC_10T_20190922.Rda"))
sce_object <- as.SingleCellExperiment(object)
jpeg(paste0(path,"DrawHeatmap_cell.type.jpeg"), units="in", width=10, height=7,
     res=600)
print(plotScoreHeatmap(singler, show.labels = TRUE,
                 annotation_col=data.frame(cell.type=sce_object$cell.type,
                                           row.names=rownames(singler))))
dev.off()

#S4F
Idents(object) = "conditions"
sub_object = subset(sub_object, idents = c("EC", "EC+T-ALL"))
FeaturePlot.1(sub_object, features = "IGFBP7",split.by = "conditions",
              reduction = "tsne",do.print = T,width=6, height=3,
              cols =c("lightgrey", "blue"))

#S4G (only RO2)
Idents(object) = "date"
sub_object = subset(object, idents = "2018-10-18")
Idents(sub_object) = "conditions"
sub_object = subset(sub_object, idents = c("EC", "EC+T-ALL"))

jpeg(paste0(path,"Fig_S4G_Dotplot.jpeg"),units="in", width=6, height=4,res=600)
DotPlot(sub_object, features = c("ETS2","ETS1","DLL1","JAG1","IGFBP7","CD59")) + 
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# S5J
Idents(object) =  "Doublets"
object %<>% subset(idents = "Singlet")
object@meta.data$cell.type_samples = paste0(object@meta.data$cell.type,"_",
                                            object@meta.data$orig.ident)
(cell.type_samples = unique(object@meta.data$cell.type_samples))
Idents(object) = "cell.type_samples"
TC <- subset(object, idents = c("T cells_3119-plus-EC2",
                                "T cells_3119" ))

TC %<>% RenameIdents("T cells_3119"= "3119 T_ALL",
                     "T cells_3119-plus-EC2"="3119 T_ALL + EC")

DefaultAssay(TC) = "SCT"
TC_markers <- FindAllMarkers.UMI(TC,assay = "SCT",slot = "data",
                                    test.use = "MAST",
                                    logfc.threshold = 0.1,
                                    only.pos = T)
Top_n =50
top = TC_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
TC %<>% ScaleData(features=top$gene)
table(rownames(TC@assays$SCT@scale.data) %in% top$gene)
DoHeatmap.1(TC, features = top$gene, Top_n = Top_n,
            do.print=T, angle = 0, group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=F, cex.row=4, width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"differentially expressed genes in 3119 T-ALL cells"))
remove(TC);GC()