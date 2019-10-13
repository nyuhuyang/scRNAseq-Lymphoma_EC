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

#========== new analysis ======================
(load(file="data/Lymphoma_EC_10_20190922.Rda"))
Idents(object) =  "Doublets"
object %<>% subset(idents = "Singlet")
# cell.type_samples
object@meta.data$cell.type_samples = paste0(object@meta.data$cell.type,"_",
                                            object@meta.data$orig.ident)
(cell.type_samples = unique(object@meta.data$cell.type_samples))
Idents(object) = "date"
object %<>% subset(idents = c("2018-10-18", "2018-12-30"))
Idents(object) = "cell.type"
EC <- subset(object,idents = "Endothelial cells")
TC <- subset(object,idents = "T cells")

Idents(object) = "cell.type_samples"
DefaultAssay(EC) = "SCT"
EC <- FindVariableFeatures(object = EC, selection.method = "vst",
                           num.bin = 20,
                           mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

df_TEC = EC_exp$integrated
TEC_cor = cor(df_TEC,method = "spearman")
TEC_cor = reshape2::melt(TEC_cor)
colnames(TEC_cor)[3] = "correlation"

jpeg(paste0(path,"correlation",paste(cluster,collapse = "_"),".jpeg"), units="in", width=10, height=7,res=600)
ggplot(TEC_cor, aes(x = Var1, y = Var2, fill = correlation)) + geom_raster()+
        geom_text(aes(label = round(correlation, 3)))+
        scale_fill_gradientn(colors = ggsci::pal_gsea()(12)[5:12])+
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x = element_text(angle = 45, vjust=1,
                                         size =12, hjust=1))
dev.off()

# hcluster
Idents(object) = "cell.type_samples"
Idents(EC) = "cell.type_samples"
Idents(TC) = "cell.type_samples"

object_exp <- AverageExpression(object)
EC_exp <- AverageExpression(EC)
TC_exp <- AverageExpression(TC)

colnames(object_exp$SCT) = gsub("Endothelial cells_","EC:",colnames(object_exp$SCT))
colnames(object_exp$integrated) = gsub("Endothelial cells_","EC:",colnames(object_exp$integrated) )
colnames(object_exp$SCT) = gsub("T cells_","T:",colnames(object_exp$SCT))
colnames(object_exp$integrated) = gsub("T cells_","T:",colnames(object_exp$integrated) )

colnames(EC_exp$SCT) = gsub("Endothelial cells_","EC:_",colnames(EC_exp$SCT))
colnames(EC_exp$integrated) = gsub("Endothelial cells_","EC:",colnames(EC_exp$integrated) )

colnames(TC_exp$SCT) = gsub("T cells_","T:",colnames(TC_exp$SCT))
colnames(TC_exp$integrated) = gsub("T cells_","T:",colnames(TC_exp$integrated) )

exp_list <- list(EC_exp$SCT,EC_exp$integrated, TC_exp$SCT, TC_exp$integrated)
p <- list()
for(i in seq_along(exp_list)){
        dm<-as.dist(1-cor(exp_list[[i]]))
        hm<-hclust(dm,method = "average")
        p[[i]] <- ggdendro::ggdendrogram(hm)+ 
                labs(title= paste(ifelse((i %% 2) != 0, "Before", "After"),"Batch correction"))+
                TitleCenter()
        
}

jpeg(paste0(path,"Cluster_Dendrogram_E_SCT_interg.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(p[[1]],p[[2]])+ labs(title="Cluster Dendrogram for Endothelial cells")+
        TitleCenter(size =20)
dev.off()

jpeg(paste0(path,"Cluster_Dendrogram_T_SCT_interg.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(p[[3]],p[[4]])+ labs(title="Cluster Dendrogram for T cells")+
        TitleCenter(size =20)
dev.off()

# ========== old analysis ===========
(load(file="data/Lymphoma_EC_12_20190125.Rda"))
object = UpdateSeuratObject(object)
Idents(object) = "singler2sub"
TSNEPlot.1(object, group.by = "singler2sub", cols = ExtractMetaColor(object),
           no.legend = T,label = T, label.repel = T)
Idents(object) = "tests"
object %<>% subset(idents = c("test4", "test6"))

Idents(object) = "res.0.6"
TSNEPlot.1(object, group.by = "res.0.6", label = T)
object %<>% subset(idents = 10, invert = T)
object %<>% RenameIdents("0" = "T",
                         "1" = "EC",
                         "2" = "T",
                         "3" = "EC",
                         "4" = "T",
                         "5" = "T",
                         "6" = "T",
                         "7" = "T",
                         "8" = "T",
                         "9" = "EC",
                         "11" = "T",
                         "12" = "EC")
object@meta.data$cell.type = as.character(Idents(object))
# cell.type_samples
object@meta.data$cell.type_samples = paste0(object@meta.data$cell.type,"_",
                                            object@meta.data$orig.ident)
(cell.type_samples = unique(object@meta.data$cell.type_samples))

Idents(object) = "cell.type"
EC <- subset(object,idents = "EC")
TC <- subset(object,idents = "T")

Idents(EC) = "cell.type_samples"
Idents(TC) = "cell.type_samples"

EC <- subset(EC,idents = c("EC_TALL-only","EC_3119"), invert = T)
TC <- subset(TC,idents = c("T_EC-only","T_EC"), invert = T)

EC_exp <- AverageExpression(EC)
TC_exp <- AverageExpression(TC)

colnames(EC_exp$RNA) = gsub("EC_","EC:",colnames(EC_exp$RNA))
colnames(TC_exp$RNA) = gsub("T_","T:",colnames(TC_exp$RNA))

exp_list1 <- list(EC_exp$RNA, TC_exp$RNA)
p1 <- list()
for(i in seq_along(exp_list1)){
        dm<-as.dist(1-cor(exp_list1[[i]]))
        hm<-hclust(dm,method = "average")
        p1[[i]] <- ggdendro::ggdendrogram(hm)+TitleCenter()
        
}

jpeg(paste0(path,"Cluster_Dendrogram_old.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(p1[[1]],p1[[2]])+ labs(title="Cluster Dendrogram for all cells")+
        TitleCenter(size =20)
dev.off()