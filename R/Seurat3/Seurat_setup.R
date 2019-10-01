########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(kableExtra)
library(magrittr)
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

#========1.3 merge ===================================
# https://support.bioconductor.org/p/76745/
# Don't use ComBat on raw counts
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
object@assays$RNA@data = object@assays$RNA@data *log(2)
remove(sce_list,object_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = meta.data 
remove(meta.data);GC()

object_data = object@assays$RNA@data
remove(object);GC()
(load(file="data/Lymphoma_EC_10_20190922.Rda"))

DefaultAssay(object) = "RNA"
object_data = object_data[rownames(object),colnames(object)]

table(rownames(object_data) == rownames(object))
table(colnames(object_data) == colnames(object))
object_data[1:5,1:5]
object@assays$RNA@data = object_data
#======1.3 batch-correct using ComBat (skip)=========================
table(batch.effect <- object@meta.data[,"tests"])
names(batch.effect) = rownames(object@meta.data)

data = object@assays$RNA@data
table(rowSums(data)>0)
data = data[rowSums(data)>0,]

jpeg(paste0(path,"Combat_data.jpeg"), units="in", width=10, height=7,res=600)
data_combat = ComBat(as.matrix(data), batch.effect, prior.plots=T, par.prior=TRUE)
dev.off()

min(matrixStats::rowMins(data_combat))
data_combat[data_combat < 0] = 0
object@assays$RNA@data = as(data_combat, "sparseMatrix")
remove(data,data_combat,batch.effect);GC()


######################################

# After removing unwanted cells from the dataset, the next step is to normalize the data.
#object <- NormalizeData(object = object, normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(object), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(path,"VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()
#======1.3 1st run of pca-tsne  =========================
object <- ScaleData(object = object,features = VariableFeatures(object))
object <- RunPCA(object, features = VariableFeatures(object),verbose =F,npcs = 100)
object <- JackStraw(object, num.replicate = 20,dims = 100)
object <- ScoreJackStraw(object, dims = 1:100)

jpeg(paste0(path,"JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 60:70)+
    ggtitle("JackStrawPlot")+
    theme(text = element_text(size=15),	
          plot.title = element_text(hjust = 0.5,size = 18))
dev.off()
npcs =65
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                       dims.use = 1:npcs, print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs, check_duplicates = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
object@assays$RNA@scale.data = matrix(0,0,0)

object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
p0 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 do.return = T,no.legend = F,label.size = 4, repel = T, title = "Original")
p1 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 no.legend = F,label.size = 4, repel = T, title = "Original")

#======1.5 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "tests")
remove(object);GC()
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
options(future.globals.maxSize= object.size(object_list)*1.5)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                  anchor.features = object.features)
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

remove(anchors,object_list);GC()
object %<>% RunPCA(npcs = 100, verbose = FALSE)
npcs =65
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                         dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
p2 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 label.size = 4, repel = T,title = "Intergrated tSNE plot")
p3 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 label.size = 4, repel = T,title = "Intergrated UMAP plot")
#=======1.9 summary =======================================
jpeg(paste0(path,"S1_remove_batch_tsne.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Clustering without integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p2+ggtitle("Clustering with integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

jpeg(paste0(path,"S1_remove_batch_umap.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1+ggtitle("Clustering without integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p3+ggtitle("Clustering with integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

TSNEPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.6", 
         do.return = F, no.legend = F, title = "tSNE plot for all clusters",
         pt.size = 0.3,alpha = 1, label.size = 5, do.print = T)

UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.6", 
           do.return = F, no.legend = F, title = "UMAP plot for all clusters",
           pt.size = 0.2,alpha = 1, label.size = 5, do.print = T)

object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = "data/Lymphoma_EC_10_20190922.Rda")
object_data = object@assays$RNA@data
save(object_data, file = "data/Lymphoma.data_EC_10_20190922.Rda")
