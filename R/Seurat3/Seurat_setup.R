########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(scran)
library(kableExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(harmony)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
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

#======1.2 load  SingleCellExperiment =========================
(load(file = "./data/sce_11_20181230.Rd"))
names(sce_list)
object_Seurat <- lapply(sce_list, as.Seurat) %>%
    lapply(FindVariableFeatures, selection.method = "vst",
           nfeatures = 2000, verbose = FALSE)

for(i in 1:length(samples)){
    object_Seurat[[i]]$conditions <- conditions[i]
    object_Seurat[[i]]$tissues <- tissues[i]
}
#========1.3 merge ===================================
# Integration datasets
# dims = 10 and 50
#object.anchors <- FindIntegrationAnchors(object.list = object_Seurat, dims = 1:30)
#object <- IntegrateData(anchorset = object.anchors, dims = 1:30)
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_Seurat)

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
#DefaultAssay(object = object) <- "integrated"
remove(sce_list,object_Seurat,object.anchors);GC()

#======1.4 mito, QC, filteration =========================
(mito.features <- grep(pattern = "^MT-", x = rownames(object),value = TRUE))
percent.mito <- Matrix::colSums(GetAssayData(object, slot = 'data')[mito.features, ])/
    Matrix::colSums(GetAssayData(object, slot = 'data'))
object[["percent.mito"]] <- percent.mito
Idents(object) = factor(Idents(object),levels = samples)
(load(file = "output/20181230/g1_11_20181230.Rda"))

object %<>% subset(subset = nFeature_RNA > 700 & nCount_RNA > 900 & percent.mito < 0.05)
# FilterCellsgenerate Vlnplot before and after filteration
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mito"), function(features){
    VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)
})

save(g2,file= paste0(path,"g2_11_20181230.Rda"))
jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                    scale_y_log10(limits = c(100,10000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                    scale_y_log10(limits = c(100,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                    scale_y_log10(limits = c(500,100000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                    scale_y_log10(limits = c(500,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                    ylim(c(0,0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,0.5))))
dev.off()

FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mito")

#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../R/seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- HumanGenes(object=object ,cc.genes[1:43])
g2m.genes <- HumanGenes(object,cc.genes[44:97])
object <- CellCycleScoring(object = object, s.features = s.genes, g2m.features = g2m.genes)
"Error: Insufficient data values to produce 24 bins."
RidgePlot(object = object, features = HumanGenes(object,c("CCND1","CDK4","CCND2")))

#======1.6 NormalizeData and ScaleData =========================
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 1e4)

object <- FindVariableFeatures(object, selection.method = 'mean.var.plot',
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(0.75, Inf))
length(x = VariableFeatures(object))
object <- ScaleData(object, features= rownames(object))

#======1.6 PCA Determine statistically significant principal components=======================
# Run the standard workflow for visualization and clustering
object <- RunPCA(object, features = VariableFeatures(object),
                 npcs = 50, verbose = FALSE)

jpeg(paste0(path,"/S1_PCElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 50)
dev.off()
jpeg(paste0(path,"/S1_PCHeatmap.jpeg"), units="in", width=10, height=7,res=600)
DimHeatmap(object, dims = c(1:3,35:40), cells = 500, balanced = TRUE)
dev.off()

RunUMAP

object <- JackStraw(object, num.replicate = 100,dims = 50)
object <- ScoreJackStraw(object, dims = 1:50)
p4 <- JackStrawPlot(object = object, dims = 30:40)
jpeg(paste0(path,"/S1_JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
p4
dev.off()

object <- FindNeighbors(object, dims = 1:37)
object <- FindClusters(object, resolution = 0.6)

#======1.7 RunHarmony=======================
object <- RunTSNE(object, dims = 1:37)
p1 <- DimPlot(object, reduction = 'tsne',group.by="orig.ident",pt.size = 1,label = T,
              label.size = 4, repel = T)
system.time(object %<>% RunHarmony.1(group.by= "orig.ident", dims.use = 1:30,
                                   theta = 2, plot_convergence = TRUE,
                                   nclust = 50, max.iter.cluster = 100))
Idents(object) %<>% factor(levels = samples)


jpeg(paste0(path,"/S1_Harmony_DimHeatmap.jpeg"), units="in", width=10, height=7,res=600)
DimHeatmap(object = object, reduction.type = "harmony", cells.use = 500, 
           dim.use = c(1:3,18:20), do.balanced = TRUE)
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time({
    object %<>% RunTSNE(reduction.use = "harmony", dims.use = 1:50, do.fast = TRUE)
    object %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = 1:50,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})

p2 <- DimPlot(object, reduction = 'tsne',group.by="orig.ident",pt.size = 1,label = T,
              label.size = 4, repel = T)
jpeg(paste0(path,"/S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1, p2)
dev.off()

g_Harmony <- DimPlot(object = object, label = F, group.by = "ident",
                        do.return = TRUE, no.legend = F, 
                        #colors.use = ExtractMetaColor(object),
                        pt.size = 1,label.size = 6 )+
    ggtitle("Tsne plot of all cell types")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot-Harmony.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()

save(object, file = "./data/Lymphoma_EC_11_20181231.Rda")
