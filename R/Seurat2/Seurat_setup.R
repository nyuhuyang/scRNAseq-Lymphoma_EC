########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

library(Seurat)
library(magrittr)
library(harmony)
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
df_samples <- readxl::read_excel("doc/181230_Single_cell_TALL_sample_list.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$conditions %in% c("EC","EC+T-ALL","EC+T-ALL+anti-Notch1","T-ALL"))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample.id[sample_n]
conditions <- df_samples$conditions[sample_n]
projects <- df_samples$project[sample_n]
tissues <- df_samples$tissue[sample_n]
tests <- df_samples$tests[sample_n]

#======1.2 load  SingleCellExperiment =========================
(load(file = "./data/sce_12_20190125.Rda"))
names(sce_list)
object_Seurat <- lapply(sce_list, as.seurat) %>%
    #lapply(NormalizeData) %>%
    #lapply(ScaleData) %>%
    lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
    object_Seurat[[i]]@meta.data$conditions <- conditions[i]
    object_Seurat[[i]]@meta.data$projects <- projects[i]
    object_Seurat[[i]]@meta.data$tissues <- tissues[i]
    object_Seurat[[i]]@meta.data$tests <- tissues[i]
}
# we will take the union of the top 1k variable genes in each dataset for alignment
object_Seurat$`3139PR` = NULL

genes.use <- object_Seurat %>% 
    lapply(function(object) head(rownames(object@hvg.info), 1200)) %>%
    unlist %>% unique
length(genes.use)

#========1.3 merge ===================================

object <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), object_Seurat)
object@var.genes = genes.use
remove(sce_list,object_Seurat);GC()

(remove <- which(colnames(object@meta.data) %in%c("Sample","is_cell_control",
                                           "pct_counts_in_top_500_features_Mito")))
object@meta.data = object@meta.data[,-c(3,5:40)]

object = SetAllIdent(object, id = "orig.ident")
#======1.4 mito, QC, filteration =========================

(mito.genes <-  grep(pattern = "^MT-", x = rownames(x = object@data), value = TRUE))
percent.mito <- Matrix::colSums(object@raw.data[mito.genes, ])/Matrix::colSums(object@raw.data)
object <- AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito")

(load(file = "./output/20190125/g1_11_20190125.Rda"))

object <- FilterCells(object = object, subset.names = c("nGene","nUMI","percent.mito"),
                      low.thresholds = c(500,1000, -Inf), 
                      high.thresholds = c(Inf,Inf, 0.1))

object@ident = factor(object@ident,levels = samples)
g2 <- lapply(c("nGene", "nUMI", "percent.mito"), function(features){
    VlnPlot(object = object, features.plot = features, nCol = 3, 
            point.size.use = 0.2,size.x.use = 10, group.by = "ident",
            x.lab.rot = T, do.return = T)})
save(g2,file = paste0(path,"g2_11_20190125.Rda"))
jpeg(paste0(path,"/S1_nGene~.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                    scale_y_log10(limits = c(100,10000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                    scale_y_log10(limits = c(100,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI~.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                    scale_y_log10(limits = c(800,100000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                    scale_y_log10(limits = c(800,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito~.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                    ylim(c(0,0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,0.5))))
dev.off()
#======1.5 FindVariableGenes=======================
object <- NormalizeData(object = object)
jpeg(paste0(path,"/S1_dispersion.jpeg"), units="in", width=10, height=7,res=600)
object <- FindVariableGenes(object = object, mean.function = ExpMean, 
                            dispersion.function = LogVMR, do.plot = T, top.genes = 1000
                            x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.65)
dev.off()
length(object@var.genes)
#table(object@var.genes %in% genes.use)

#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "./data/seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- HumanGenes(object,cc.genes[1:43])
g2m.genes <- HumanGenes(object,cc.genes[44:97])
object <- CellCycleScoring(object = object, s.genes = s.genes, g2m.genes = g2m.genes, 
                           set.ident = TRUE)
RidgePlot(object = object, features.plot = HumanGenes(object,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
object@meta.data$CC.Difference <- object@meta.data$S.Score - object@meta.data$G2M.Score
object@meta.data$S.Score = object@meta.data$S.Score - min(object@meta.data$S.Score)
object@meta.data$G2M.Score = object@meta.data$G2M.Score - min(object@meta.data$G2M.Score)
tail(x = object@meta.data)


#======1.6 PCA =========================
object %<>% NormalizeData %>% ScaleData %>%
    RunPCA(pc.genes = object@var.genes, pcs.compute = 30, do.print = F)

jpeg(paste0(path,"/S1_PCElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
PCElbowPlot(object, num.pc = 50)
dev.off()

jpeg(paste0(path,"/S1_PCHeatmap.jpeg"), units="in", width=10, height=7,res=600)
PCHeatmap(object, pc.use = c(1:3,25:30), cells.use = 500, do.balanced = TRUE)
dev.off()

GC()
system.time({
    object %<>% RunTSNE(reduction.use = "pca", dims.use = 1:30, do.fast = TRUE)
    object %<>% FindClusters(reduction.type = "pca", resolution = 0.6, dims.use = 1:30,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})

p0 <- TSNEPlot(object, do.return = T, pt.size = 0.5, group.by = "orig.ident")

#======1.6 RunHarmony=======================
jpeg(paste0(path,"/S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony("orig.ident", dims.use = 1:30,
                                   theta = 2, plot_convergence = TRUE,
                                   nclust = 50, max.iter.cluster = 100))
dev.off()

object@ident %<>% factor(levels = samples)
p1 <- DimPlot(object = object, reduction.use = "harmony", pt.size = .1, group.by = "orig.ident", do.return = T)
p2 <- VlnPlot(object = object, features.plot = "Harmony1", group.by = "orig.ident", do.return = TRUE,
              x.lab.rot = T)
jpeg(paste0(path,"/S1_Harmony_vplot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1,p2)
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time({
    object %<>% RunTSNE(reduction.use = "harmony", dims.use = 1:30, do.fast = TRUE)
    object %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = 1:30,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})

system.time({
        object %<>% RunUMAP(reduction.use = "harmony", dims.use = 1:30)
        object %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = 1:30,
                                 save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                                 force.recalc = TRUE, print.output = FALSE)
})


p3 <- TSNEPlot(object, do.return = T, pt.size = 0.5, group.by = "orig.ident")
jpeg(paste0(path,"/S1_pca_vs_harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Raw data")+
              theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p3+ggtitle("After alignment")+
              theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) )
dev.off()

g_Harmony <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                        do.return = TRUE, no.legend = T, 
                        pt.size = 1,label.size = 6 )

p4 <- TSNEPlot.1(object, do.return = T, pt.size = 0.5, group.by = "orig.ident",
               do.label = T, label.repel =T,no.legend = T)
jpeg(paste0(path,"/TSNEplot_cluster_samples~.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p4+ggtitle("Tsne plot of all clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          g_Harmony+ggtitle("Tsne plot of all samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) )
dev.off()


g_UMAP <- DimPlot.1(object = object, reduction.use = "umap", 
                    do.label = T, group.by = "ident",label.repel =T,
                    do.return = TRUE, no.legend = T, 
                    pt.size = 1,label.size = 6)

p5 <- DimPlot.1(object, reduction.use = "umap", do.return = T, pt.size = 0.5, group.by = "orig.ident",
                do.label = T, label.repel =T,no.legend = T)

jpeg(paste0(path,"/UMAPplot_cluster_samples~.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p5+ggtitle("UMAP plot of all clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          g_UMAP+ggtitle("UMAP plot of all samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) )
dev.off()

jpeg(paste0(path,"/S1_tsne_umap~.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g_Harmony+ggtitle("Tsne plot of all clusters")+
                  theme(text = element_text(size=15),							
                        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          g_UMAP+ ggtitle("UMAP plot of all clusters")+
                  theme(text = element_text(size=15),							
                        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 
)
dev.off()

saveRDS(object@scale.data, file = "data/Lymphoma_EC.scale.data_Harmony_12_20190125.rds")
object@scale.data = NULL
save(object, file = "./data/Lymphoma_EC_12_20190125.Rda")
