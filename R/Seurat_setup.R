########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(scran)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/181018_Single_cell_sample list.xlsx")
samples <- df_samples$samples
projects <- df_samples$projects
conditions <- df_samples$conditions
df_samples %>% kable() %>% kable_styling()

CancerCell_raw <- list()
CancerCell_Seurat <- list()
for(i in 1:length(samples)){
    CancerCell_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                             samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
    colnames(CancerCell_raw[[i]]) <- paste0(samples[i],
                                            "_",colnames(CancerCell_raw[[i]]))
    CancerCell_Seurat[[i]] <- CreateSeuratObject(CancerCell_raw[[i]],
                                                 min.cells = 3,
                                                 min.genes = 200,
                                                 project = projects[i],
                                                 names.delim = "_")
    CancerCell_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
CancerCell <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), CancerCell_Seurat)
remove(CancerCell_raw,CancerCell_Seurat);GC()
CancerCell <- FilterCells(CancerCell, subset.names = "nGene",
                    low.thresholds = 200,
                    high.thresholds = Inf) %>%
    NormalizeData() %>%
    ScaleData(display.progress = FALSE) %>%
    FindVariableGenes(do.plot = FALSE, display.progress = FALSE)
save(CancerCell, file = "./data/CancerCell_20181024.Rda")

#======1.2 QC, pre-processing and normalizing the data=========================
# 1.2.1 Calculate median UMI per cell
Iname = load(file = "./data/CancerCell_20181024.Rda")
CancerCell_raw_data <- as.matrix(x = CancerCell@raw.data)
mean(colSums(CancerCell_raw_data))
median(colSums(CancerCell_raw_data))
min(colSums(CancerCell_raw_data))
remove(CancerCell_raw_data);GC()

# 1.2.3 calculate mitochondria percentage
mito.genes <- grep(pattern = "^MT-", x = rownames(x = CancerCell@data), value = TRUE)
percent.mito <- Matrix::colSums(CancerCell@raw.data[mito.genes, ])/Matrix::colSums(CancerCell@raw.data)
CancerCell <- AddMetaData(object = CancerCell, metadata = percent.mito, col.name = "percent.mito")

CancerCell@ident = factor(CancerCell@ident,levels = samples)

g1 <- VlnPlot(object = CancerCell, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

CancerCell <- FilterCells(object = CancerCell, subset.names = c("nGene","nUMI","percent.mito"),
                    low.thresholds = c(500,1000, -Inf), 
                    high.thresholds = c(Inf,Inf, 0.1))

g2 <- VlnPlot(object = CancerCell, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                    scale_y_log10(limits = c(500,10000)),#+ylim(c(0,1000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                    scale_y_log10(limits = c(500,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                    scale_y_log10(limits = c(1000,100000)),#+ylim(c(0,1000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                    scale_y_log10(limits = c(1000,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                    ylim(c(0,0.2)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,0.2))))
dev.off()
######################################

# After removing unwanted cells from the dataset, the next step is to normalize the data.
CancerCell <- NormalizeData(object = CancerCell, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
CancerCell <- FindVariableGenes(object = CancerCell, mean.function = ExpMean, 
                          dispersion.function = LogVMR, do.plot = FALSE, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(CancerCell@var.genes)
#======1.3 1st run of pca-tsne  =========================
CancerCell <- ScaleData(object = CancerCell) %>%
    RunPCA() %>%
    FindClusters(dims.use = 1:20, force.recalc = T, print.output = FALSE) %>%
    RunTSNE()

p1 <- TSNEPlot(object = CancerCell, do.label = F, group.by = "orig.ident", 
         do.return = TRUE, no.legend = F, #colors.use = singler.colors,
         pt.size = 1,label.size = 8 )+
    ggtitle("Original")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

save(CancerCell, file = "./data/CancerCell_20181024.Rda")
Iname = load("./data/CancerCell_20181024.Rda")

#======1.4 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "./data/seurat_resources/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- HumanGenes(CancerCell,cc.genes[1:43])
g2m.genes <- HumanGenes(CancerCell,cc.genes[44:97])
# Assign Cell-Cycle Scores
CancerCell <- CellCycleScoring(object = CancerCell, s.genes = s.genes, g2m.genes = g2m.genes, 
                         set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = CancerCell, features.plot = HumanGenes(CancerCell,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
# regressing out the difference between the G2M and S phase scores
CancerCell@meta.data$CC.Difference <- CancerCell@meta.data$S.Score - CancerCell@meta.data$G2M.Score
# view cell cycle scores and phase assignments
head(x = CancerCell@meta.data)

#======1.5 Add project id =========================
#------
batchname = CancerCell@meta.data$orig.ident
batch.effect = as.numeric(factor(batchname,levels = samples))
names(batch.effect) = rownames(CancerCell@meta.data)
CancerCell <- AddMetaData(object = CancerCell, metadata = batch.effect, col.name = "batch.effect")
table(CancerCell@meta.data$batch.effect)
head(x = CancerCell@meta.data)

#======1.6 vars.to.regress ScaleData =========================
jpeg(paste0(path,"/S2_CC.Difference.jpeg"), units="in", width=10, height=7,res=600)
SingleFeaturePlot.1(CancerCell,"nUMI",threshold=10000)
SingleFeaturePlot.1(CancerCell,"nGene",threshold=2000)
SingleFeaturePlot.1(CancerCell,"batch.effect",threshold=3.0)
SingleFeaturePlot.1(CancerCell,"percent.mito",threshold=0.05)
SingleFeaturePlot.1(CancerCell,"CC.Difference",threshold=0.05)
dev.off()
CancerCell <- ScaleData(object = CancerCell, 
                  model.use = "linear", do.par=T, do.center = T, do.scale = T,
                  #vars.to.regress = c("nUMI","percent.mito","batch.effect","CC.Difference"),
                  display.progress = T)
#======1.7 Performing MNN-based correction =========================
#https://bioconductor.org/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html#4_performing_mnn-based_correction
set.seed(100)
original <- lapply(samples, function(x) CancerCell@scale.data[CancerCell@var.genes, 
                                                 (CancerCell@meta.data$orig.ident %in% x)])
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.order=T,
                                             approximate=TRUE)))
dim(mnn.out$corrected)
rownames(mnn.out$corrected) = CancerCell@cell.names
colnames(mnn.out$corrected) = paste0("MNN_",1:ncol(mnn.out$corrected))
#Storing a new MNN
CancerCell <- SetDimReduction(object = CancerCell, reduction.type = "MNN", slot = "cell.embeddings",
                       new.data = mnn.out$corrected)
CancerCell <- SetDimReduction(object = CancerCell, reduction.type = "MNN", slot = "key", 
                       new.data = "MNN_")
remove(original);GC()
CancerCell <- SetAllIdent(CancerCell,id = "orig.ident")
DimPlot(object = CancerCell, reduction.use = "MNN", pt.size = 0.5)

#======1.7 unsupervised clustering based on MNN =========================
CancerCell <- RunPCA(object = CancerCell, pc.genes = CancerCell@var.genes, pcs.compute = 100, 
               do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = CancerCell)
PCElbowPlot(object = CancerCell, num.pc = 100)
PCHeatmap(CancerCell, pc.use = c(1:3, 25:30), cells.use = 500, do.balanced = TRUE)

DimElbowPlot.1(object = CancerCell, reduction.type = "MNN", 
             dims.plot = 50,slot = "cell.embeddings")

CancerCell <- RunTSNE(object = CancerCell, reduction.use = "MNN", dims.use = 1:50, 
                do.fast = TRUE, perplexity= 30)

CancerCell <- FindClusters(object = CancerCell, reduction.type = "MNN", 
                    dims.use = 1:50, resolution = 0.6, 
                     k.param = 30,force.recalc = T,
                     save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)

p2 <- TSNEPlot.1(object = CancerCell, do.label = F, group.by = "orig.ident", 
           do.return = TRUE, no.legend = T, 
           pt.size = 1,label.size = 4 )+
    ggtitle("Corrected")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"remove_batch.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1 +theme(legend.position="none"),p2)
dev.off()

jpeg(paste0(path,"tsneplot.jpeg"), units="in", width=10, height=7,res=600)
p2
dev.off()

save(CancerCell, file = "./data/CancerCell_20181024.Rda")
