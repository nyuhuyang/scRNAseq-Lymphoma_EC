library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "./data/Lymphoma_EC_11_20181230.Rda"))
GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();
singler = CreateSinglerObject(as.matrix(object@assays$RNA@data), annot = NULL,project.name = "3119_EC", 
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)

# singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object@assays$RNA@data)
singler$meta.data$orig.ident = CancerCell@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = CancerCell@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = CancerCell@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_Lymphoma_EC_11_20181230.Rda")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./output/singler_CancerCell_20181024.RData")
lnames = load(file = "../SingleR/data/Hpca.RData")
lnames = load(file = "../SingleR/data/Blueprint_encode.RData");lnames

singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
length(singler_colors1);length(singler_colors2)

singler$seurat = CancerCell
SingleR.DrawScatter(sc_data = singler$seurat@data,cell_id = 10, 
                    ref = Hpca, sample_id = 232)

# Step 2: Multiple correlation coefficients per cell types are aggregated 
# to provide a single value per cell type per single-cell. 
# In the examples below we use the 80% percentile of correlation values.
# for visualization purposes we only present a subset of cell types (defined in labels.use)
out = SingleR.DrawBoxPlot(sc_data = singler$seurat@data,cell_id = 10, 
                          ref = immgen,main_types = T,
                          labels.use=c('B cells','T cells','DC','Macrophages','Monocytes','NK cells',
                                       'Mast cells','Neutrophils','Fibroblasts','Endothelial cells'))
print(out$plot)

##############################
# Human Primary Cell Atlas (HPCA)
###############################
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf)
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single.main, top.n = 50))
dev.off()
jpeg(paste0(path,"DrawHeatmap_normF~.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single.main,top.n = 50,normalize = F))
dev.off()
#Next, we can use the fine-tuned labels to color the t-SNE plot:

# HPCA-------
ncolor = length(unique(singler$singler[[1]]$SingleR.single$labels));ncolor
out1 = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single,
                          singler$meta.data$xy,do.label=F,
                          do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                          label.size = 5, dot.size = 2,do.legend = T,alpha = 1,
                          label.repel = F,force=2,
                          title = "Supervised sub cell type labeling by HPCA",
                          colors = singler_colors2[1:ncolor])
jpeg(paste0(path,"/PlotTsne_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
out1
dev.off()
# Blueprint+ENCOD-------
out2 = SingleR.PlotTsne.1(singler$singler[[2]]$SingleR.single,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[2]]$SingleR.single$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=4,
                         title = "Supervised sub cell type labeling by Blueprint+ENCOD",)
#label.repel = T,force=2)
jpeg(paste0(path,"/PlotTsne_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
out2
dev.off()


#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[2]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,
            singler$singler[[1]]$SingleR.single.main$labels)) %>%
        kable_styling()

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub"=singler$singler[[1]]$SingleR.single$labels,
                       "singler1main"=singler$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub"=singler$singler[[2]]$SingleR.single$labels,
                       "singler2main"=singler$singler[[2]]$SingleR.single.main$labels,
                       "cell.names" = rownames(singler$singler[[1]]$SingleR.single$labels))
rownames(singlerDF) = singlerDF[,"cell.names"]
singlerDF = singlerDF[,-which(colnames(singlerDF) == "cell.names")]
#ident.DF[is.na(ident.DF)] <- "unknown"
CancerCell <- AddMetaData(object = CancerCell,
                   metadata = as.data.frame(singlerDF))
CancerCell <- SetAllIdent(object = CancerCell, id = "singler2sub")
##############################
# process singler.color
##############################
apply(singlerDF,2,function(x) length(unique(x)))

CancerCell = AddMetaColor(object = CancerCell, colors = singler_colors1[1:16], 
                         label = "singler2sub")


##############################
# draw tsne plot
##############################
CancerCell <- SetAllIdent(object = CancerCell, id = "singler2sub")
df_colors <- ExtractMetaColor(object = CancerCell)
p3 <- DimPlot.1(object = CancerCell, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, 
          group.by = "ident", do.return = TRUE,
          cols.use = df_colors, no.legend = F,
          do.label =F,label.size=4, label.repel = T,force=1)+
        ggtitle("Supervised sub cell type labeling by Blueprint+ENCOD")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5))

jpeg(paste0(path,"PlotTsne_sub1~~.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(CancerCell,file="./output/CancerCell_20181024.RData")
##############################
# subset Seurat
###############################
table(CancerCell@meta.data$orig.ident)
table(CancerCell@ident)

p4 <- SplitTSNEPlot(CancerCell,group.by = "ident", split.by = "conditions",
                    no.legend = F,do.label =F,label.size=3,
                    return.plots =T, label.repel = T,force=1)
jpeg(paste0(path,"splitetsne~.jpeg"), units="in", width=10, height=7,
     res=600)
do.call(plot_grid, p4[[1]])
dev.off()


p5 <- SplitTSNEPlot(CancerCell,group.by = "ident", split.by = "orig.ident",
                    no.legend = T,do.label =T,label.size=3,
                    return.plots =T, label.repel = T,force=1)
jpeg(paste0(path,"splitetsne.jpeg"), units="in", width=10, height=7,
     res=600)
do.call(plot_grid, p5)
dev.off()
