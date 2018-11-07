########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(SingleR)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
lnames = load(file="./output/CancerCell_20181024.RData");lnames
lnames = load(file = "./output/singler_CancerCell_20181024.RData")
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
length(singler_colors1);length(singler_colors2)

merge.ident <- gsub("TALL-plus-EC-2","TALL-plus-EC",CancerCell@meta.data$orig.ident)
merge.identDF = data.frame("merge.ident" = merge.ident,
                           row.names = rownames(singler$singler[[1]]$SingleR.single$labels))
CancerCell <- AddMetaData(object = CancerCell, metadata = merge.identDF)


Counts.Barplot <- function(labels, samples, file.name, mar = c(5, 4, 4, 12), 
                           main="Total numbers of cell types in each sample",
                           col = singler.colors, args.legend = NULL){
        path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        counts <- table(as.vector(labels), samples)
        counts <- counts[rev(sort(unique(labels))),]
        jpeg(paste0(path,file.name), units="in", width=10, height=7,
             res=600)
        par(mfrow=c(1, 1), mar = mar)
        barplot(counts, main=main,
                xlab="Samples", ylab="Cell numbers",
                col= rev(col),
                legend = rownames(counts),
                args.legend = args.legend)
        dev.off()
}
#----------- CancerCell@meta.data$singler2sub -------------
Counts.Barplot(labels = CancerCell@meta.data$singler2sub,
               samples = CancerCell@meta.data$merge.ident, 
               file.name = "1_bar_chart.jpeg", col = singler_colors1[1:16], 
               args.legend = list(x = "topright", bty = "n",inset=c(-0.3, 0)))
#----------- CancerCell@meta.data$singler2sub -------------

Counts.Barplot(labels = CancerCell@meta.data$singler1sub,
               samples = CancerCell@meta.data$merge.ident, mar=c(5, 4, 8, 20),
               file.name = "1_bar_chart～.jpeg", col = singler_colors2[1:34], 
               args.legend = list(x = "topright", bty = "n",inset=c(-0.8, -0.35)))


counts <- table(CancerCell@meta.data$singler1sub, CancerCell@meta.data$merge.ident)
kable(counts) %>% kable_styling()

#-pathway tsneplot-----------------------------------------------------------------------
h.all = read.delim("../seurat_resources/h.all.v6.2.symbols.gmt",row.names =1,header = F,
                   stringsAsFactors = F)
h.all <- h.all[,-1]
h.all.df <- as.data.frame(t(h.all))
h.all.list <- df2list(h.all.df)

PI3K_AKT <- h.all.list[["HALLMARK_PI3K_AKT_MTOR_SIGNALING"]]

CancerCell <-AddModuleScore(CancerCell, genes.list = list(PI3K_AKT),
                            ctrl.size = 5,enrich.name = "PI3K-AKT")
min_score =  min(CancerCell@meta.data$`PI3K-AKT1`)
CancerCell@meta.data$`PI3K-AKT1` = CancerCell@meta.data$`PI3K-AKT1` - min_score
CancerCell_list <- SplitSeurat(CancerCell, split.by = "merge.ident")

p <- list()
levels <- CancerCell_list[[length(CancerCell_list)]];levels
for(i in 1:length(levels)){
        p[[i]] <- SingleFeaturePlot.1(object = CancerCell_list[[i]],
                                      feature = "PI3K.AKT1",
                                      title=levels[i], threshold = 0.25)
}

jpeg(paste0(path,"PI3K-AKT.jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid, p))
dev.off()

# Differential analysis for EC==================================
g <- SplitTSNEPlot(CancerCell, split.by = "merge.ident")
g1 <- TSNEPlot.1(CancerCell, no.legend = F, colors.use = ExtractMetaColor(object = CancerCell))
jpeg(paste0(path,"SplitTSNEPlot~.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g[[1]],g[[3]]))
dev.off()

table(CancerCell@ident)
EC <- SubsetData(CancerCell, ident.use = "Endothelial cells")
SplitTSNEPlot(EC, split.by = "merge.ident",select.plots = c(1,3), do.print = T)
EC <- SetAllIdent(EC, id = "merge.ident")

EC_markers <- FindAllMarkers.UMI(EC, test.use = "MAST",
                                 only.pos = FALSE, logfc.threshold = 0.1)
EC_markers %>% head(20) %>% kable() %>% kable_styling()
top <-  EC_markers %>% group_by(cluster) %>% top_n(200, avg_logFC)
write.csv(top[!duplicated(top$gene),],paste0(path,"EC_markers_singler2sub.csv"))

EC_samples <- c("EC-only","TALL-plus-EC")
top <- group_top_mutate(df = EC_markers, EC_samples)
top %>% head(20) %>% kable() %>% kable_styling()
g2 <- DoHeatmap.1(EC,top,Top_n = 50, 
                    group.order = EC_samples,ident.use = "Endothelial cells",
                    group.label.rot = F,cex.row = 5,remove.key =F)
jpeg(paste0(path,"/heatmap_EC.jpeg"), units="in", width=10, height=7,res=600)
g2
dev.off()

# 4.1.3 GO enrichment analysis ==============
table(abs(top$avg_logFC) > 0.1)
top_markers <- top[abs(top$avg_logFC) > 0.1,"gene"]
length(unique(top_markers))
library(biomaRt)
# collect gene names from biomart
mart <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
# Get ensembl gene ids and GO terms
getBM <- biomaRt::getBM(attributes = c("external_gene_name", "go_id"),
                        filters= "external_gene_name",
                        values = rownames(EC@data),
                        mart = mart)

EC.Go <-  topGOterms(int.genes = unique(top_markers),
                     bg.genes = rownames(EC@data),
                     organism =  "Human",
                     ontology = c("BP","CC","MF"),
                     getBM = getBM)
EC.Go %>% kable() %>% kable_styling()
write.csv(EC.Go,paste0(path,"EC_go_analysis.csv"))


# Differential analysis for T cells ==================================
g <- SplitTSNEPlot(CancerCell, split.by = "merge.ident")
g1 <- TSNEPlot.1(CancerCell, no.legend = F, colors.use = ExtractMetaColor(object = CancerCell))
jpeg(paste0(path,"SplitTSNEPlot2_T.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g[[2]],g[[3]]))
dev.off()

CancerCell <- SetAllIdent(CancerCell, id = "singler2sub")
not_T_cells <- c("DC","Endothelial cells","Erythrocytes","Macrophages",
                 "Macrophages M1","Macrophages M2","mv Endothelial cells")
T_cells <- SubsetData(CancerCell, ident.remove = not_T_cells)
FeaturePlot(T_cells,features.plot = "CD3D", do.identify = T)
SplitTSNEPlot(T_cells, split.by = "merge.ident", do.print = T)
remove <- FeaturePlot(T_cells,features.plot = "EMCN", do.identify = T)
cells.use <- T_cells@cell.names[!(T_cells@cell.names %in% remove)]
T_cells <- SubsetData(CancerCell, cells.use = cells.use)

T_cells <- SetAllIdent(T_cells, id = "merge.ident")
T_cells_markers <- FindAllMarkers.UMI(T_cells, test.use = "MAST",
                                 only.pos = FALSE, logfc.threshold = 0.1)
T_cells_markers %>% head(20) %>% kable() %>% kable_styling()
top <-  T_cells_markers %>% group_by(cluster) %>% top_n(200, avg_logFC)
write.csv(top[!duplicated(top$gene),],paste0(path,"T_cells_markers_singler2sub.csv"))

T_cells_samples <- c("TALL-only","TALL-plus-EC")
top <- group_top_mutate(df = T_cells_markers, T_cells_samples)
top %>% head(20) %>% kable() %>% kable_styling()
g3 <- DoHeatmap.1(T_cells,top,Top_n = 50, 
                  group.order = T_cells_samples,ident.use = "T cells",
                  group.label.rot = F,cex.row = 5,remove.key =F)
jpeg(paste0(path,"/heatmap_T_cells.jpeg"), units="in", width=10, height=7,res=600)
g3
dev.off()

# 4.1.3 GO enrichment analysis ==============
table(abs(top$avg_logFC) > 0.2)
top_markers <- top[abs(top$avg_logFC) > 0.2,"gene"]
length(unique(top_markers))
library(biomaRt)
# collect gene names from biomart
mart <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
# Get ensembl gene ids and GO terms
getBM <- biomaRt::getBM(attributes = c("external_gene_name", "go_id"),
                        filters= "external_gene_name",
                        values = rownames(T_cells@data),
                        mart = mart)

T_cells.Go <-  topGOterms(int.genes = unique(top_markers),
                     bg.genes = rownames(T_cells@data),
                     organism =  "Human",
                     ontology = c("BP","CC","MF"),
                     getBM = getBM)
T_cells.Go %>% kable() %>% kable_styling()
write.csv(T_cells.Go,paste0(path,"T_cells_go_analysis.csv"))
