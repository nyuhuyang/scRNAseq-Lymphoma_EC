########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(SingleR)
library(kableExtra)
library(eulerr)
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
                           row.names = rownames(CancerCell@meta.data))
CancerCell <- AddMetaData(object = CancerCell, metadata = merge.identDF)


markers  <- HumanGenes(CancerCell,c("CD59","IGFBP7","IGF2R","HMGB1"))
df <- data.frame("gene" = c("IGF2R"),
                 "alias" = c("CD222","CD59"))
SplitSingleFeaturePlot(CancerCell,group.by = "ident",split.by = "merge.ident",
                       no.legend = T, alias = Alias(df = df, gene = markers),
                       label.size=3,do.print =T,markers = markers,
                       threshold = NULL,ncol=2)
# grouped bar chart
# create a dataset
EC <- SubsetData(CancerCell, ident.use = c("Endothelial cells","mv Endothelial cells"))
EC <- SetAllIdent(EC, id = "orig.ident")

UMI.1 <- Matrix::rowMeans(x = EC@raw.data[, WhichCells(object = EC,ident = "EC-only")])
UMI.2 <- Matrix::rowMeans(x = EC@raw.data[, WhichCells(object = EC,ident = c("TALL-plus-EC-2",
                                                                             "TALL-plus-EC"))])
sample=c(rep("EC-only" , 4) , rep("TALL-plus_EC" , 4) )
genes=rep(markers , 2)
value=c(UMI.1[markers],UMI.2[markers])
data=data.frame(sample,genes,value)

# Grouped
jpeg(paste0(path,"/EC_group_barchart.jpeg"), units="in", width=10, height=7,res=600)
ggplot(data, aes(fill=genes, y=value, x=sample)) + 
        geom_bar(position="dodge", stat="identity")+
        ylab("UMI")
dev.off()


# Counts.Barplot
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

# Differential analysis for EC ==================================
CancerCell = AddMetaColor(object = CancerCell, colors = singler_colors1[1:16], 
                          label = "singler2sub")
g <- SplitTSNEPlot(CancerCell, split.by = "orig.ident",do.label = F,force =3)
g1 <- TSNEPlot.1(CancerCell, no.legend = F, colors.use = ExtractMetaColor(object = CancerCell))
jpeg(paste0(path,"SplitTSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid,g))
dev.off()

table(CancerCell@ident)
EC <- SubsetData(CancerCell, ident.use = c("Endothelial cells","mv Endothelial cells"))
SplitTSNEPlot(EC, split.by = "orig.ident",do.print = T)
EC <- SetAllIdent(EC, id = "orig.ident")

EC_markers <- FindAllMarkers.UMI(EC, test.use = "MAST",
                                 only.pos = FALSE, logfc.threshold = 0.1)
EC_markers %>% head(20) %>% kable() %>% kable_styling()
top <-  EC_markers %>% group_by(cluster) %>% top_n(300, avg_logFC)
write.csv(top[!duplicated(top$gene),],paste0(path,"EC_markers_singler2sub.csv"))

EC_samples <- c("EC-only","TALL-plus-EC","TALL-plus-EC-2")
top <- group_top_mutate(df = EC_markers, EC_samples)
top %>% head(20) %>% kable() %>% kable_styling()
g2 <- DoHeatmap.1(EC,top,Top_n = 40, 
                    group.order = EC_samples,ident.use = "Endothelial cells",
                    group.label.rot = T,cex.row = 4,remove.key =F)
jpeg(paste0(path,"/heatmap_EC.jpeg"), units="in", width=10, height=7,res=600)
g2
dev.off()
# make color bar ========
MakeCorlorBar(EC, top,Top_n=40)

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

# 4.1.4 Gene analysis ==============
subset.EC <- SplitSeurat(EC, split.by = "orig.ident")
(lvl <- subset.EC[[length(subset.EC)]])

markers <- HumanGenes(EC, c("VWF"))#,"IFI27","HLA-B","IFI6"))##,"CD3D"))
for(marker in markers){
        g <- list()
        g[[1]] <- SingleFeaturePlot.1(object = subset.EC[[1]], 
                                      feature = marker,title = lvl[1],threshold=3)
        g[[2]] <- SingleFeaturePlot.1(object = subset.EC[[2]], 
                                      feature = marker,title = lvl[2],threshold=3)
        g[[3]] <- SingleFeaturePlot.1(object = subset.EC[[3]], 
                                      feature = marker,title = lvl[3],threshold=3)
        jpeg(paste0(path,"Splited_EC_",marker,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        print(paste0(which(markers == marker),":",length(markers)))
        dev.off()
}
# 4.1.5 Gene sets analysis ==============
for(j in 1:length(gene_sets)){
        p <- list()
        for(i in 1:length(lvl)){
                p[[i]] <- SingleFeaturePlot.1(object = subset.EC[[i]],
                                              feature = gene_sets_names[j],
                                              title=lvl[i], threshold = 0.1)
        }
        jpeg(paste0(path,gene_sets_names[j],".jpeg"), units="in", width=10, height=7,res=600)
        print(do.call(plot_grid, p))
        dev.off()
}
# 4.1.6 cell cycle analysis ==============
df1 <- cbind.data.frame(subset.EC[[1]]@meta.data[,c("G2M.Score","S.Score")],
                        "sample" = rep(lvl[1],length(subset.EC[[1]]@cell.names)))
df2 <- cbind.data.frame(subset.EC[[2]]@meta.data[,c("G2M.Score","S.Score")],
                        "sample" = rep(lvl[2],length(subset.EC[[2]]@cell.names)))
df3 <- cbind.data.frame(subset.EC[[3]]@meta.data[,c("G2M.Score","S.Score")],
                        "sample" = rep(lvl[3],length(subset.EC[[3]]@cell.names)))

cell_cycle <- rbind.data.frame(df1,df2,df3)
jpeg(paste0(path,"EC_cell_cycle.jpeg"), units="in", width=10, height=7,
     res=600)
ggplot(cell_cycle, aes(G2M.Score, S.Score, color= sample)) +geom_point()
dev.off()

# 4.2 Differential analysis for T cells ==================================
table(CancerCell@ident)
T_cells <- SubsetData(CancerCell, ident.remove = c("Endothelial cells","mv Endothelial cells"))
remove <- FeaturePlot(T_cells,features.plot = "CD3D", do.identify = T)
cells.use <- T_cells@cell.names[!(T_cells@cell.names %in% remove)]
T_cells <- SubsetData(CancerCell, cells.use = cells.use)
SplitTSNEPlot(T_cells, split.by = "orig.ident", do.print = F,do.return =F)

T_cells <- SetAllIdent(T_cells, id = "orig.ident")
T_cells_markers <- FindAllMarkers.UMI(T_cells, test.use = "MAST",
                                 only.pos = FALSE, logfc.threshold = 0.1)
T_cells_markers %>% head(20) %>% kable() %>% kable_styling()
top <-  T_cells_markers %>% group_by(cluster) %>% top_n(200, avg_logFC)
write.csv(top[!duplicated(top$gene),],paste0(path,"T_cells_markers_singler2sub.csv"))

T_cells3 <- SubsetData(T_cells, ident.remove = "TALL-plus-EC")
T_cells3_markers <- FindAllMarkers.UMI(T_cells3, test.use = "MAST",
                                      only.pos = T, logfc.threshold = 0.1)

T_cells_samples <- c("TALL-only","TALL-plus-EC-2")
top <- group_top_mutate(df = T_cells3_markers, T_cells_samples)
top %>% tail(20) %>% kable() %>% kable_styling()
write.csv(top[!duplicated(top$gene),],paste0(path,"T_cells_A_C_markers.csv"))

g4 <- DoHeatmap.1(T_cells1,top,Top_n = 40, 
                  group.order = T_cells_samples,ident.use = "T cells",
                  group.label.rot = T,cex.row = 4,remove.key =F)
jpeg(paste0(path,"/heatmap_T_cells_EC.jpeg"), units="in", width=10, height=7,res=600)
g4
dev.off()

MakeCorlorBar(T_cells, top,Top_n=40)

# 4.2.1 venn diagram for T cells ==================================

# 4.2.3 GO enrichment analysis ==============
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

# 4.2.4 Gene analysis ==============
T_cells@meta.data$orig.ident <- gsub("TALL-plus-EC-2","TALL-plus-EC",
                                     T_cells@meta.data$orig.ident)
subset.T_cells <- SplitSeurat(T_cells, split.by = "orig.ident")
(lvl <- subset.T_cells[[length(subset.T_cells)]])

markers <- HumanGenes(EC, c("IFI27","HES4","CD7","RUNX3"))
for(marker in markers){
        g <- list()
        g[[1]] <- SingleFeaturePlot.1(object = subset.T_cells[[1]], 
                                      feature = marker,title = lvl[1],threshold=0.1)
        g[[2]] <- SingleFeaturePlot.1(object = subset.T_cells[[2]], 
                                      feature = marker,title = lvl[2],threshold=0.1)
        g[[3]] <- SingleFeaturePlot.1(object = subset.T_cells[[3]], 
                                      feature = marker,title = lvl[3],threshold=0.1)
        jpeg(paste0(path,"Splited_T.cells_",marker,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        print(paste0(which(markers == marker),":",length(markers)))
        dev.off()
}

# 4.1.5 Gene sets analysis ==============
for(j in 1:length(gene_sets)){
        p <- list()
        for(i in 1:length(lvl)){
                p[[i]] <- SingleFeaturePlot.1(object = subset.T_cells[[i]],
                                              feature = gene_sets_names[j],
                                              title=lvl[i], threshold = 0.1)
        }
        jpeg(paste0(path,gene_sets_names[j],"_Tcells.jpeg"), units="in", width=10, height=7,res=600)
        print(do.call(plot_grid, p))
        dev.off()
}

# 4.2.6 cell cycle analysis ==============
df1 <- cbind.data.frame(subset.T_cells[[1]]@meta.data[,c("G2M.Score","S.Score")],
                        "sample" = rep(lvl[1],length(subset.T_cells[[1]]@cell.names)))
df2 <- cbind.data.frame(subset.T_cells[[2]]@meta.data[,c("G2M.Score","S.Score")],
                        "sample" = rep(lvl[2],length(subset.T_cells[[2]]@cell.names)))
df3 <- cbind.data.frame(subset.T_cells[[3]]@meta.data[,c("G2M.Score","S.Score")],
                        "sample" = rep(lvl[3],length(subset.T_cells[[3]]@cell.names)))

cell_cycle <- rbind.data.frame(df1,df2,df3)
jpeg(paste0(path,"T_cell_cycle.jpeg"), units="in", width=10, height=7,
     res=600)
ggplot(cell_cycle, aes(G2M.Score, S.Score, color= sample)) +geom_point()
dev.off()



data.use <- T_cells@data["HES4",] #%>% t() %>% #scale(center = F) %>%
        t() %>% as.data.frame() %>% gather(key = Subtypes.markers, value = ave.expr)
jpeg(paste0(path,"/density~.jpeg"), units="in", width=10, height=7,res=600)
ggplot(data.use, aes(x = ave.expr, fill = Subtypes.markers)) +
        geom_density(alpha = .5) + scale_y_sqrt() +
        ggtitle("Normalized lineage scores in human muscle-invasive bladder cancer")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 15, face = "bold")) 
dev.off()


# 4.3 vlnplot ==================================
marker.genes <- HumanGenes(T_EC,c("JARID2","MTF2","EZH2", "EED", "SUZ12",
                                  "CRIP2","TRAP1"))
g1 <- lapply(marker.genes, function(marker.gene){ 
        VlnPlot.1(T_EC, features.plot = marker.gene, nCol = 1,
                  point.size.use = 0.001,x.lab = "",
                  x.lab.rot = T, do.return = T)+
                coord_flip()+
                theme(axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank())
})
jpeg(paste0(path,"/PRC2.jpeg"), units="in", width=10, height=7,res=600)
do.call(plot_grid, c(g1,nrow = 1,align="hv"))
dev.off()

# 4.4 heatmap ==================================
markers <- HumanGenes(T_EC,c("TP53","CDNK2A","EZH2", "EED", "SUZ12","CRIP2",
                             "TRAP1","BCL2L11","PTEN","NOTCH1", "FBXW7"))
samples <- c("EC-only","EC-plus-TALL","EC2-plus-TALL",
           "TALL-only","TALL-plus-EC","TALL-plus-EC-2")
T_EC_markers1 <- T_EC_markers[c(16973:50915,1:16972,
                                50916:67888,84861:101822, 67889:84860),]
CancerCell <- SetAllIdent(CancerCell, id = "old.ident")
CancerCell@ident <- factor(CancerCell@ident, levels = samples)
table(CancerCell@ident)

jpeg(paste0(path,"/heatmap.jpeg"), units="in", width=10, height=7,res=600)
DoHeatmap.1(CancerCell, T_EC_markers1, add.genes = markers,cex.row = 3,
            Top_n = 25, ident.use = "TALL and EC")
dev.off()



# compare EC 2019-01-25
object@scale.data = readRDS("./data/Lymphoma_EC.scale.data_Harmony_12_20190125.rds")
object <- SetAllIdent(object, id='orig.ident')
EC_only <- SubsetData(object, ident.use = c("EC","EC-only"))

Angiocrine_factors <- read.delim("data/seurat_resources/Angiocrine_factors.txt",
                                 stringsAsFactors =F)
Angiocrine_factors <- HumanGenes(EC_only, Angiocrine_factors[,1])
blueprint_encode <- read.csv("../SingleR/output/blueprint_encode_main.csv", stringsAsFactors =F)
EC_mark <- HumanGenes(EC_only,c(Angiocrine_factors,blueprint_encode$Endothelial_cells), unique = T)

y = EC_only@scale.data[EC_mark,]
library(gplots)

hc <- hclust(as.dist(1-cor(as.matrix(y), method="spearman")), method="complete")
cc = gsub("_.*","",hc$labels)
cc = gsub("EC-only","#195016",cc)
cc = gsub("EC","#B3DE69",cc)
table(cc)
jpeg(paste0(path,"/Heatmap2_EC.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          Colv = as.dendrogram(hc), Rowv= FALSE,
          ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",
          key.xlab = "nUMI",
          cexRow = 0.001,
          margins = c(2,5),
          #scale = "row",
          breaks = seq(-3,3,length.out = 101),
          col = bluered,
          main = paste("compare EC cells only"))
par(lend = 1)           # square line ends for the color legend
legend(0, 0.8,       # location of the legend on the heatmap plot
       legend = c("EC-only", "EC"), # category labels
       col = c("#195016", "#B3DE69"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
          
dev.off()
