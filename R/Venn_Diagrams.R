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
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="./output/CancerCell_20181024.RData"))

# 4.1 venn diagram for T cells ==================================
table(CancerCell@ident)
T_cells <- SubsetData(CancerCell, ident.remove = c("Endothelial cells","mv Endothelial cells"))
remove <- FeaturePlot(T_cells,features.plot = "CD3D", do.identify = T)
cells.use <- T_cells@cell.names[!(T_cells@cell.names %in% remove)]
T_cells <- SubsetData(CancerCell, cells.use = cells.use)

T_cells <- SetAllIdent(T_cells, id = "orig.ident")
#T_cells <- RenameIdent(T_cells, "TALL-plus-EC-2","TALL-plus-EC")
T_cells_markers <- FindAllMarkers.UMI(T_cells, return.thresh = 0.05,
                                      test.use = "MAST",min.pct = -Inf,
                                      min.cells.gene = -Inf, min.cells.group = -Inf,
                                      only.pos = T, logfc.threshold = -Inf)
write.csv(T_cells_markers,paste0(path,"T_cells_markers1.csv"))

eulerr.1 <- function(df, shape =  "circle", split.by = NULL,
                     cut_off = "avg_logFC", cut_off_value = 0.05){
        df$cluster <- as.vector(df$cluster)
        if(!is.null(split.by)) df <- df[(df$cluster %in% split.by),]
        df_list <- split(df,df$cluster)
        
        if(cut_off == "avg_logFC"){
                pos.share_genes <- sapply(df_list, function(df) df[(df$avg_logFC > -cut_off_value),"gene"])
        }  
        if(any(cut_off %in% c("p_val","p_val_adj"))){
                pos_genes <- lapply(df_list, function(df) df[(df$avg_logFC > 0),"gene"])
                shared_genes <- lapply(df_list, function(df) df[(abs(df[,cut_off]) > cut_off_value),"gene"])
                pos.share_genes <- mapply(function(x,y) unique(c(x, y)), pos_genes, shared_genes)
        }
                return(eulerr::euler(pos.share_genes,shape = shape))
}

##-----------------------------------
jpeg(paste0(path,"/Venn_T_TE_pval.jpeg"), units="in", width=10, height=7,res=600)
p_val = c(0.5,0.1,0.05,0.01)
g <- list()
for(i in 1:length(p_val)){
        fit <- eulerr.1(T_cells_markers,split.by = c("TALL-only", "TALL-plus-EC-2"),
                        cut_off = "p_val", cut_off_value = p_val[i])
        g[[i]] <- plot(fit, quantities = TRUE, fill = c("#377EB8","#A6D854"),lty = 1:2,
                       labels = list(labels = c("TALL", "TALL\n+EC")), legend = F,
                       main = paste("p_val:", p_val[i]))
}
do.call(cowplot::plot_grid, c(g,ncol =2))
dev.off()

##-----------------------------------
jpeg(paste0(path,"/Venn_T_TE_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
cut_off_values <- c(c(0.1,0.05,0.01),c(0.99,0.5,0.1),c(0.1,0.05,0.01))
names(cut_off_values) <- c(rep("p_val",3),rep("p_val_adj",3),rep("avg_logFC",3))
g <- list()
for(i in 1:length(cut_off_values)){
        fit <- eulerr.1(T_cells_markers, split.by = c("TALL-only", "TALL-plus-EC"),
                        cut_off = names(cut_off_values[i]),
                        cut_off_value = cut_off_values[i])
        g[[i]] <- plot(fit, quantities = TRUE, fill = c("#377EB8","#A6D854"),lty = 1:2,
                       main = paste(names(cut_off_values[i]),":",cut_off_values[i]))
}
do.call(cowplot::plot_grid, c(g,ncol =3))
dev.off()

##-----------------------------------
jpeg(paste0(path,"/Venn_T_3_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
colors_list <- list(c("#377EB8","#A6D854"),c("#377EB8","#E41A1C"),
                    c("#A6D854","#E41A1C"))
labels_list <- list(c("TALL-only", "TALL-plus-EC"),c("TALL-only", "TALL-plus-EC-2"),
                    c("TALL-plus-EC","TALL-plus-EC-2"))
cut_off <- c("p_val","p_val_adj","avg_logFC")
fit <- list(); g <- list()
m=1
for(j in 1:length(cut_off)){
        for(i in 1:length(labels_list)){
                fit[[i]] <- eulerr.1(T_cells_markers, split.by = labels_list[[i]],
                                     cut_off = cut_off[j],
                                     cut_off_value = 0.05)
                g[[m]] <- plot(fit[[i]], quantities = TRUE, lty = 1:2,
                               fill = colors_list[[i]],
                               legend = F,main =  paste(cut_off[j],": 0.05"))
                m <- m+1
        }    
}
do.call(cowplot::plot_grid, c(g,ncol =3))
dev.off()

##-----------------------------------
jpeg(paste0(path,"/Venn_T_all_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
g <- lapply(c("p_val","p_val_adj","avg_logFC"), function(cut_off) {
        eulerr.1(T_cells_markers, cut_off = cut_off, cut_off_value = 0.05) %>%
                plot(quantities = TRUE, lty = 1:2,
                     fill = c("#377EB8","#A6D854","#E41A1C"),
                     labels = list(labels = c("TALL", "TALL+EC2","\n\n\nTAL+EC")), 
                     legend = F,main = paste(cut_off,": 0.05"))
})
do.call(cowplot::plot_grid, c(g,ncol =2))
dev.off()

##-----------------------------------
jpeg(paste0(path,"/Venn_T_all_p_val_adj.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr.1(T_cells_markers,
              cut_off = "p_val_adj", cut_off_value = 0.05) %>%
        plot(quantities = TRUE, lty = 1:3,
             fill = c("#377EB8","#A6D854","#E41A1C"),
             legend = T,main = paste("p_val_adj : 0.05"))
print(g)
dev.off()

##-----------------------------------
jpeg(paste0(path,"/Venn_T_all_logFC.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr.1(T_cells_markers,
              cut_off = "avg_logFC", cut_off_value = 0.01) %>%
                plot(quantities = TRUE, lty = 1:3,
                     fill = c("#377EB8","#A6D854","#E41A1C"),
                     legend = T,main = paste("avg_logFC: 0.01"))
print(g)
dev.off()

# 4.2 venn diagram for EC cells ==================================
EC <- SubsetData(CancerCell, ident.use = c("Endothelial cells","mv Endothelial cells"))
SplitTSNEPlot(EC, do.print = T, do.return = F)
EC <- SetAllIdent(EC, id = "orig.ident")
#EC <- RenameIdent(EC, "TALL-plus-EC-2","TALL-plus-EC")
EC@meta.data$orig.ident <- gsub("TALL-plus-EC-2","EC2-plus-TALL",EC@meta.data$orig.ident)
EC@meta.data$orig.ident <- gsub("TALL-plus-EC","EC-plus-TALL",EC@meta.data$orig.ident)

EC_markers <- FindAllMarkers.UMI(EC,return.thresh = 0.05,
                                      test.use = "MAST",min.pct = -Inf,
                                      min.cells.gene = -Inf, min.cells.group = -Inf,
                                      only.pos = T, logfc.threshold = -Inf)
write.csv(EC_markers,paste0(path,"EC_markers1.csv"))

#-----
jpeg(paste0(path,"/Venn_EC_all_p_val_adj.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr.1(EC_markers,
              cut_off = "p_val_adj", cut_off_value = 0.05) %>%
        plot(quantities = TRUE, lty = 1:3,
             fill = c("#377EB8","#A6D854","#E41A1C"),
             legend = T,main = paste("p_val_adj: 0.05"))
print(g)
dev.off()

#-----
jpeg(paste0(path,"/Venn_EC_all_p_val.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr.1(EC_markers,
              cut_off = "p_val", cut_off_value = 0.05) %>%
        plot(quantities = TRUE, lty = 1:3,
             fill = c("#377EB8","#A6D854","#E41A1C"),
             legend = T,main = paste("p_val: 0.05"))
print(g)
dev.off()

#-----
jpeg(paste0(path,"/Venn_EC_all_logFC.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr.1(EC_markers,
              cut_off = "avg_logFC", cut_off_value = 0.01) %>%
        plot(quantities = TRUE, lty = 1:3,
             fill = c("#377EB8","#A6D854","#E41A1C"),
             legend = T,main = paste("avg_logFC: 0.01"))
print(g)
dev.off()

# 4.2 venn diagram for T_EC cells ==================================
T_EC <- MergeSeurat(T_cells, EC)
T_EC <- SetAllIdent(T_EC, id = "orig.ident")
table(T_EC@ident)

T_EC_markers <- FindAllMarkers.UMI(T_EC,return.thresh = 1,
                                 test.use = "MAST",min.pct = -Inf,
                                 min.cells.gene = -Inf, min.cells.group = -Inf,
                                 only.pos = F, logfc.threshold = -Inf)
write.csv(T_EC_markers,paste0(path,"T_EC_markers.csv"))

#-----
jpeg(paste0(path,"/Venn_TEC_all_p_val.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr.1(T_EC_markers,shape =  "ellipse",
              cut_off = "p_val", cut_off_value = 0.05) %>%
        plot(quantities = TRUE, lty = 1:6,
             fill = c("#1e601a","#A6D854","#FDC086",
                      "#377EB8","#E41A1C","#6A3D9A"),
             legend = T,main = paste("p_val: 0.05"))
print(g)
dev.off()

#-----
jpeg(paste0(path,"/Venn_TEC_all_p_val_adj.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr.1(T_EC_markers,shape =  "ellipse",
              cut_off = "p_val_adj", cut_off_value = 0.05) %>%
        plot(quantities = TRUE, lty = 1:6,
             fill = c("#1e601a","#A6D854","#FDC086",
                      "#377EB8","#E41A1C","#6A3D9A"),
             legend = T,main = paste("p_val_adj: 0.05"))
print(g)
dev.off()

#-----
jpeg(paste0(path,"/Venn_TEC_all_avg_logFC.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr.1(T_EC_markers,shape =  "ellipse",
              cut_off = "avg_logFC", cut_off_value = 0.01) %>%
        plot(quantities = TRUE, lty = 1:6,
             fill = c("#1e601a","#A6D854","#FDC086",
                      "#377EB8","#E41A1C","#6A3D9A"),
             legend = T,main = paste("avg_logFC: 0.01"))
print(g)
dev.off()

