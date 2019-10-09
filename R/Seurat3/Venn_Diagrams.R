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
library(magrittr)
library(tidyr)
library(ggplot2)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 3.1.1 load data
# test normalization methods
dates <- c("2018-10-18","2018-12-30")
cluster <- c("EC+RO2", "EC+3119")
cell.type = c("EC","TALL")
k=2
TEC <- list()
for(i in seq_along(c("2018-10-18","2018-12-30"))){
        TEC[[i]] = read.csv(file=paste0("Yang/6. Differential analysis/",cell.type[k],
                                        "/",cell.type[k],"_markers_",dates[i],"_RNA.csv"),
                                row.names = 1,stringsAsFactors=F)
        #TEC[[i]] = TEC[[i]][TEC[[i]]$cluster %in% "TEC",]
        TEC[[i]]$gene = rownames(TEC[[i]])
        TEC[[i]]$cluster = cluster[i]
}

TEC = bind_rows(TEC)
g_ec <- mapply(function(x,y) eulerr(TEC, cut_off = x, cut_off_value = y,
                                    do.print = F,do.lenged = F), 
               c("p_val","p_val_adj","avg_logFC"),
               c(0.05,0.05,0.1),SIMPLIFY = FALSE)

jpeg(paste0(path,"Venn_",cell.type[k],"_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
do.call(cowplot::plot_grid, c(g_ec,ncol =2, scale = 0.8))
dev.off()
#========
# test EC only
assay = "SCT"
slot = "data"
cell.type = "EC"
cluster <- c("RO2-plus-EC-2","RO2-plus-EC","EC2","EC+3119")
(dates <- c(paste("Endothelial cells",cluster[1:3],"vs_EC",sep = "_"),"2018-12-30"))
TEC <- list()
for(i in seq_along(cluster)){
        assay_slot = paste0(path,assay,"_",slot,"/")
        TEC[[i]] = read.csv(file = paste0("Yang/6. Differential analysis/",
                                          basename(assay_slot),"/", cell.type,"_markers_",
                                          dates[i],"_",basename(assay_slot),".csv"),
                            row.names = 1,stringsAsFactors=F)
        TEC[[i]]$gene = rownames(TEC[[i]])
        TEC[[i]]$cluster = cluster[i]
}

TEC = bind_rows(TEC)
g_ec <- mapply(function(x,y) eulerr(TEC, cut_off = x, cut_off_value = y,
                                    do.print = F,do.lenged = T), 
               c("p_val","p_val_adj","avg_logFC"),
               c(0.1,0.5,0.01),SIMPLIFY = FALSE)

jpeg(paste0(path,"Venn_",paste(cluster,collapse = "_"),"L.jpeg"), units="in", width=10, height=7,res=600)
do.call(cowplot::plot_grid, c(g_ec,ncol =2, scale = 0.8))
dev.off()

#==== correlation ==== 
TEC <- list()
for(i in seq_along(cluster)){
        assay_slot = paste0(path,assay,"_",slot,"/")
        TEC[[i]] = read.csv(file = paste0("Yang/6. Differential analysis/",
                                          basename(assay_slot),"/", cell.type,"_markers_",
                                          dates[i],"_",basename(assay_slot),".csv"),
                            row.names = 1,stringsAsFactors=F)
        TEC[[i]]$cluster = cluster[i]
        TEC[[i]] = TEC[[i]][,-grep("p_val|pct.1|pct.2|p_val_adj|avg_UMI.1|avg_UMI.2",colnames(TEC[[i]]))]
}
df_TEC = bind_cols(TEC)
df_TEC = df_TEC[,-grep("cluster",colnames(df_TEC))]
colnames(df_TEC) = cluster
rownames(df_TEC)= rownames(TEC[[i]])
