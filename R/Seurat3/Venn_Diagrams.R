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
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 3.1.1 load data
dates <- c("2018-10-18","2018-12-30")
cluster <- c("EC+RO2", "EC+3119")
cell.type = c("EC","TALL")
k=2
TEC <- list()
for(i in seq_along(c("2018-10-18","2018-12-30"))){
        TEC[[i]] = read.csv(file=paste0("output/20190926/",cell.type[k],"_markers_",dates[i],".csv"),
                                row.names = 1,stringsAsFactors=F)
        #TEC[[i]] = TEC[[i]][TEC[[i]]$cluster %in% "TEC",]
        TEC[[i]]$gene = rownames(TEC[[i]])
        TEC[[i]]$cluster = cluster[i]
}

TEC = bind_rows(TEC)
g_ec <- mapply(function(x,y) eulerr(TEC, cut_off = x, cut_off_value = y,
                                    do.print = T,do.lenged = F), 
               c("p_val","p_val_adj","avg_logFC"),
               c(0.05,0.05,0.1),SIMPLIFY = FALSE)

jpeg(paste0(path,"Venn_",cell.type[k],"_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
do.call(cowplot::plot_grid, c(g_ec,ncol =2, scale = 0.8))
dev.off()

