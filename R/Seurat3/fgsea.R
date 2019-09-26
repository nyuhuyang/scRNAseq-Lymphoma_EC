########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(ggplot2)
library(fgsea)
library(tibble)
library(ggpubr)
library(ggsci)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 3.1.1 load data
dates <- c("2018-10-18","2018-12-30")
cluster <- c("EC+RO2", "EC+3119")
cell.type = c("EC","TALL")
k=1
res <- list()
for(i in seq_along(c("2018-10-18","2018-12-30"))){
        res[[i]] = read.csv(file=paste0("output/20190924/",cell.type[k],"_markers_",dates[i],".csv"),
                            row.names = 1,stringsAsFactors=F)
        #res[[i]] = res[[i]][res[[i]]$cluster %in% "TEC",]
        res[[i]]$gene = rownames(res[[i]])
        res[[i]]$cluster = cluster[i]
}

res = bind_rows(res)
res = res[order(res["p_val_adj"]),]
hist(res$avg_logFC,breaks = 200, xlim = c(-1,1))
head(res, 20)

hallmark <- gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
biocarta <- gmtPathways("../seurat_resources/msigdb/c2.cp.biocarta.v6.2.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v6.2.symbols.gmt")
tft <- gmtPathways("../seurat_resources/msigdb/c3.tft.v6.2.symbols.gmt")
c6 <- gmtPathways("../seurat_resources/msigdb/c6.all.v6.2.symbols.gmt")
GO <- gmtPathways("../seurat_resources/msigdb/c5.all.v6.2.symbols.gmt")
allpathways <- c(hallmark,biocarta,kegg)

hallmark %>% head() %>% lapply(head)
biocarta %>% head() %>% lapply(head)

names(hallmark) = gsub("HALLMARK","",names(hallmark))
names(hallmark) = gsub("_"," ",names(hallmark))
names(allpathways) = gsub("_"," ",names(allpathways))

hallmark_fgesa <- FgseaDotPlot(stats=res, pathways=hallmark, nperm=1000,
             avg_logFC =0.01 ,padj = 0.25,pval = 0.05,
             order.by = c("EC+RO2","NES"),decreasing = F,
             size = "-log10(pval)", fill = "NES",do.return = T,
             sample = paste("Educated",cell.type[k]), 
             pathway.name = "Hallmark",rotate.x.text = T,
             font.xtickslab=14, font.main=18,font.ytickslab = 10,
             font.legend = list(size = 14),font.label = list(size = 14),
             width=8, height=7,hjust=0.8)
write.csv(hallmark_fgesa,paste0(path,cell.type[k],"_hallmark_fgesa.csv"))

allpathways_fgesa <- FgseaDotPlot(stats=res, pathways=allpathways, nperm=1000,
             avg_logFC =0.01 ,padj = 0.25,pval = 0.05,
             order.by = c("EC+RO2","NES"),decreasing = F,
             size = "-log10(pval)", fill = "NES",do.return = T,
             sample = paste("Educated",cell.type[k]), 
             rotate.x.text = T, pathway.name = "Hallmark, biocarta,and KEGG",
             font.xtickslab=14, font.main=12,font.ytickslab = 8,
             font.legend = list(size = 12),font.label = list(size = 12),
             width=6, height=7,hjust=0.67)
write.csv(allpathways_fgesa,paste0(path,cell.type[k],"allpathways_fgesa.csv"))
