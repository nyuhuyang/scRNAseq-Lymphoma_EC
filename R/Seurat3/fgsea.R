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
res <- list()
for(i in seq_along(c("2018-10-18","2018-12-30"))){
        res[[i]] = read.csv(file=paste0("output/20190924/EC_markers_",dates[i],".csv"),
                       stringsAsFactors=F)
        res[[i]] = res[[i]][res[[i]]$cluster %in% "TEC",]
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
             avg_logFC =0.001 ,padj = 0.25,pval = 0.05,
             order.by = c("EC+RO2","NES"),decreasing = F,
             size = "-log10(pval)", fill = "NES",do.return = T,
             sample = "Educated endothelial cells", 
             pathway.name = "Hallmark",rotate.x.text = T,
             font.xtickslab=14, font.main=18,font.ytickslab = 10,
             font.legend = list(size = 14),font.label = list(size = 14),
             width=7, height=7,hjust=0.8)
write.csv(hallmark_fgesa,paste0(path,"hallmark_fgesa.csv"))

allpathways_fgesa <- FgseaDotPlot(stats=res, pathways=allpathways, nperm=1000,
             avg_logFC =0.05 ,padj = 0.2,pval = 0.05,
             order.by = c("EC+RO2","NES"),decreasing = F,
             size = "-log10(pval)", fill = "NES",do.return = T,
             sample = "Educated endothelial cells", 
             rotate.x.text = T, pathway.name = "Hallmark, biocarta,and KEGG",
             font.xtickslab=14, font.main=12,font.ytickslab = 10,
             font.legend = list(size = 12),font.label = list(size = 12),
             width=6, height=7,hjust=0.67)
write.csv(allpathways_fgesa,paste0(path,"allpathways_fgesa.csv"))

df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
groups = c("Untreated","Pt-17","Pt-25")
for(i in 1:length(groups)){
        res_B = read.csv(file = paste0("output/20190622/B/B_MCL_DE/B_",groups[i],".csv"))
        res_T = read.csv(file = paste0("output/20190622/T/T_NK_DE/T_",groups[i],".csv"))
        
        (samples = df_samples$sample[df_samples$sample %in% unique(res_B$cluster)])
        #res_B$cluster %<>% factor(levels = samples)
        #res_T$cluster %<>% factor(levels = samples)
        
        FgseaDotPlot(stats=res_B, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
                     order.by = c(4,"NES"),decreasing = F,
                     size = "-log10(pval)", fill = "NES",
                     sample = paste(groups[i],"B_MCL clusters"), 
                     pathway.name = "Hallmark",rotate.x.text = F)
        (samples = df_samples$sample[df_samples$sample %in% unique(res_T$cluster)])
        FgseaDotPlot(stats=res_T, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
                     order.by = c(4,"NES"),decreasing = F,
                     size = "-log10(pval)", fill = "NES",
                     sample = paste(groups[i],"T_NK clusters"), 
                     pathway.name = "Hallmark",rotate.x.text = F)

}

#============ T cells  =====================
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",4)
sample_n = which(df_samples$tests %in% tests)
df <- as.data.frame(df_samples[sample_n,])
(samples <- c("Normal",unique(df$sample)))
cell.type <- c("T_cells:CD4+","T_cells:CD8+","NK_cells")

res_T_list <- list()
for(i in 2:length(samples)){
        res_T_list[[i-1]] = read.csv(file = paste0(path,samples[i],"_vs_Normal",".csv"))
}
res_T <- do.call("rbind.data.frame", res_T_list)
res_T = res_T[grep("CD8+",res_T$cluster1.vs.cluster2),]
res_T$cluster1.vs.cluster2 %<>% as.character %>% gsub('\\..*',"",.)

res_T$cluster1.vs.cluster2 %<>% factor(levels = samples[2:5])
        
for(i in 2:length(samples)) FgseaBarplot(pathways=hallmark, stats=res_T, nperm=1000,
                                          cluster = samples[i],no.legend = F,
                                          cut.off = "padj",cut.off.value = 0.25,
                                          sample="CD8+ T cells of",pathway.name = "Hallmark", hjust=0.5,
                                          width=10, height=7)

FgseaDotPlot(stats=res_T, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
                     order.by = c(4,"NES"),decreasing = F,
                     size = "-log10(pval)", fill = "NES",
                     sample = paste("Pt-17's CD8+ T cells"), 
                     pathway.name = "Hallmark",rotate.x.text = F)
