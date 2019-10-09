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
# 3.1.1 load pathway
set.seed(100)
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

# 3.1.2 read DE files and generate GSEA

# test normalization methods
assay = c("RNA","RNA","SCT","SCT")
slot = c("data","data_MNN","data","scale.data")
dates <- c("2018-10-18","2018-12-30")
cluster <- c("EC+RO2", "EC+3119")
cell.type = c("EC","TALL")
# test EC only
assay = "SCT"
slot = "data"
cell.type = "EC"
cluster <- c("RO2-plus-EC-2","RO2-plus-EC","EC2","EC+3119")
(dates <- c(paste("Endothelial cells",cluster[1:3],"vs_EC",sep = "_"),"2018-12-30"))

for(d in seq_along(assay)){
        assay_slot = paste0(path,assay[d],"_",slot[d],"/")
        if(!dir.exists(assay_slot)) dir.create(assay_slot, recursive = T)
        for(k in seq_along(cell.type)){
                dataset <- list()
                for(i in seq_along(dates)){
                        
                        dataset[[i]] = read.csv(file = paste0("Yang/6. Differential analysis/",
                                                              basename(assay_slot),"/", cell.type[k],"_markers_",
                                                              dates[i],"_",basename(assay_slot),".csv"),
                                                row.names = 1,stringsAsFactors=F)
                        if(basename(assay_slot) == "SCT_scale.data") {
                                dataset[[i]]$avg_logFC = tanh(dataset[[i]]$avg_logFC)
                                dataset[[i]]$avg_logFC = dataset[[i]]$avg_logFC*10
                        }
                        if(basename(assay_slot) == "RNA_data_MNN"){
                                dataset[[i]]$avg_logFC = dataset[[i]]$avg_logFC*200
                        }
                        dataset[[i]]$cluster = cluster[i]
                        dataset[[i]]$gene = rownames(dataset[[i]])
                }
                res = bind_rows(dataset)
                res = res[order(res["p_val_adj"]),]
                avg_logFC = 0; padj = 1; pval = 0.05
                hallmark_fgesa <- FgseaDotPlot(stats=res, pathways=hallmark, nperm=1000,
                                               order.by = c(cluster[1],"NES"),decreasing = F,
                                               size = "-log10(pval)", fill = "NES",do.return = T,
                                               sample = paste0("Educated_",cell.type[k]), 
                                               pathway.name = "Hallmark",rotate.x.text = T,
                                               font.xtickslab=14, font.main=18,font.ytickslab = 10,
                                               font.legend = list(size = 14),font.label = list(size = 14),
                                               width=6, height=7,hjust=0.8,
                                               save_path = paste(paste0(assay_slot,"Dotplot"),cell.type[k],"hallmark",
                                                                  avg_logFC,padj,pval,
                                                                 basename(assay_slot),".jpeg",sep = "_"))
                write.csv(hallmark_fgesa,paste0(assay_slot,cell.type[k],"_hallmark_fgesa",basename(assay_slot),".csv"))
                
                tryCatch(
                        {FgseaDotPlot(stats=res, pathways=allpathways, nperm=1000,
                                                  order.by = c(cluster[1],"NES"),decreasing = F,
                                                  size = "-log10(pval)", fill = "NES",do.return = T,
                                                  sample = paste("Educated",cell.type[k]), 
                                                  rotate.x.text = T, pathway.name = "Hallmark, biocarta,and KEGG",
                                                  font.xtickslab=14, font.main=12,font.ytickslab = 8,
                                                  font.legend = list(size = 12),font.label = list(size = 12),
                                                  width=6, height=7,hjust=0.67,
                                                  save_path = paste(paste0(assay_slot,"Dotplot"),cell.type[k],"allpathways",
                                                                    avg_logFC,padj,pval,
                                                                    basename(assay_slot),".jpeg",sep = "_"))
                        },
                        error = function(e) e
                )
        }
        svMisc::progress(d/length(assay)*100)
}
