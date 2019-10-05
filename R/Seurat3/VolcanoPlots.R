########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(ggplot2)
library(ggrepel)
rm(list=ls());gc()
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# set group
assay = c("RNA","RNA","SCT","SCT")
slot = c("data","data_MNN","data","scale.data")
dates <- c("2018-10-18","2018-12-30")
cluster <- c("EC+RO2", "EC+3119")
cell.type = c("EC","TALL")
# read data and generate VolcanoPlots
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
                                dataset[[i]]$avg_logFC = dataset[[i]]$avg_logFC/log2(exp(1))
                                dataset[[i]]$avg_logFC = tanh(dataset[[i]]$avg_logFC)
                        }
                        if(basename(assay_slot) == "RNA_data_MNN"){
                                dataset[[i]]$avg_logFC = dataset[[i]]$avg_logFC*100
                        }
                        dataset[[i]]$cluster = cluster[i]
                        dataset[[i]]$gene = rownames(dataset[[i]])
                        jpeg(paste0(assay_slot,"Volcano_plot_",cell.type[k],"_",cluster[i],"_",
                                    basename(assay_slot),".jpeg"), units="in", width=10, height=7,res=600)
                        print(VolcanoPlots(dataset[[i]],alpha = 0.9,size = 1.5,cut_off_logFC = 0.25,
                                           cut_off_pvalue = 10^(-5))+
                                      ggtitle(paste(cluster[i]," vs.",cell.type[k]))+
                                      theme(plot.title = element_text(size=15, hjust = 0.5,face="plain")))
                        dev.off()
                }
        }
        svMisc::progress(d/length(assay)*100)
}

   

 