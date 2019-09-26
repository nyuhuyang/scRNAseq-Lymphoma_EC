########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(ggplot2)
library(ggrepel)
rm(list=ls());gc()
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# 3.1.1 load data ==================
dates <- c("2018-10-18","2018-12-30")
cluster <- c("EC+RO2", "EC+3119")
cell.type = c("EC","TALL")
k=1
dataset <- list()
for(i in seq_along(dates)){
        dataset[[i]] = read.csv(file=paste0("output/20190924/",cell.type[k],"_markers_",dates[i],".csv"),
                            row.names = 1,stringsAsFactors=F)
        dataset[[i]]$cluster = cluster[i]
        dataset[[i]]$gene = rownames(dataset[[i]])
        jpeg(paste0(path,"Volcano_plot_",cell.type[k],"_",cluster[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(VolcanoPlots(dataset[[i]],alpha = 0.9,size = 1.5,cut_off_logFC = 0.5)+
                      ggtitle(paste(cluster[i]," vs.",cell.type[k]))+
                      theme(plot.title = element_text(size=15, hjust = 0.5,face="plain")))
        dev.off()
}

