library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
marker_path <- paste0(path,"markers/")
if(!dir.exists(marker_path))dir.create(marker_path, recursive = T)

#====== 2.1 Identify cell types ==========================================
(load(file="data/Lymphoma_EC_10_20190922.Rda"))
DefaultAssay(object) <- 'RNA'
df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.Main")

colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) x[1:18]) %>% 
    lapply(function(x) FilterGenes(object,x)) %>% 
    lapply(function(x) x[!is.na(x)]) %>% 
    lapply(function(x) x[1:min(length(x),12)])
marker.list <- marker.list[!is.na(marker.list)]
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()
Idents(object) <- "integrated_snn_res.0.6"

for(i in 1:length(marker.list)){
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, features = marker,pt.size = 0.5, label=T,
                    reduction = "tsne")+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(marker_path,names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    svMisc::progress(i/length(marker.list)*100)
}

#======== rename ident =================
object %<>% RenameIdents("0" = "T cells",
                         "1" = "Endothelial cells",
                         "2" = "T cells",
                         "3" = "T cells",
                         "4" = "Endothelial cells",
                         "5" = "T cells",
                         "6" = "T cells",
                         "7" = "T cells",
                         "8" = "Endothelial cells",
                         "9" = "T cells",
                         "10" = "T cells")
object@meta.data$cell.type = as.character(Idents(object))
HSC <- grepl("HSC",object@meta.data$singler1main)
object@meta.data[HSC,"cell.type"] ="HSC"
Idents(object) = "cell.type"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell.type", 
                       colors = c("#66A61E","#FF7F00","#386CB0"))
TSNEPlot.1(object, group.by = "cell.type",cols = ExtractMetaColor(object),label = T,
           label.repel = T, pt.size = 0.5,label.size = 4, repel = T,no.legend = T,
           do.print = T,do.return = T,title = "Cell types")
save(object,file="data/Lymphoma_EC_10_20190922.Rda")

# subset
Idents(object) = "date"
(dates <- unique(Idents(object)))
for(date in dates){
    subset_object <- subset(object, idents = date)
    Idents(subset_object) = "cell.type"
    TSNEPlot.1(subset_object, group.by = "cell.type",split.by = "orig.ident",
               cols = ExtractMetaColor(subset_object),label = T,
               label.repel = T, pt.size = 0.5,label.size = 4, repel = T,no.legend = T,
               do.print = T,do.return = T,border = T,
               title = paste("tSNE plots of all samples on",unique(subset_object$date)),
               unique.name = "date")
}
remove(subset_object);GC()
