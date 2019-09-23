library(Seurat)
library(ggpubr)
library(dplyr)
library(tidyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# load data
(load(file = "data/Lymphoma_EC_12_20190125.Rda"))
object %<>% SetAllIdent("res.0.6")
TSNEPlot(object,do.label = T)

# === T cells only=====
T_cells <- SubsetData(object, ident.remove = c(1,3,9,10,12))
T_cells@meta.data$orig.ident = gsub("TALL-plus-EC-2","TALL-plus-EC",T_cells@meta.data$orig.ident)
T_cells <- SetAllIdent(T_cells, id = "orig.ident")
TSNEPlot(T_cells, do.label = T)
T_cells %<>% SetAllIdent("conditions")
T_cells <- SubsetData(T_cells, ident.use = c("EC+T-ALL","T-ALL"))
genes <- c("IGF1R","HES4","DTX1")
T_cells@meta.data[,genes] = t(T_cells@data[genes,])
T_cells_meta.data = T_cells@meta.data[,c("orig.ident","conditions","tests",genes)]
T_cells_meta.data = T_cells_meta.data[!T_cells_meta.data$tests %in% c("test1","test5"),]
T_cells_meta.data$batch = T_cells_meta.data$tests
T_cells_meta.data$batch = gsub("test4","batch1",T_cells_meta.data$batch)
T_cells_meta.data$batch = gsub("test6","batch2",T_cells_meta.data$batch)
table(T_cells_meta.data$orig.ident,T_cells_meta.data$batch)
T_cells_meta.data$conditions = gsub("EC\\+T-ALL","T-ALL\\+EC",T_cells_meta.data$conditions)
# Create bar plots of means
for(gene in genes){
        jpeg(paste0(path,"barplot_T_cells_",gene,".jpeg"), units="in", width=10, height=7,res=600)
        g <- ggbarplot(T_cells_meta.data, x = "batch", y = gene, 
                  add = c("point", "jitter"),ylab = paste(gene, "(UMI)"),
                  color = "conditions", palette = c("#00AFBB", "#E7B800"),
                  position = position_dodge(0.8))
        print(g)
        dev.off()  
}

T_cells_meta.data_sum <- gather(T_cells_meta.data[,c("conditions",genes)],
                                key = "genes",value = "UMI", -conditions)
jpeg(paste0(path,"barplot_T_cells.jpeg"), units="in", width=10, height=7,res=600)
g <- ggbarplot(T_cells_meta.data_sum, x = "genes", y = "UMI", 
               add = c("point", "jitter"),
               color = "conditions", palette = c("#00AFBB", "#E7B800"),
               position = position_dodge(0.8))
print(g)
dev.off()

# === EC cells only=====
EC <- SubsetData(object, ident.use = c(1,3,9,12))
EC@meta.data$orig.ident = gsub("TALL-plus-EC-2","TALL-plus-EC",EC@meta.data$orig.ident)
EC %<>% SetAllIdent("conditions")
EC %<>% SubsetData(ident.use = c("EC","EC+T-ALL"))
table(EC@meta.data$orig.ident,EC@meta.data$conditions)
table(EC@meta.data$orig.ident,EC@meta.data$tests)
EC %<>% SetAllIdent("tests")
EC %<>% SubsetData(ident.use = c("test4","test6"))
genes <- FilterGenes(EC,c("STAT1","IFI27","ANXA2","ISG15","CAV1"))
EC@meta.data[,genes] = t(EC@data[genes,])
EC_meta.data = EC@meta.data[,c("orig.ident","conditions","tests",genes)]
EC_meta.data$batch = EC_meta.data$tests
EC_meta.data$batch = gsub("test4","batch1",EC_meta.data$batch)
EC_meta.data$batch = gsub("test6","batch2",EC_meta.data$batch)

for(gene in genes){
        jpeg(paste0(path,"barplot_EC_",gene,".jpeg"), units="in", width=10, height=7,res=600)
        g <- ggbarplot(EC_meta.data, x = "batch", y = gene, 
                       add = c("point", "jitter"),ylab = paste(gene, "(UMI)"),
                       color = "conditions", palette = c("#00AFBB", "#E7B800"),
                       position = position_dodge(0.8))
        print(g)
        dev.off()  
}

EC_meta.data_sum <- gather(EC_meta.data[,c("conditions",genes)],
                                key = "genes",value = "UMI", -conditions)
jpeg(paste0(path,"barplot_EC.jpeg"), units="in", width=10, height=7,res=600)
g <- ggbarplot(EC_meta.data_sum, x = "genes", y = "UMI", 
               add = c("point", "jitter"),
               color = "conditions", palette = c("#00AFBB", "#E7B800"),
               position = position_dodge(0.8))
print(g)
dev.off()
