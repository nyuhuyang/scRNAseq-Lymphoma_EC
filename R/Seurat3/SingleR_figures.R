library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggthemes)z z zzzz
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file = "data/Lymphoma_EC_10_20190922.Rda"))
(load(file="output/singler_Lymphoma_EC_10T_20190922.Rda"))

# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object)
if(length(singler$singler[[1]]$SingleR.single$labels) < ncol(object)){
        all.cell = colnames(object);length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = subset(object, cells = know.cell)
}

table(rownames(singler$singler[[1]]$SingleR.single$labels) == colnames(object))
##############################
# add singleR label to Seurat
###############################

singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident" = gsub("\\_.*","",rownames(singler$singler[[1]]$SingleR.single$labels)),
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))
head(singlerDF)
singlerDF = singlerDF[colnames(object),]
table(rownames(singlerDF) == colnames(object))
apply(singlerDF,2,function(x) length(unique(x)))

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR::SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR::SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

table(singlerDF$singler1sub) %>% kable %>% kable_styling()
##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
apply(singlerDF[,c("singler1sub","singler1main")],2,function(x) length(unique(x)))
singlerDF[,c("singler1main")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaData(object = object,metadata = singlerDF)
object <- AddMetaColor(object = object, label= "singler1main", colors = singler_colors1)
Idents(object) <- "singler1main"
object %<>% subset(idents = "Erythrocytes", invert = TRUE)

UMAPPlot.1(object, cols = ExtractMetaColor(object),label = T, label.repel = T,pt.size = 1,
         label.size = 4, repel = T,no.legend = T,do.print = T,
         title = "Cell types identified by Blueprint + Encode")
TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T, label.repel = T,pt.size = 1,
           label.size = 4, repel = T,no.legend = T,do.print = T,
           title = "Cell types identified by Blueprint + Encode")
save(object,file="data/Lymphoma_EC_10_20190922.Rda")