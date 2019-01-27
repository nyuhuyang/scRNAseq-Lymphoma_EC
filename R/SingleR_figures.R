library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(wordcloud)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file = "output/20181231/singler_Lymphoma_EC_12_20190125.Rda"))
(load(file = "data/Lymphoma_EC_12_20190125.Rda"))
##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler2sub"=singler$singler[[1]]$SingleR.single$labels,
                       "singler2main"=singler$singler[[1]]$SingleR.single.main$labels,
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% object@cell.names)

#knowDF = data.frame("cell.names"= object@cell.names)
#ident.DF = full_join(singlerDF,knowDF, by="cell.names")
#ident.DF<- apply(ident.DF,2,as.character)
#rownames(ident.DF) = ident.DF[,"cell.names"]
#ident.DF = ident.DF[,-which(colnames(ident.DF) == "cell.names")]
apply(singlerDF,2,function(x) length(unique(x)))
#ident.DF[is.na(ident.DF)] <- "unknown"
object <- AddMetaData(object = object,
                   metadata = singlerDF)
object <- SetAllIdent(object = object, id = "singler2sub")

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"/DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"DrawHeatmap_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n = 50,normalize = F))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

# cell number
table(object@meta.data$singler2sub) %>% kable %>% kable_styling()
table(object@meta.data$singler2sub,object@meta.data$orig.ident) %>% kable %>% kable_styling

# cell percentage
prop.table(x = table(object@meta.data$singler2sub,
                     object@meta.data$orig.ident),margin = 2) %>%
  kable()  %>% kable_styling

singler$meta.data$orig.ident %>% table() %>% kable() %>% kable_styling()
object@meta.data$singler2main %>% table() %>% kable() %>% kable_styling()

#--------------- wordcloud ---------------------
singler2main <- object@meta.data$singler2main %>% table %>% as.data.frame
colnames(singler2main) = c("main_types","Freq")

set.seed(1234)
wordcloud(words = singler2main$main_types, freq = singler2main$Freq,scale=c(5,.5),
          min.freq = 1, max.words=20000, random.order=TRUE, rot.per=0.35, 
          colors= singler.colors[1:20])



##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors2[duplicated(singler_colors2)]
length(singler_colors2)
apply(object@meta.data[,c("singler1sub","singler1main","singler2sub","singler2main")],
      2,function(x) length(unique(x)))
object@meta.data[,c("singler2sub")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaColor(object = object, label= "singler2sub", colors = singler_colors2)
object <- SetAllIdent(object = object, id = "singler2sub")
TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F)

##############################
# draw tsne plot
##############################
p4 <- DimPlot.1(object = object, do.label = T, group.by = "ident", 
                reduction.use= "tsne",
                 do.return = TRUE, no.legend = T, label.repel =T,
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 3,force = 2)+
  ggtitle("Cell type labeling by Blueprint + Encod")+
  theme(text = element_text(size=10),							
        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"Plottsne_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
print(p4)
dev.off()

save(object,file="./data/Lymphoma_EC_12_20190125.Rda")
##############################
# subset Seurat
###############################
table(object@meta.data$orig.ident)
table(object@ident)

df_samples <- readxl::read_excel("doc/181230_Single_cell_TALL_sample_list.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower

DF <- object@meta.data
DF$tests = DF$orig.ident
DF$tests = gsub("3119PR","test1",DF$tests)
DF$tests = gsub("EC-only|TALL-only|TALL-plus-EC|TALL-plus-EC-2","test4",DF$tests)
DF$tests = gsub("A-3119|E-3119|S-3119","test5",DF$tests)
DF$tests = gsub("EC|3119-plus-EC|3119","test6",DF$tests)

object@meta.data = DF

tests <- paste0("test",4:6)
for(test in tests){
        cell.use <- rownames(object@meta.data)[object@meta.data$tests %in% test]
        subset.object <- SubsetData(object, cells.use = cell.use)
        
        g <- SplitDimPlot(subset.object,group.by = "ident",split.by = "orig.ident",
                          reduction.use= "tsne",
                           no.legend = T,do.label =F,label.size=3,
                           return.plots =T, label.repel = T,force=2)
        jpeg(paste0(path,test,"_conditions_tsnePlot.jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        dev.off()
}

save(object, file = "./data/Lymphoma_EC_12_20190125.Rda")
