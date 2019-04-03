library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(magrittr)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 2.1 identify phenotype for each cluster  ==========================================
(load(file="./output/object_20181024.RData"))
# T cell
T_Cell <- FilterGenes(object,c("CD3G","CD3D","CD3E","CD8A","MALAT1","CD1B",
                               "HLA-B","TUBB","PCNA","TGFB1","ZNF90","CXCR4"))


EC <- FilterGenes(object,c("FN1","IGFBP4","TGM2","EMCN","SRGN",
                           "VEGFA","VWF","KDR","FLT1","CDH5","FTL","VIM",
                           "CCND1"))

markers <-FilterGenes(object,c(T_Cell,EC))
markers <- c("FTL","HLA-B","MALAT1","TUBB")
for(i in 1:length(markers)) {
    jpeg(paste0(path,markers[i],".jpeg"), units="in", width=10, height=7,
        res=600)
    p1 <- SingleFeaturePlot.1(object = object, feature = markers[i],
                              gradient.use = c("lightgrey","blue4"),
                              threshold= NULL)
    print(p1)
    print(paste0(i,":",length(markers)))
    dev.off()
}
object1 <- SplitSeurat(object, split.by = "orig.ident")


# 2.2 add GeneSets score ==============
# Read in a list of Angiocrine_factors Nature 2016
Angiocrine_factors <- read.delim("data/seurat_resources/Angiocrine_factors.txt",
                                 stringsAsFactors =F)
Angiocrine_factors <- Angiocrine_factors[,1]
# Read other genesets
(gene_setsfiles <- list.files("data/seurat_resources/"))
(gene_setsfiles <- gene_setsfiles[c(2,3,4,5)])
gene_sets <- lapply(paste0("data/seurat_resources/",gene_setsfiles),
                    function(x) read.delim(x, stringsAsFactors =F)) %>%
    lapply(function(x) x[-1,1]) # extract first column, remove first row
gene_sets[[(length(gene_sets)+1)]] = Angiocrine_factors
gene_sets_names <- c(sub(".txt","",gene_setsfiles),"Angiocrine_factors")

gene_sets <- lapply(gene_sets, function(x) FilterGenes(object,x))

for(i in 1:length(gene_sets)){
    object <- .AddModuleScore(object, genes.list = gene_sets[i],
                      ctrl.size = 5,enrich.name = gene_sets_names[i])
}

# 2.3 Split into 6 sample groups ==============
table(object@meta.data$singler2sub)
EC_index <- which(object@meta.data$singler2sub %in% c("Endothelial cells",
                                                        "mv Endothelial cells"))
EC_plus_TALL <- intersect(EC_index,
                          which(object@meta.data$old.ident %in% "TALL-plus-EC"))
EC2_plus_TALL <- intersect(EC_index,
                          which(object@meta.data$old.ident %in% "TALL-plus-EC-2"))
object@meta.data[EC_plus_TALL,"old.ident"] = "EC-plus-TALL"
object@meta.data[EC2_plus_TALL,"old.ident"] = "EC2-plus-TALL"
table(object@meta.data$old.ident)
save(object, file="./output/object_20181024.RData")


# 2.4 Find conservative markers ==============
head(object@meta.data)
object %<>% SetAllIdent("res.0.6")
TSNEPlot(object)
res_markers <- FindAllMarkers.UMI(object,logfc.threshold=1,min.pct = 0.8)
write.csv(res_markers,paste0(path,"res_markers.csv"))
