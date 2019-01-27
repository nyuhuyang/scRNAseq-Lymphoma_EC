library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file="./output/CancerCell_20181024.RData")
lnames
# T cell
T_Cell <- HumanGenes(CancerCell,c("CD2","CD3G","CD3D","CD4","CD8A","IL2RA","FOXP3",
                           "IL7R","SELL","IL2RG","GIMAP5"))

Treg <- HumanGenes(CancerCell,c("FOXP3","CD4","IL2RA","CTLA4","PDCD1","ENTPD1","CD38",
                         "ICOS","TNFSF9","TNFRSF9"))
CD4_Naive_T <- HumanGenes(CancerCell,c("CD4","IL7R","GIMAP5","SELL","IL2RG"))
Endothelium <- HumanGenes(CancerCell,c("Cdh5","Pecam1","Flt1","Plvap","Kdr","ptprb",
                                   "Vwf","EMCN","Car4","VEGFA"))
EC <- HumanGenes(CancerCell,c("IFI27","CD3D","HLA-B","IFI6"))

for(i in 1:length(markers)) {
        
    jpeg(paste0("output/",EC[i],".jpeg"), units="in", width=10, height=7,
        res=600)
    p1 <- SingleFeaturePlot.1(object = CancerCell, feature = EC[i])
    print(p1)
    print(paste0(i,":",length(EC)))
    dev.off()
}
CancerCell1 <- SplitSeurat(CancerCell, split.by = "orig.ident")


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

gene_sets <- lapply(gene_sets, function(x) HumanGenes(CancerCell,x))

for(i in 1:length(gene_sets)){
    CancerCell <- .AddModuleScore(CancerCell, genes.list = gene_sets[i],
                      ctrl.size = 5,enrich.name = gene_sets_names[i])
}

# 2.3 Split into 6 sample groups ==============
table(CancerCell@meta.data$singler2sub)
EC_index <- which(CancerCell@meta.data$singler2sub %in% c("Endothelial cells",
                                                        "mv Endothelial cells"))
EC_plus_TALL <- intersect(EC_index,
                          which(CancerCell@meta.data$old.ident %in% "TALL-plus-EC"))
EC2_plus_TALL <- intersect(EC_index,
                          which(CancerCell@meta.data$old.ident %in% "TALL-plus-EC-2"))
CancerCell@meta.data[EC_plus_TALL,"old.ident"] = "EC-plus-TALL"
CancerCell@meta.data[EC2_plus_TALL,"old.ident"] = "EC2-plus-TALL"
table(CancerCell@meta.data$old.ident)
save(CancerCell, file="./output/CancerCell_20181024.RData")




