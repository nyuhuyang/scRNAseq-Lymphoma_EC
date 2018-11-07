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

test.markers <- unique(c(T_Cell,Treg,CD4_Naive_T,Endothelium))
for(i in 1:length(test.markers)) {
    jpeg(paste0(path,"/markers/",test.markers[i],".jpeg"), units="in", width=10, height=7,
        res=600)
    p1 <- SingleFeaturePlot.1(object = CancerCell, feature = test.markers[i])
    print(p1)
    print(paste0(i,":",length(test.markers)))
    dev.off()
}
