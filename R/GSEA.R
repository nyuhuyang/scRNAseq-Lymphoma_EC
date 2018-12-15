########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(MAST)
library(sva)
library(kableExtra)
library(readr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 5.1 load data ==============
(lnames = load(file="./output/CancerCell_20181024.RData"))

# Differential analysis for EC ==================================
table(CancerCell@ident)
EC <- SubsetData(CancerCell, ident.use = c("Endothelial cells","mv Endothelial cells"))
EC@meta.data$orig.ident <- gsub("TALL-plus-EC-2","TALL-plus-EC",EC@meta.data$orig.ident)
EC@meta.data$orig.ident <- gsub("TALL-plus-EC","EC-plus-TALL",EC@meta.data$orig.ident)
EC <- SetAllIdent(EC, id = "orig.ident")

PrepareGSEA(EC, k = 15)

# Run GSEA
GSEA_EC <- read_delim("output/20181215/gsea_report_for_EC-plus-TALL_1544851750346.xls",
                        "\t", escape_double = FALSE, trim_ws = TRUE)
GSEA_EC %>% head(20) %>% kable() %>% kable_styling()
(gsea_path <- paste0("~/gsea_home/output/",tolower(format(Sys.Date(), "%b%d")), 
                     "/my_analysis.Gsea.","1544851750346"))
(h.all <- sapply(GSEA_EC$NAME[1:9], function(name) {
                                paste0("enplot_",name, "_.*\\.png")}) %>%
                sapply(function(x) list.files(path = gsea_path, pattern =x)) %>%
                paste(gsea_path, ., sep = "/"))
CombPngs(h.all, ncol = 3)


# T_cells GSEA
T_cells <- SubsetData(CancerCell, ident.remove = c("Endothelial cells","mv Endothelial cells"))
remove <- FeaturePlot(T_cells,features.plot = "CD3D", do.identify = T)
cells.use <- T_cells@cell.names[!(T_cells@cell.names %in% remove)]
T_cells <- SubsetData(CancerCell, cells.use = cells.use)

T_cells <- SetAllIdent(T_cells, id = "orig.ident")
T_cells@meta.data$orig.ident <- gsub("TALL-plus-EC-2","TALL-plus-EC",
                                     T_cells@meta.data$orig.ident)
T_cells <- SetAllIdent(T_cells, id = "orig.ident")
GSEA_T_expr <- PrepareGSEA(T_cells)
GSEA_T_expr %>% head(20) %>% kable() %>% kable_styling()

GSEA_T_c2 <- read_delim("output/20181115/gsea_report_for_TALL-plus-EC_1542317084239.xls",
                          "\t", escape_double = FALSE, trim_ws = TRUE)

GSEA_T_c2 %>% head(20) %>% kable() %>% kable_styling()
