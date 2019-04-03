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
library(magrittr)
library(sva)
library(kableExtra)
library(readr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 5.1 load data ==============
(load(file = "data/Lymphoma_EC_12_20190125.Rda"))

# GSEA for EC ==================================
object <- SetAllIdent(object, id = "tests")
CancerCell <- SubsetData(object, ident.use = c("test4"))
CancerCell <- SetAllIdent(CancerCell, id = "res.0.6")
table(CancerCell@ident)
TSNEPlot(CancerCell, do.label = T)
EC <- SubsetData(CancerCell, ident.use = c(1,3,9,12))
#EC@meta.data$orig.ident = gsub("TALL-plus-EC-2","TALL-plus-EC",
#                               EC@meta.data$orig.ident)
EC <- SetAllIdent(EC, id = "orig.ident")
table(EC@ident)
EC <- SubsetData(EC, ident.use = c("EC-only","TALL-plus-EC"))
EC %<>% NormalizeData
PrepareGSEA(EC, k = 100)

# Run GSEA
GSEA_EC <- read_delim("output/20190131/gsea_report_for_EC_only_vs_TALL_plus_EC_1548993830806.xls",
                        "\t", escape_double = FALSE, trim_ws = TRUE)
GSEA_EC %>% head(10) %>% kable() %>% kable_styling()
(gsea_path <- paste0("~/gsea_home/output/",tolower(format(Sys.Date(), "%b%d")), 
                     "/h.all.EC_vs_3119.plus.EC.Gsea.1548978407433"))
(h.all.EC_vs_3119.plus.EC <- sapply(GSEA_EC$NAME[1:9], function(name) {
                                paste0("enplot_",name, "_([0-9]+)*\\.png$")}) %>%
                sapply(function(x) list.files(path = gsea_path, pattern =x)) %>%
                .[sapply(.,length)>0] %>% #Remove empty elements from list with character(0)
                paste(gsea_path, ., sep = "/")) 
CombPngs(h.all.EC_vs_3119.plus.EC, ncol = 3)


# GSEA for T_cells ==================================
object <- SetAllIdent(object, id = "tests")
CancerCell <- SubsetData(object, ident.use = c("test5","test6"))
CancerCell <- SetAllIdent(CancerCell, id = "res.0.6")
TSNEPlot(CancerCell, do.label = T)
table(CancerCell@ident)
T_cells <- SubsetData(CancerCell, ident.remove = c(1,3,9,10,12))
#T_cells@meta.data$orig.ident = gsub("TALL-plus-EC-2","TALL-plus-EC",T_cells@meta.data$orig.ident)
T_cells <- SetAllIdent(T_cells, id = "orig.ident")
T_cells <- SubsetData(T_cells, ident.use = c("3119","E-3119"))
TSNEPlot(T_cells, do.label = T)
T_cells %<>% NormalizeData
PrepareGSEA(T_cells, k = 100)

# Run GSEA
GSEA_T <- read_delim("output/20190201/gsea_report_for_3119_versus_E-3119.xls",
                      "\t", escape_double = FALSE, trim_ws = TRUE)
GSEA_T %>% head(16) %>% kable() %>% kable_styling()
(gsea_path <- paste0("~/gsea_home/output/",tolower(format(Sys.Date(), "%b%d")), 
                     "/c2.cp.biocarta.T.EC_vs_T.Gsea.1548456752605"))
(c2.cp.biocarta.T.EC_vs_T <- sapply(GSEA_T$NAME[1:9], function(name) {
        paste0("enplot_",name, "_.*\\.png$")}) %>%
                sapply(function(x) list.files(path = gsea_path, pattern =x)) %>%
                .[sapply(.,length)>0] %>% #Remove empty elements from list with character(0)
                paste(gsea_path, ., sep = "/")) 
CombPngs(c2.cp.biocarta.T.EC_vs_T, ncol = 3)

#====== 2.1 pathway analysis ==========================================

gene_set = read.delim("./doc/h.all.v6.2.symbols.gmt",row.names =1,header = F,
                      stringsAsFactors = F)
gene_set <- gene_set[,-1]
gene_set %>% kable() %>% kable_styling()
gene_set.df <- as.data.frame(t(gene_set))

gene_set.list <- df2list(gene_set.df)
gene_set_list <- lapply(gene_set.list, function(x) HumanGenes(T_cells,x))

for(i in 1:length(gene_set_list)){
        T_cells <- .AddModuleScore(T_cells, genes.list = gene_set_list[i],
                                         ctrl.size = 5,enrich.name = names(gene_set_list[i]))
}

SplitSingleFeaturePlot(T_cells, markers = names(gene_set_list),
                       group.by = "ident",split.by = "orig.ident",
                       no.legend = T,label.size=3,do.print =T,
                       threshold = 0.1)

table(T_cells@meta.data$HALLMARK_NOTCH_SIGNALING >0.05,
      T_cells@meta.data$orig.ident) %>% prop.table(margin = 2)

cell.use = T_cells@meta.data$orig.ident %in% "TALL-only" %>%
        rownames(T_cells@meta.data)[.]
length(cell.use)
remove = sample(cell.use, size = (1415-895))
cells.use <- T_cells@cell.names[!(T_cells@cell.names %in% remove)]
new_T_cells <- SubsetData(T_cells, cells.use = cells.use)

SplitSingleFeaturePlot(new_T_cells, markers = "HALLMARK_NOTCH_SIGNALING",
                       group.by = "ident",split.by = "orig.ident",
                       no.legend = T,label.size=3,do.print =T,
                       threshold = 0.05)
