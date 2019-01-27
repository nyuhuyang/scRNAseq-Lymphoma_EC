########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(SingleR)
library(kableExtra)
library(eulerr)
library(magrittr)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="./data/Lymphoma_EC_12_20190125.Rda"))

# 4.1 venn diagram for T cells ==================================
object <- SetAllIdent(object, id = "tests")
CancerCell <- SubsetData(object, ident.use = c("test6"))
CancerCell <- SetAllIdent(CancerCell, id = "res.0.6")
table(CancerCell@ident)
T_cells <- SubsetData(CancerCell, ident.remove = c(1,3,9,10,12))
T_cells <- SetAllIdent(T_cells, id="orig.ident")
T_cells <- SubsetData(T_cells, ident.remove = "EC")
table(T_cells@ident)
T_cells %<>% NormalizeData
T_cells_markers <- FindAllMarkers.UMI(T_cells, return.thresh = 1,
                                      test.use = "MAST",min.pct = -Inf,
                                      min.cells.gene = -Inf, min.cells.group = -Inf,
                                      only.pos = F, # don't change it!
                                      logfc.threshold = -Inf)
write.csv(T_cells_markers,paste0(path,"T_cells_markers1.csv"))
T_cells_markers = read.csv(paste0(path,"T_cells_markers1.csv"), header = 1, row.names = 1)
table(T_cells_markers$cluster)


##=====================================
# Generage venn_diagrams

T_cells_markers$cluster = gsub("3119-plus-EC$","T cells after co-culture with EC",T_cells_markers$cluster)
T_cells_markers$cluster = gsub("3119$","T cells",T_cells_markers$cluster)

g_t <- mapply(function(x,y) eulerr(T_cells_markers, cut_off = x, cut_off_value = y,
                                 do.print = FALSE,do.lenged = FALSE), 
            c("p_val","p_val_adj","avg_logFC"),
            c(0.05,0.05,0.1),SIMPLIFY = FALSE)

jpeg(paste0(path,"Venn_T_all_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
do.call(cowplot::plot_grid, c(g_t,ncol =2, scale = 0.8))
dev.off()

# 4.2 venn diagram for EC cells ################################
object <- SetAllIdent(object, id = "tests")
CancerCell <- SubsetData(object, ident.use = c("test6"))
CancerCell <- SetAllIdent(CancerCell, id = "res.0.6")
table(CancerCell@ident)
EC <- SubsetData(CancerCell, ident.use = c(1,3,9,12))
EC <- SetAllIdent(EC, id = "orig.ident")
table(EC@ident)
EC <- SubsetData(EC, ident.remove = "3119")
EC %<>% NormalizeData
EC_markers <- FindAllMarkers.UMI(EC,return.thresh = 1,
                                      test.use = "MAST",min.pct = -Inf,
                                      min.cells.gene = -Inf, min.cells.group = -Inf,
                                      only.pos = F, # don't change
                                 logfc.threshold = -Inf)
write.csv(EC_markers,paste0(path,"EC_markers1.csv"))

##=====================================
# Generage venn_diagrams
EC_markers$cluster = gsub("3119-plus-EC$","E cells after co-culture with T-ALL",EC_markers$cluster)
EC_markers$cluster = gsub("EC$","E cells",EC_markers$cluster)

g_e <- mapply(function(x,y) eulerr(EC_markers, cut_off = x, cut_off_value = y,
                                   do.print = FALSE,do.lenged = F), 
              c("p_val","p_val_adj","avg_logFC"),
              c(0.05,0.05,0.1),SIMPLIFY = FALSE)

jpeg(paste0(path,"Venn_EC_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
do.call(cowplot::plot_grid, c(g_e,ncol =2, scale = 0.8))
dev.off()

# 4.3 venn diagram for EC cells between 2 bateches ################################
object <- SetAllIdent(object, id = "tests")
CancerCell <- SubsetData(object, ident.use = c("test4","test6"))
CancerCell <- SetAllIdent(CancerCell, id = "res.0.6")
table(CancerCell@ident)
EC2 <- SubsetData(CancerCell, ident.use = c(1,3,9,12))

EC2 <- SetAllIdent(EC2, id = "orig.ident")
table(EC2@ident)
EC2 <- SubsetData(EC2, ident.use = c("EC","EC-only"))
EC2 %<>% NormalizeData
# Combat
m = as.matrix(EC2@data)
m = m[rowSums(m)>0,]
tests = EC2@meta.data$tests
com = sva::ComBat(m, tests, prior.plots=FALSE, par.prior=TRUE)
EC2@data = Matrix(com)
EC_only <- FindAllMarkers.UMI(EC2,return.thresh = 1,
                                 test.use = "MAST",min.pct = -Inf,
                                 min.cells.gene = -Inf, min.cells.group = -Inf,
                                 only.pos = F, # don't change
                                 logfc.threshold = -Inf)
write.csv(EC_only,paste0(path,"EC_only_markers.csv"))

##=====================================
# Generage venn_diagrams
EC_only$cluster = gsub("EC-only$","EC in October",EC_only$cluster)
EC_only$cluster = gsub("EC$","EC in December",EC_only$cluster)

g_ec <- mapply(function(x,y) eulerr(EC_only, cut_off = x, cut_off_value = y,
                                   do.print = FALSE,do.lenged = F), 
              c("p_val","p_val_adj","avg_logFC"),
              c(0.05,0.05,0.1),SIMPLIFY = FALSE)

jpeg(paste0(path,"Venn_EC-only_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
do.call(cowplot::plot_grid, c(g_ec,ncol =2, scale = 0.8))
dev.off()

jpeg(paste0(path,"Venn_EC-only_pval_2batech.jpeg"), units="in", width=10, height=7,res=600)
do.call(cowplot::plot_grid, c(list(g_e$p_val,g_ec$p_val),ncol =2, scale = 0.7))
dev.off()
# 4.4 venn diagram for EC cells between 2 bateches before and after co-culture #####################
object <- SetAllIdent(object, id = "tests")
CancerCell <- SubsetData(object, ident.use = c("test4","test6"))
CancerCell <- SetAllIdent(CancerCell, id = "res.0.6")
table(CancerCell@ident)
EC4 <- SubsetData(CancerCell, ident.use = c(1,3,9,12))
table(EC4@meta.data$orig.ident)
EC4@meta.data$orig.ident = gsub("TALL-plus-EC-2$","TALL-plus-EC",EC4@meta.data$orig.ident)
EC4 %<>% SetAllIdent(id = "orig.ident")
table(EC4@ident)
EC4 %<>% SubsetData(ident.remove =  c("3119","TALL-only"))
EC4 %<>% NormalizeData
# Combat
m = as.matrix(EC4@data)
m = m[rowSums(m)>0,]
tests = EC4@meta.data$tests
com = sva::ComBat(m, tests, prior.plots=FALSE, par.prior=TRUE)
EC4@data = Matrix(com)

EC_4 <- FindAllMarkers.UMI(EC4,return.thresh = 1,
                              test.use = "MAST",min.pct = -Inf,
                              min.cells.gene = -Inf, min.cells.group = -Inf,
                              only.pos = F, # don't change
                              logfc.threshold = -Inf)
write.csv(EC_4,paste0(path,"EC_4_comparision.csv"))

# Generage venn_diagrams

EC_4$cluster = gsub("TALL-plus-EC$","co-cultured EC in October",EC_4$cluster)
EC_4$cluster = gsub("3119-plus-EC$","co-cultured EC in December",EC_4$cluster)
EC_4$cluster = gsub("EC-only$","EC in October",EC_4$cluster)
EC_4$cluster = gsub("EC$","EC in December",EC_4$cluster)
table(EC_4$cluster)

g_ec4 <- mapply(function(x,y) eulerr(EC_4, cut_off = x, cut_off_value = y,
                                    do.print = FALSE,do.lenged = T), 
               c("p_val","p_val_adj","avg_logFC"),
               c(0.05,0.05,0.1),SIMPLIFY = FALSE)

eulerr(EC_4, cut_off = "p_val", cut_off_value = 0.05, do.print = T,do.lenged = T,do.return = F)

jpeg(paste0(path,"Venn_EC4_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
do.call(cowplot::plot_grid, c(g_ec4,ncol =2))
dev.off()


# 20181204 ======
##-----------------------------------
jpeg(paste0(path,"/Venn_T_TE_pval.jpeg"), units="in", width=10, height=7,res=600)
p_val = 0.05#c(0.5,0.1,0.05,0.01)
g <- list()
for(i in 1:length(p_val)){
        fit <- eulerr(T_cells_markers,key = c("3119", "3119-plus-EC"),
                        cut_off = "p_val", cut_off_value = p_val[i])
        g[[i]] <- plot(fit, quantities = TRUE, fill = c("#377EB8","#A6D854"),lty = 1:2,
                       labels = list(labels = c("TALL", "TALL\n+EC")), legend = F,
                       main = paste("p_val:", p_val[i]))
}
do.call(cowplot::plot_grid, c(g,ncol =2))
dev.off()

##-----------------------------------
jpeg(paste0(path,"/Venn_T_TE_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
cut_off_values <- c(c(0.1,0.05,0.01),c(0.99,0.5,0.1),c(0.1,0.05,0.01))
names(cut_off_values) <- c(rep("p_val",3),rep("p_val_adj",3),rep("avg_logFC",3))
g <- list()
for(i in 1:length(cut_off_values)){
        fit <- eulerr(T_cells_markers, key = c("3119", "TALL-plus-EC"),
                        cut_off = names(cut_off_values[i]),
                        cut_off_value = cut_off_values[i])
        g[[i]] <- plot(fit, quantities = TRUE, fill = c("#377EB8","#A6D854"),lty = 1:2,
                       main = paste(names(cut_off_values[i]),":",cut_off_values[i]))
}
do.call(cowplot::plot_grid, c(g,ncol =3))
dev.off()

##-----------------------------------
jpeg(paste0(path,"/Venn_T_3_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
colors_list <- list(c("#377EB8","#A6D854"),c("#377EB8","#E41A1C"),
                    c("#A6D854","#E41A1C"))
labels_list <- list(c("T3119", "T3119-plus-EC"),c("T3119", "T3119-plus-EC"),
                    c("T3119-plus-EC","T3119-plus-EC"))
cut_off <- c("p_val","p_val_adj","avg_logFC")
fit <- list(); g <- list()
m=1
for(j in 1:length(cut_off)){
        for(i in 1:length(labels_list)){
                fit[[i]] <- eulerr(T_cells_markers, key = labels_list[[i]],
                                     cut_off = cut_off[j],
                                     cut_off_value = 0.05)
                g[[m]] <- plot(fit[[i]], quantities = TRUE, lty = 1:2,
                               fill = colors_list[[i]],
                               legend = F,main =  paste(cut_off[j],": 0.05"))
                m <- m+1
        }    
}
do.call(cowplot::plot_grid, c(g,ncol =3))
dev.off()

##-----------------------------------
jpeg(paste0(path,"/Venn_T_all_pval_adj_logFC.jpeg"), units="in", width=10, height=7,res=600)
g <- lapply(c("p_val","p_val_adj","avg_logFC"), function(cut_off) {
        eulerr(T_cells_markers, cut_off = cut_off, cut_off_value = 0.05) %>%
                plot(quantities = TRUE, lty = 1:2,
                     fill = c("#377EB8","#A6D854","#E41A1C"),
                     labels = list(labels = c("TALL", "TALL+EC2","\n\n\nTAL+EC")), 
                     legend = F,main = paste(cut_off,": 0.05"))
})
do.call(cowplot::plot_grid, c(g,ncol =2))
dev.off()


# T cell and EC
T_EC <- MergeSeurat(T_cells, EC)
T_EC <- SetAllIdent(T_EC, id = "orig.ident")
table(T_EC@ident)

T_EC_markers <- FindAllMarkers.UMI(T_EC,return.thresh = 1,
                                   test.use = "MAST",min.pct = -Inf,
                                   min.cells.gene = -Inf, min.cells.group = -Inf,
                                   only.pos = F, logfc.threshold = -Inf)
write.csv(T_EC_markers,paste0(path,"T_EC_markers.csv"))

#-----

jpeg(paste0(path,"/Venn_TEC_all_p_val.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr(T_EC_markers,shape =  "ellipse",
            cut_off = "p_val", cut_off_value = 0.05) %>%
        plot(quantities = TRUE, lty = 1:6,
             fill = c("#1e601a","#A6D854","#FDC086",
                      "#377EB8","#E41A1C","#6A3D9A"),
             legend = T,main = paste("p_val: 0.05"))
print(g)
dev.off()

#-----
jpeg(paste0(path,"/Venn_TEC_all_p_val_adj.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr(T_EC_markers,shape =  "ellipse",
            cut_off = "p_val_adj", cut_off_value = 0.05) %>%
        plot(quantities = TRUE, lty = 1:6,
             fill = c("#1e601a","#A6D854","#FDC086",
                      "#377EB8","#E41A1C","#6A3D9A"),
             legend = T,main = paste("p_val_adj: 0.05"))
print(g)
dev.off()

#-----
jpeg(paste0(path,"/Venn_TEC_all_avg_logFC.jpeg"), units="in", width=10, height=7,res=600)
g <- eulerr(T_EC_markers,shape =  "ellipse",
            cut_off = "avg_logFC", cut_off_value = 0.01) %>%
        plot(quantities = TRUE, lty = 1:6,
             fill = c("#1e601a","#A6D854","#FDC086",
                      "#377EB8","#E41A1C","#6A3D9A"),
             legend = T,main = paste("avg_logFC: 0.01"))
print(g)
dev.off()

