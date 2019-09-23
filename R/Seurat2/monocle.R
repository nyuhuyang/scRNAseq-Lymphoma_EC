#Identification of important genes
list.of.packages <- c("devtools","dplyr","pheatmap","VGAM", "irlba",
                      "matrixStats", "igraph", "combinat", "fastICA",
                      "grid", "reshape2", "plyr", "parallel", "methods")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#check package
library("devtools")
library("DDRTree")
library("pheatmap")
library("M3Drop")
library("monocle")
library("reshape")
source("../R/Seurat_functions.R")
# 5.1 Importing data from Seurat object=================
(load(file="./output/CancerCell_20181024.RData"))
table(EC@ident)
CancerCell_Mono <- importCDS(EC, import_all = TRUE)

# 5.1.1 Estimate size factors and dispersions
# estimateSizeFactors() and estimateDispersions() will only work,
# and are only needed, if you are working with a CellDataSet 
# with a negbinomial() or negbinomial.size() expression family.
CancerCell_Mono <- estimateSizeFactors(CancerCell_Mono)
CancerCell_Mono <- estimateDispersions(CancerCell_Mono)

# 5.1.2 Filtering low-quality cells
CancerCell_Mono <- detectGenes(CancerCell_Mono, min_expr = 0.1)
print(head(fData(CancerCell_Mono)))
print(head(pData(CancerCell_Mono)))

# 5.1.3 If you are using RPC values to measure expresion, 
# as we are in this vignette, it's also good to look at the distribution
# of mRNA totals across the cells:
pData(CancerCell_Mono)$Total_mRNAs <- Matrix::colSums(exprs(CancerCell_Mono))
upper_bound <- 10^(mean(log10(pData(CancerCell_Mono)$Total_mRNAs)) + 
                     2*sd(log10(pData(CancerCell_Mono)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(CancerCell_Mono)$Total_mRNAs)) - 
                     2*sd(log10(pData(CancerCell_Mono)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(CancerCell_Mono), color = conditions, geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

# 5.2 Classifying and counting cells of different types
# 5.2.1 Classifying cells with CellTypeHierarchy
CancerCell <- SetAllIdent(CancerCell, id = "old.ident")
all(row.names(pData(CancerCell_Mono)) == names(EC@ident))
pData(CancerCell_Mono)$CellType <- EC@ident
table(pData(CancerCell_Mono)$CellType)
pie <- ggplot(pData(CancerCell_Mono), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# 5.3 Constructing Single Cell Trajectories
# 5.3.1 choosing genes that define progress
expressed_genes <- row.names(subset(fData(CancerCell_Mono), num_cells_expressed >= 10))
length(expressed_genes)
diff_test_res <- differentialGeneTest(CancerCell_Mono[expressed_genes,],
                                      fullModelFormulaStr = "~ CellType",
                                      cores = 4) #takes long time
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
CancerCell_Mono <- setOrderingFilter(CancerCell_Mono, ordering_genes)
plot_ordering_genes(CancerCell_Mono)

# 5.3.2 reduce data dimensionality
#Now we're ready to try clustering the cells:.
plot_pc_variance_explained(CancerCell_Mono, return_all = F) # norm_method = 'log',
CancerCell_Mono <- reduceDimension(CancerCell_Mono, max_components = 2,
                                  method = 'DDRTree')
#Trajectory step 3: order cells along the trajectory
CancerCell_Mono <- orderCells(CancerCell_Mono)

g1 <- plot_cell_trajectory(CancerCell_Mono, color_by = "CellType",cell_size = 3)
g2 <- plot_cell_trajectory(CancerCell_Mono, color_by = "State",cell_size = 3)

jpeg(paste0(path,"/trajectory_EC.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(g1,g2)
dev.off()
