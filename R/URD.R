library(Seurat)
library(URD)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 6.1 convert from seurat to URD  ==========================================
# load Seurat
(load(file = "data/Lymphoma_EC_12_20190125.Rda"))
colors <- ExtractMetaColor(object)
URD_object <- seuratToURD(object)
remove(object);GC()
pcSDPlot(URD_object)
plotDim(URD_object,"singler2sub",plot.title = "tSNE: Stage",
        discrete.colors=colors,legend = F)

#======= 6.2 Calculate Diffusion Map =============
# In this case, knn=100 (larger than sqrt(n.cells)) works well because there are not many cell types.
# Sigma 16 is slightly smaller than the sigma auto-determined by using NULL parameter.
system.time(URD_object <- calcDM(URD_object, knn = 100, sigma=16))
plotDimArray(URD_object, reduction.use = "dm", dims.to.plot = 1:8, 
             outer.title = "Diffusion Map (Sigma 16, 100 NNs): Stage", 
             label="singler2sub", plot.title="", legend=F)
jpeg(paste0(path,"Developmental_stage.jpeg"), units="in", width=10, height=7,
     res=600)
plotDim(URD_object, "singler2sub", transitions.plot = 10000, 
        plot.title="Developmental stage (with transitions)",
        discrete.colors=colors,legend = F)+
        theme(text = element_text(size=15),			
              plot.title = element_text(hjust = 0.5))
dev.off()
# ========6.3 Calculate pseudotime ===================
# Here we use all cells from the first stage as the root
root.cells <- cellsInCluster(URD_object, "orig.ident", c("EC-only","EC"))

# Then we run 'flood' simulations
URD_floods <- floodPseudotime(URD_object, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

# The we process the simulations into a pseudotime
URD_object <- floodPseudotimeProcess(URD_object, URD_floods, floods.name="pseudotime")

pseudotimePlotStabilityOverall(URD_object)

jpeg(paste0(path,"pseudotime.jpeg"), units="in", width=10, height=7,
     res=600)
plotDim(URD_object, "pseudotime")
dev.off()

plotDists(URD_object, "pseudotime", "singler2sub", plot.title="Pseudotime by stage")

# =======6.4 Find tips =======================
# Create a subsetted object of just those cells from the final stage
axial.6somite <- urdSubset(URD_object, cells.keep=cellsInCluster(URD_object, "stage", "L-6S"))

# Use the variable genes that were calculated only on the final group of stages (which
# contain the last stage).
axial.6somite@var.genes <- var.by.stage[[4]]

# Calculate PCA and tSNE
axial.6somite <- calcPCA(axial.6somite, mp.factor = 1.5)
pcSDPlot(axial.6somite)
