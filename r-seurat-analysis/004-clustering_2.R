#-------------------------------------------------------------------------------
# Merscope data clustering

library(magrittr)
library(tidyverse)
library(Seurat)

set.seed(2025)
options(scipen = 9999)
options(ggrepel.max.overlaps = Inf)
options(future.globals.maxSize= 2097152000)

data_dir <- file.path("..", "..", "merscope-data-caspase-SARS-Cov-2")

#-------------------------------------------------------------------------------
# Perform clustering with 4 genes omitted: Ifnb1, Fosb, Ptx3, Ptgs2

# Load filtered data
merged_obj <- readRDS(file.path(data_dir, "data", "seurat", "merged_obj", "spatial_filtered.rds"))

# Remove 4 genes
features2keep <- Features(merged_obj)
features2keep <- features2keep[!(features2keep %in% c("Ifnb1", "Fosb", "Ptx3", "Ptgs2"))]
merged_obj <- subset(merged_obj, features = features2keep)

# 4 samples are integrated
table(merged_obj$sample)

# Split the samples to perform normalization and PCA without integration
merged_obj[["Vizgen"]] <- split(merged_obj[["Vizgen"]], f = merged_obj$sample)
merged_obj <- SCTransform(merged_obj, assay = "Vizgen")
merged_obj <- RunPCA(merged_obj, npcs = 30, features = rownames(merged_obj))
merged_obj <- RunUMAP(merged_obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

#-------------------------------------------------------------------------------
# Perform integration with Harmony
merged_obj <- IntegrateLayers(object = merged_obj, 
                              assay = "SCT", 
                              method = HarmonyIntegration, 
                              orig.reduction = "pca", 
                              new.reduction = "integrated.harmony")

merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.harmony", 
                            dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.2, cluster.name = "harmony_clusters_r0.2")
merged_obj <- RunUMAP(merged_obj, dims = 1:30, reduction = "integrated.harmony", 
                      reduction.name = "umap.integrated")

# Plot clusters
DimPlot(merged_obj, reduction = "umap.integrated",
        group.by = c("sample"))
DimPlot(merged_obj, reduction = "umap.integrated",
        group.by = "harmony_clusters_r0.2", 
        label = T)

fovs <- c("Ctrl.uninfected", "R3_KO.infected", "C8_R3_KO.infected", "WT.infected")
for (i in 1:4) {
  fov <- fovs[i]
  print(ImageDimPlot(merged_obj, fov = fov, group.by = "harmony_clusters_r0.2",
             cols = "polychrome", size = 0.75))
}

# Check cluster proportions by sample
sample.cluster.table <- table(merged_obj$sample, merged_obj$harmony_clusters_r0.2)
round(sweep(sample.cluster.table, 1, rowSums(sample.cluster.table), FUN="/") * 100, digits=2)

# Now that integration is complete, rejoin layers
merged_obj[["Vizgen"]] <- JoinLayers(merged_obj[["Vizgen"]])

# Set slots to "Vizgene" to allow PrepSCTFindMarkers() to run
# Only the first is "Vizgen", the rest are "RNA"
for (i in 1:length(SCTResults(object=merged_obj, slot="umi.assay"))) {
  slot(object = merged_obj@assays$SCT@SCTModel.list[[i]], name="umi.assay") <- "Vizgen"
}
merged_obj <- PrepSCTFindMarkers(merged_obj)

# Save integrated clustered data
saveRDS(merged_obj, file.path(data_dir, "data", "seurat", "merged_obj", "spatial_integrated_2.rds"))
