#-------------------------------------------------------------------------------
# Merscope data cluster annotation

library(magrittr)
library(tidyverse)
library(Seurat)

set.seed(2025)
options(scipen = 9999)
options(ggrepel.max.overlaps = Inf)
options(future.globals.maxSize= 2097152000)

data_dir <- file.path("..", "..", "merscope-data-caspase-SARS-Cov-2")

#-------------------------------------------------------------------------------

# Load clustered data
merged_obj <- readRDS(file.path(data_dir, "data", "seurat", "merged_obj", "spatial_integrated.rds"))

# Load clustered data with four genes removed: Ifnb1, Fosb, Ptx3, Pgrs2
merged_obj2 <- readRDS(file.path(data_dir, "data", "seurat", "merged_obj", "spatial_integrated_2.rds"))

#-------------------------------------------------------------------------------
# Identify cluster markers

# Set identity classes
Idents(merged_obj) <- "harmony_clusters_r0.2"

ident.names <- unique(Idents(merged_obj))

# Find all markers genes for clusters
all_markers <- FindAllMarkers(merged_obj, 
                              assay = "SCT", 
                              only.pos = T)
saveRDS(all_markers, file.path(data_dir, "data", "seurat", "markers", "harmony_clusters_r0.2_all_markers.rds"))

#-------------------------------------------------------------------------------
# Plot all cluster marker genes to compare clusters
pdf(file.path("output", "seurat", "dotplots", paste0("harmony_r0.2_all_cluster_marker_gene_dotplot.pdf")),
      width = 11.69, height = 8.27)
for (i in seq_along(ident.names)) {
  
  ident.name <- ident.names[i]
  markers.to.plot <- all_markers$gene[all_markers$cluster == ident.name]
  
  print(DotPlot(merged_obj, features = markers.to.plot, 
                cluster.idents = T) +
        RotatedAxis() + ggtitle(paste("Cluster", ident.name)))
}
dev.off()

#-------------------------------------------------------------------------------
# Load annotated scRNA-seq reference data from:
# Guo, M., Morley, M.P., Jiang, C. et al. Guided construction of single cell 
# reference for human and mouse lung. Nat Commun 14, 4566 (2023). 
# https://doi.org/10.1038/s41467-023-40173-5
# Downloaded from https://www.lungmap.net/dataset/?experiment_id=LMEX0000004397
lungmap.data <- readRDS(file.path(data_dir, "extdata", "LungMAP_MouseLung_CellRef.v1.1_seurat.rds"))
Idents(lungmap.data) <- "celltype_level3"
lungmap.data <- NormalizeData(lungmap.data)

lungmap_celltypes <- unique(lungmap.data[[]]$celltype_level3)
lungmap_celltypes_fullname <- lungmap.data[[]][match(lungmap_celltypes, lungmap.data[[]]$celltype_level3), 
                                               c("celltype_level3_fullname", "lineage_level1", "lineage_level2")]
# Full names of the annotated cell types
lungmap_celltypes.df <- data.frame(lungmap_celltypes, lungmap_celltypes_fullname)


panel.genes <- Features(merged_obj)

# 1 panel gene is not found in the reference data:
# "Slfn10-ps" is a pseudogene # Not found in reference data
"Slfn10-ps" %in% Features(lungmap.data)

# 7 panel genes are known by other names in the reference data:
# Names in the Merscope data
panel.genes.not.found <- c("Adgre1", "Fcmr", "Mturn", "Siglecf", "Ifit1bl1", "Dipk2b", "Acod1")
panel.genes.not.found %in% Features(lungmap.data)

# Names in the reference data
panel.genes.aka <- c("Emr1", "Faim3", "2410066E13Rik", "Siglec5", "Gm14446", "4930578C19Rik", "Irg1")
names(panel.genes.aka) <- panel.genes.not.found
panel.genes.aka %in% Features(lungmap.data)


#-------------------------------------------------------------------------------
# Plot unbiased cluster marker genes in the reference data
pdf(file.path("output", "seurat", "dotplots", paste0("lungmap_harmony_r0.2_all_cluster_marker_gene_dotplot.pdf")),
      width = 11.69, height = 8.27)
for (i in seq_along(ident.names)) {
  
  ident.name <- ident.names[i]
  markers.to.plot <- all_markers$gene[all_markers$cluster == ident.name]#[1:3]
  
  markers.to.plot <- markers.to.plot[!(markers.to.plot %in% "Slfn10-ps")]
  replace.marker.name.i <- markers.to.plot %in% panel.genes.not.found
  markers.to.plot[replace.marker.name.i] <- panel.genes.aka[match(markers.to.plot[replace.marker.name.i], panel.genes.not.found)]
  
  print(DotPlot(lungmap.data, features = markers.to.plot, 
                cluster.idents = T) +
        RotatedAxis() + ggtitle(paste("Cluster", ident.name)))

}
dev.off()

#-------------------------------------------------------------------------------
# Identify marker genes of reference cell types considering only genes in the Merscope panel
lm_features <- c(panel.genes, panel.genes.aka)[c(panel.genes, panel.genes.aka) %in% Features(lungmap.data)]
lungMAP_markers <- FindAllMarkers(lungmap.data, 
                                  features = lm_features,
                                  only.pos = T)

# Plot the reference cell type marker genes in the Merscope data
pdf(file.path("output", "seurat", "dotplots", paste0("harmony_r0.2_ref_marker_gene_dotplot.pdf")),
      width = 11.69, height = 8.27)
for (i in seq_along(lungmap_celltypes)) {
  
  ident.name <- lungmap_celltypes[i]
  markers.to.plot <- lungMAP_markers$gene[lungMAP_markers$cluster == ident.name]
  
  replace.marker.name.i <- markers.to.plot %in% panel.genes.aka
  markers.to.plot[replace.marker.name.i] <- panel.genes.not.found[match(markers.to.plot[replace.marker.name.i], panel.genes.aka)]
  
  print(DotPlot(merged_obj, features = markers.to.plot, 
                cluster.idents = T) +
        RotatedAxis() + ggtitle(ident.name))

}
dev.off()
#-------------------------------------------------------------------------------
# Add cluster names to metadata

harmony_clusters_r0.2_name <- c(
  "Endothelial cells", 
  "Myeloid cells", 
  "T-cells/natural killer cells", 
  "Epithelial alveolar cells", 
  "Pericytes", 
  "B-cells", 
  "Ciliated/Deuterosomal cells", 
  "Mesenchymal cells", 
  "Fosb+ cells", 
  "Ptx3+ cells", 
  "Ptgs2+ cells", 
  "Mature dendritic cells", 
  "Ciliated/Deuterosomal + Secretory cells", 
  "T-cells/interstitial macrophages", 
  "Ifnb1+ cells"
)

# Note order of cluster levels
merged_obj$harmony_clusters_r0.2[1:10]

merged_obj$harmony_clusters_r0.2_name <- factor(harmony_clusters_r0.2_name[c(0:1, 10:14, 2:9)[as.integer(merged_obj$harmony_clusters_r0.2)] + 1])
merged_obj$harmony_clusters_r0.2_name <- factor(merged_obj$harmony_clusters_r0.2_name, levels = c(
  "Myeloid cells", 
  "T-cells/natural killer cells", 
  "T-cells/interstitial macrophages", 
  "B-cells", 
  "Mature dendritic cells", 
  "Endothelial cells", 
  "Pericytes", 
  "Mesenchymal cells", 
  "Ciliated/Deuterosomal cells", 
  "Ciliated/Deuterosomal + Secretory cells", 
  "Epithelial alveolar cells", 
  "Fosb+ cells", 
  "Ifnb1+ cells",
  "Ptgs2+ cells", 
  "Ptx3+ cells"), ordered = T
)


table(merged_obj$harmony_clusters_r0.2)
table(merged_obj$harmony_clusters_r0.2_name)

# Plot UMAP with cell type labels
Idents(merged_obj) <- "harmony_clusters_r0.2_name"
DimPlot(merged_obj, reduction = "umap.integrated",
        group.by = c("harmony_clusters_r0.2_name"), 
        label = T)

#-------------------------------------------------------------------------------
# Repeat analysis of cluster marker genes with Ifnb1, Ptx3, Ptgs2 and Fosb removed.
#-------------------------------------------------------------------------------

# Set identity classes
Idents(merged_obj2) <- "harmony_clusters_r0.2"

ident.names2 <- unique(Idents(merged_obj2))

# Find all markers genes for clusters
all_markers_2 <- FindAllMarkers(merged_obj2, 
                              assay = "SCT", 
                              only.pos = T)
saveRDS(all_markers_2, file.path(data_dir, "data", "seurat", "markers", "harmony_clusters_r0.2_v2_all_markers.rds"))

#-------------------------------------------------------------------------------
# Plot all cluster marker genes to compare clusters
pdf(file.path("output", "seurat", "dotplots", paste0("harmony_r0.2_all_cluster_marker_gene_dotplot_v2.pdf")),
      width = 11.69, height = 8.27)
for (i in seq_along(ident.names2)) {
  
  ident.name <- ident.names2[i]
  markers.to.plot <- all_markers_2$gene[all_markers_2$cluster == ident.name]
  
  print(DotPlot(merged_obj2, features = markers.to.plot, 
                cluster.idents = T) +
        RotatedAxis() + ggtitle(paste("Cluster", ident.name)))
}
dev.off()

#-------------------------------------------------------------------------------
# Plot unbiased cluster marker genes in the reference data
pdf(file.path("output", "seurat", "dotplots", paste0("lungmap_harmony_r0.2_all_cluster_marker_gene_dotplot_v2.pdf")),
      width = 11.69, height = 8.27)
for (i in seq_along(ident.names2)) {
  
  ident.name <- ident.names2[i]
  markers.to.plot <- all_markers_2$gene[all_markers_2$cluster == ident.name]
  
  markers.to.plot <- markers.to.plot[!(markers.to.plot %in% "Slfn10-ps")]
  replace.marker.name.i <- markers.to.plot %in% panel.genes.not.found
  markers.to.plot[replace.marker.name.i] <- panel.genes.aka[match(markers.to.plot[replace.marker.name.i], panel.genes.not.found)]
  
  print(DotPlot(lungmap.data, features = markers.to.plot, 
                cluster.idents = T) +
        RotatedAxis() + ggtitle(paste("Cluster", ident.name)))

}
dev.off()

#-------------------------------------------------------------------------------
# Plot the reference cell type marker genes in the Merscope data
pdf(file.path("output", "seurat", "dotplots", paste0("harmony_r0.2_ref_marker_gene_dotplot_v2.pdf")),
      width = 11.69, height = 8.27)
for (i in seq_along(lungmap_celltypes)) {
  
  ident.name <- lungmap_celltypes[i]
  markers.to.plot <- lungMAP_markers$gene[lungMAP_markers$cluster == ident.name]
  
  replace.marker.name.i <- markers.to.plot %in% panel.genes.aka
  markers.to.plot[replace.marker.name.i] <- panel.genes.not.found[match(markers.to.plot[replace.marker.name.i], panel.genes.aka)]
  
  print(DotPlot(merged_obj2, features = markers.to.plot, 
                cluster.idents = T) +
        RotatedAxis() + ggtitle(ident.name))

}
dev.off()
#-------------------------------------------------------------------------------
# Add cluster names to metadata

harmony_clusters_r0.2_name_v2 <- c(
  "Endothelial cells", 
  "Myeloid cells", 
  "T-cells/natural killer cells", 
  "Epithelial alveolar cells", 
  "Pericytes", 
  "Ciliated/Deuterosomal cells", 
  "B-cells", 
  "Mesenchymal cells",
  "Mature dendritic cells", 
  "Ciliated/Deuterosomal + Secretory cells", 
  "T-cells/interstitial macrophages", 
  "T-cells/interstitial macrophages"
)


# Note order of cluster levels
merged_obj2$harmony_clusters_r0.2[1:10]
merged_obj2$harmony_clusters_r0.2_name <- factor(harmony_clusters_r0.2_name_v2[c(0:1, 10:11, 2:9)[as.integer(merged_obj2$harmony_clusters_r0.2)] + 1])
merged_obj2$harmony_clusters_r0.2_name <- factor(merged_obj2$harmony_clusters_r0.2_name, levels = c(
  "Myeloid cells", 
  "T-cells/natural killer cells", 
  "T-cells/interstitial macrophages", 
  "B-cells", 
  "Mature dendritic cells", 
  "Endothelial cells", 
  "Pericytes", 
  "Mesenchymal cells", 
  "Ciliated/Deuterosomal cells", 
  "Ciliated/Deuterosomal + Secretory cells", 
  "Epithelial alveolar cells"), ordered = T
)

table(merged_obj2$harmony_clusters_r0.2)
table(merged_obj2$harmony_clusters_r0.2_name)

# Plot UMAP with cell type labels
Idents(merged_obj2) <- "harmony_clusters_r0.2_name"
DimPlot(merged_obj2, reduction = "umap.integrated",
        group.by = c("harmony_clusters_r0.2_name"), 
        label = T)

#-------------------------------------------------------------------------------
# Dot plots of cell type cluster marker genes
# Plot top 3 markers for each cluster 

Idents(merged_obj) <- "harmony_clusters_r0.2_name"
DefaultAssay(merged_obj) <- "SCT"
cell_type_order <- c("1", "2", "13", "5", "11", "0", "4", "7", "6", "12", "3", "8", "14", "10", "9")

  
markers.to.plot <- unique(do.call(c, lapply(cell_type_order, function(x) all_markers$gene[all_markers$cluster == x][1:3])))

pdf(file.path("output", "seurat", "dotplots", paste0("harmony_r0.2_cluster_marker_dotplot_top3_genes.pdf")),
      width = 11.69, height = 8.27)
print(DotPlot(merged_obj, 
              features = markers.to.plot) +
        RotatedAxis() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)))
dev.off()


# Version with 4 genes omitted
Idents(merged_obj2) <- "harmony_clusters_r0.2_name"
DefaultAssay(merged_obj2) <- "SCT"

cell_type_order_2 <- c("1", "2", "10", "6", "8", "0", "4", "7", "5", "9", "3")

markers.to.plot2 <- unique(do.call(c, lapply(cell_type_order_2, function(x) all_markers_2$gene[all_markers_2$cluster == x][1:3])))

# Figure S2E
pdf(file.path("output", "seurat", "dotplots", paste0("harmony_r0.2_cluster_marker_dotplot_top3_genes_v2.pdf")),
      width = 11.69, height = 8.27)
print(DotPlot(merged_obj2, features = markers.to.plot2) +
        RotatedAxis() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)))
dev.off()


#-------------------------------------------------------------------------------
# Save data
saveRDS(merged_obj, file.path(data_dir, "data", "seurat", "merged_obj", "spatial_annotated.rds"))
saveRDS(merged_obj2, file.path(data_dir, "data", "seurat", "merged_obj", "spatial_annotated_2.rds"))
