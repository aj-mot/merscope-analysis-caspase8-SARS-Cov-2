#-------------------------------------------------------------------------------
# Preprocess the MERSCOPE data after custom Cellpose segmentation

library(magrittr)
library(tidyverse)
library(sf)
library(sfarrow)
library(BiocParallel)
library(Seurat)
library(readr)
library(arrow)
library(sf)

set.seed(2025)
options(scipen = 9999)
options(ggrepel.max.overlaps = Inf)

# Project data must be downloaded from DOI: 10.5281/zenodo.15719666 and placed 
# at the following location:
data_dir <- file.path("..", "..", "merscope-data-caspase-SARS-Cov-2")

#-------------------------------------------------------------------------------
# MERSCOPE data:
region_names <- c("Ctrl-uninfected", "R3_KO-infected", "C8_R3_KO-infected", "WT-infected")

# Cellpose custom segmenation data
reseg_dir <- file.path(data_dir, "data", "cellpose")
reseg_regions <- paste0("Bader_SARS_region", 0:3)
reseg_indirs <- file.path(reseg_dir, reseg_regions, "VPT_NextFlow", "finalFiles")

## LOAD IN TRANSCRIPT AND METADATA FILES ----
count_files <- sapply(reseg_indirs, function(x) list.files(x, "cell_by_gene_repartitioned.csv", full.names = T))
transcript_files <- sapply(reseg_indirs, function(x) list.files(x, "detected_transcripts.csv", full.names = T))
meta_files <- sapply(reseg_indirs, function(x) list.files(x, "cell_metadata_resegmented.csv", full.names = T))

# Read in files
counts <- lapply(count_files, function(XX) read_csv(XX, col_types = c(cell = "c")))

transcripts <- lapply(transcript_files, function(XX) {
 read_csv(XX, col_types = c(transcript_id = "c", cell_id = "c")) })

metadata <- lapply(meta_files, function(XX) {
  tmp_meta <- read.delim(XX, sep = ",", colClasses = c(EntityID = "character"))
  rownames(tmp_meta) <- tmp_meta$EntityID
  colnames(tmp_meta)[match("EntityID", colnames(tmp_meta))] <- "cell_id"
  tmp_meta })

# Rename files in lists
names(counts) <- region_names
names(transcripts) <- region_names
names(metadata) <- region_names


# Spatial coordinate offsets to allow for joint analysis
coord_adjust_list <- list(
  `Ctrl-uninfected`   = c(-2800, 1500),
  `R3_KO-infected`    = c(500, 4000),
  `C8_R3_KO-infected` = c(2500, 0),
  `WT-infected`       = c(0, -1000))


#### CREATE SEURAT OBJECTS ----
sample_ids <- region_names
obj_list <- list()
obj_list <- sapply(sample_ids, function(XX) {
  message(paste("Creating Seurat object for sample", XX))
  
  # Create a Seurat object containing the RNA info 
  in_counts <- t(counts[[XX]][, -1])
  colnames(in_counts) <- counts[[XX]]$cell
  
  gene_counts <- in_counts[!grepl("Blank", rownames(in_counts)), ]
  blank_counts <- in_counts[grepl("Blank", rownames(in_counts)), ]
  
  sobj <- CreateSeuratObject(counts = gene_counts, 
                             assay = "Vizgen")
  # Assay is "Vizgen", not "RNA"
  
  # Add Blank probe counts
  sobj[["BlankProbes"]] <- CreateAssayObject(counts = blank_counts)
    
  # Add metadata
  sobj <- AddMetaData(sobj, metadata = metadata[[XX]])
  sobj$sample <- XX
  
  # Add centroids fov
  centroids <- data.frame(x = sobj$center_x, y = sobj$center_y, 
                          cell = colnames(x = sobj), 
                stringsAsFactors = FALSE)
  cents <- CreateCentroids(centroids)
  segmentations.data <- list(centroids = cents)
  
  # Add detected transcripts micron-level molecule coordinates
  microns <- data.frame(x = transcripts[[XX]]$global_x, y = transcripts[[XX]]$global_y, 
                gene = transcripts[[XX]]$gene, stringsAsFactors = FALSE)
  
  coords <- CreateFOV(coords = segmentations.data, 
                      type = c("centroids"), 
                      molecules = microns,
                      assay = "Vizgen")

  coords <- subset(x = coords, cells = intersect(x = Cells(x = coords[["centroids"]]), 
        y = Cells(x = sobj)))
  sobj[[XX]] <- coords
  
  # Calculate fraction of blank probes per cell
  sobj$nCount_ALL <- sobj$nCount_Vizgen + sobj$nCount_BlankProbes
  sobj$frac_BlankProbes <- sobj$nCount_BlankProbes / sobj$nCount_ALL
  
  # Remove cells with 0 nCount_Xenium
  sobj <- subset(sobj, subset = nCount_Vizgen != 0)
  
  # Rename cells to add sample ID as prefix
  sobj <- RenameCells(sobj, add.cell.id = XX)
  
  # Adjusted coordinates so sample don't overlap in joint analysis
  sobj$adj_x_centroid <- sobj$center_x + coord_adjust_list[[XX]][1]
  sobj$adj_y_centroid <- -1 * (sobj$center_y + coord_adjust_list[[XX]][2])
  
  # Add spatial coordinates as dimension reduction objects
  position_xy <- cbind(sobj$adj_x_centroid, sobj$adj_y_centroid)
  row.names(position_xy) <- row.names(sobj@meta.data)
  colnames(position_xy) <- c("SP_1", "SP_2")
  sobj[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                       assay = DefaultAssay(sobj))
  obj_list[[XX]] <- sobj
})  

# Merge objects (cannot do spatial DimPlots for this)
merged_spatial_unfiltered <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)])

# Add spatial dimension reduction object separately
position_xy <- cbind(merged_spatial_unfiltered$adj_x_centroid,
                     merged_spatial_unfiltered$adj_y_centroid)
row.names(position_xy) <- row.names(merged_spatial_unfiltered@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
merged_spatial_unfiltered[["sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "SP_", assay = DefaultAssay(merged_spatial_unfiltered))

# View the merged sample coordinates  
DimPlot(merged_spatial_unfiltered, reduction = "sp", group.by = "sample", label = TRUE)

# JoinLayers since this is Seurat v5
merged_spatial_unfiltered <- JoinLayers(merged_spatial_unfiltered)

# Save file
saveRDS(merged_spatial_unfiltered, file.path(data_dir, "data", "seurat", "merged_obj", "spatial_unfiltered.rds"))

#-------------------------------------------------------------------------------