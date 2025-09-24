#-------------------------------------------------------------------------------
# Merscope data analysis
# Add SARS-CoV-2 infection region metadata

library(magrittr)
library(tidyverse)
library(Seurat)
library(reticulate)
library(anndata)
library(dplyr)

set.seed(2025)
options(scipen = 9999)
options(ggrepel.max.overlaps = Inf)
options(future.globals.maxSize= 2097152000)

data_dir <- file.path("..", "..", "merscope-data-caspase-SARS-Cov-2")

#-------------------------------------------------------------------------------

# Load Merscope data with annotated clusters
merged_obj <- readRDS(file.path(data_dir, "data", "seurat", "merged_obj", "spatial_annotated.rds"))

# Load Merscope data with annotated clusters with four genes removed: Ifnb1, Fosb, Ptx3, Pgrs2
merged_obj2 <- readRDS(file.path(data_dir, "data", "seurat", "merged_obj", "spatial_annotated_2.rds"))

#-------------------------------------------------------------------------------
# Load data with cell distances to Covid infected regions
region_names <- c("Ctrl-uninfected", "R3_KO-infected", "C8_R3_KO-infected", "WT-infected")

# Files in AnnData format with distance from each cell to nearest SARS-CoV-2 infected region
covid_dist_dir <- file.path(data_dir, "data", "cellpose", "cell_to_infected_region_distances")
covid_anndata_files <- sapply(0:3, function(i) list.files(covid_dist_dir, paste0(i, ".h5ad"), full.names = T))

# Set up a vector with distance to nearest SARS-CoV-2 protein (microns)
distance_to_nearest_covid <- rep(NA, nrow(merged_obj[[]]))
names(distance_to_nearest_covid) <- rownames(merged_obj[[]])

for (i in 2:4) {

  region <- region_names[i]
  covid_anndata <- read_h5ad(covid_anndata_files[[i]])
  
  this.distance_to_nearest_covid <- covid_anndata[["obs"]]$distance_to_nearest_covid
  names(this.distance_to_nearest_covid) <- paste(region, rownames(covid_anndata[["obs"]]), sep = "_")
  
  # Some cells have been filtered - don't try to match these in the Seurat object
  keep_cell.i <- names(this.distance_to_nearest_covid) %in% names(distance_to_nearest_covid)
  distance_to_nearest_covid[match(names(this.distance_to_nearest_covid[keep_cell.i]), names(distance_to_nearest_covid))] <- this.distance_to_nearest_covid[keep_cell.i]
}
merged_obj[[]]$distance_to_nearest_covid <- distance_to_nearest_covid

# Cells are identical in merged_obj and merged_obj2
all(Cells(merged_obj) == Cells(merged_obj2))
merged_obj2[[]]$distance_to_nearest_covid <- distance_to_nearest_covid

# Cells outside the cropped region have distance set to NA
table(merged_obj[[]]$sample, is.na(merged_obj[[]]$distance_to_nearest_covid))

table(merged_obj2[[]]$sample, is.na(merged_obj2[[]]$distance_to_nearest_covid))

#-------------------------------------------------------------------------------

# Bin cells by distance to nearest Covid
covid.dist.deciles <- quantile(merged_obj[[]]$distance_to_nearest_covid[merged_obj[[]]$distance_to_nearest_covid > 0], 
                               probs = seq(0, 1, 0.1), na.rm = T)
covid_dist_decile <- cut(merged_obj[[]]$distance_to_nearest_covid, 
                                        breaks = c(-Inf, covid.dist.deciles), 
                                        labels = 0:10)
merged_obj[[]]$covid_dist_decile <- covid_dist_decile
merged_obj2[[]]$covid_dist_decile <- covid_dist_decile
# Cells inside infected regions are labelled 0
# Cells outside regions are binned into 10 deciles

table(merged_obj[[]]$covid_dist_decile)
table(merged_obj[[]]$covid_dist_decile, merged_obj[[]]$sample)

# Compare the number of cells close to Covid regions between samples
ggplot(subset(merged_obj[[]], !is.na(covid_dist_decile)), 
       aes(x = covid_dist_decile, fill = sample)) +
  geom_bar(position = position_dodge())

ggplot(subset(merged_obj[[]], !is.na(covid_dist_decile)), 
       aes(x = covid_dist_decile, y = distance_to_nearest_covid, fill = sample)) + 
  geom_violin()

# Categorise cells by whether they are within 10 microns of an infected region
merged_obj[[]]$in_covid_region <- merged_obj[[]]$distance_to_nearest_covid <= 10
merged_obj[[]]$in_covid_region[is.na(merged_obj[[]]$in_covid_region) & merged_obj[[]]$sample == "Ctrl-uninfected"] <- F
merged_obj[[]]$in_covid_region <- c("Cell within 10 um of SARS-CoV2 infected region", "Non-infected region")[2 - as.integer(merged_obj[[]]$in_covid_region)]
merged_obj[[]]$in_covid_region <- factor(merged_obj[[]]$in_covid_region, levels = c("Cell within 10 um of SARS-CoV2 infected region", "Non-infected region"), ordered = T)

merged_obj2[[]]$in_covid_region <- merged_obj[[]]$in_covid_region

#-------------------------------------------------------------------------------

# Save data
saveRDS(merged_obj, file.path(data_dir, "data", "seurat", "merged_obj", "spatial_annotated.rds"))
saveRDS(merged_obj2, file.path(data_dir, "data", "seurat", "merged_obj", "spatial_annotated_2.rds"))
