#-------------------------------------------------------------------------------
# Merscope data analysis
# Identification of spatial niches

library(magrittr)
library(tidyverse)
library(Seurat)
library(Cairo)
library(randomcoloR)

set.seed(2025)
options(scipen = 9999)
options(ggrepel.max.overlaps = Inf)
options(future.globals.maxSize= 2097152000)

data_dir <- file.path("..", "..", "merscope-data-caspase-SARS-Cov-2")


# Load Merscope data with annotated clusters with four genes removed: Ifnb1, Fosb, Ptx3, Pgrs2
merged_obj2 <- readRDS(file.path(data_dir, "data", "seurat", "merged_obj", "spatial_annotated_2.rds"))

# Use 11 clusters obtained with 4 genes omitted to find niches
ct_var = merged_obj2$harmony_clusters_r0.2_name
ct_name_var = "harmony_clusters_r0.2_name"

ct_name_colours <- c("#D4DFC6", "#DA9161", "#775BDE", "#C484CF", "#7FE85D", 
                     "#E0D285", "#80A3DD", "#A2E095", "#DAB8CD", "#E25984",
                     "#8BD5E3")
cluster_labels <- c(
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
  "Epithelial alveolar cells"
)

names(ct_name_colours) <- cluster_labels

# Update sample names
sample_names <- c("WT uninfected", "WT infected", "R3 KO infected", "C8 R3 KO infected")
names(sample_names) <- c("Ctrl-uninfected", "WT-infected", "R3_KO-infected", "C8_R3_KO-infected")

# Niche analysis parameters
num_neighbors <- 30
num_niches <- 10

# Set directory names for saving files
niche_name_var <-  paste0(ct_name_var, "_nichek", num_niches, "_neighborsk",
                          num_neighbors, "_2")
plot_save_dir <- file.path("output", "seurat", "niches", niche_name_var)
if (!dir.exists(plot_save_dir)) dir.create(plot_save_dir)

# Create FOV
coords <- merged_obj2@meta.data[, c("adj_x_centroid", "adj_y_centroid")]
merged_obj2[["fov"]]  <- CreateFOV(coords, assay = "Xenium", type = "centroids")

# Run Niche Analysis
set.seed(2025)
merged_obj2 <- BuildNicheAssay(object = merged_obj2, fov = "fov", group.by = ct_name_var,
                               niches.k = num_niches, neighbors.k = num_neighbors)

niche_var <- merged_obj2$niches
merged_obj2@meta.data[[niche_name_var]] <- merged_obj2$niches

# Rename and order niches
niche_names <- paste0("N", 1:10)
niche_order <- c(1, 10, 9, 6, 7, 5, 8, 2, 3, 4)
niche_name <- factor(niche_names[match(niche_var, niche_order)], 
                     levels = niche_names, 
                     ordered = T)

niche_colours <- c("#D4DFC6", "#DA9161", "#7FE85D", "#775BDE", "#C484CF",
                   "#E0D285", "#8BD5E3", "#A2E095", "#DAB8CD", "#E25984")
names(niche_colours) <- niche_names

# Plot niches

# Niche composition by cell type
# Figure 2H
CairoPDF(file.path(plot_save_dir, paste0(niche_name_var, "_by_ct.pdf")), width = 10, height = 6)
a0 <- table(ct_var, niche_name) %>%
  as.data.frame()
a0$ct_var <- as.character(a0$ct_var)
a0$ct_var <- factor(a0$ct_var, levels = cluster_labels, ordered = T)

a <- a0 %>%
  ggplot(aes(x = niche_name, y = Freq, fill = ct_var)) +
  geom_col(position = position_fill()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = ct_name_colours) +
  labs(x = "Cell Niche", y = "Cell Type Proportion", fill = "", title = niche_name_var) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_blank())
print(a)
dev.off()

# Niche proportion by sample
CairoPDF(file.path(plot_save_dir, paste0(niche_name_var, "_by_sample.pdf")), width = 8, height = 5)
a1 <- table(merged_obj2$sample, niche_name) %>%
  as.data.frame()
# Change and order sample names
a1$sample_name <- sample_names[match(as.character(a1$Var1), names(sample_names))]
a1$sample_name <- factor(a1$sample_name, levels = sample_names, ordered = T)

a <- a1 %>%
  ggplot(aes(x = sample_name, y = Freq, fill = niche_name)) +
  geom_col(position = position_fill()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = niche_colours) +
  theme_classic(base_size = 14) +
  labs(x = "", y = "Niche Proportions", title = niche_name_var, fill = "Niche") +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_blank())
print(a)
dev.off()

# Niche by proportion by sample within Covid regions
# Figure 2I
CairoPDF(file.path(plot_save_dir, paste0(niche_name_var, "_covid_regions_by_sample.pdf")), width = 8, height = 5)
a1 <- table(merged_obj2$sample, niche_name, merged_obj2$in_covid_region) %>%
  as.data.frame()
# Change and order sample names
a1$sample_name <- sample_names[match(as.character(a1$Var1), names(sample_names))]
a1$sample_name <- factor(a1$sample_name, levels = sample_names, ordered = T)
levels(a1$Var3) <- c("SARS-CoV2 infected region", "Non-infected region")
a <- a1 %>%
  ggplot(aes(x = sample_name, y = Freq, fill = niche_name)) +
  geom_col(position = position_fill()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = niche_colours) +
  theme_classic(base_size = 14) +
  labs(x = "", y = "Niche Proportions", title = niche_name_var, fill = "Niche") +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_blank()) +
  facet_wrap(vars(Var3))
print(a)
dev.off()


# Spatial plot of niches in each sample
# Figures 2K & S2F
merged_obj2$niche_names <- niche_name
fovs <- c("Ctrl.uninfected", "WT.infected", "R3_KO.infected", "C8_R3_KO.infected")
for (k in 1:4) {
  fov <- fovs[k]
  p1 <- ImageDimPlot(merged_obj2, 
                     fov = fov, 
                     group.by = "niche_names", 
                     cols = niche_colours,
                     dark.background = F) +
    ggtitle(sample_names[k]) +
    labs(fill = "Niche")
  
  CairoPDF(file.path(plot_save_dir, paste0(fov, "_niche_image_dim_plot.pdf")), 
           width = 11.69, height = 8.27)
  print(p1)
  dev.off()
}

#-------------------------------------------------------------------------------
# Save data
saveRDS(merged_obj2, file.path(data_dir, "data", "seurat", "merged_obj", "spatial_niches_2.rds"))
