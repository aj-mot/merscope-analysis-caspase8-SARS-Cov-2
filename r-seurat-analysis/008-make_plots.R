#-------------------------------------------------------------------------------
# Merscope data analysis

library(magrittr)
library(tidyverse)
library(Seurat)
library(Cairo)

set.seed(2025)
options(scipen = 9999)
options(ggrepel.max.overlaps = Inf)
options(future.globals.maxSize= 2097152000)

data_dir <- file.path("..", "..", "merscope-data-caspase-SARS-Cov-2")

#-------------------------------------------------------------------------------

# Load clustered and annotated data
merged_obj <- readRDS(file.path(data_dir, "data", "seurat", "merged_obj", "spatial_annotated.rds"))

# Load processed data with four genes removed: Ifnb1, Fosb, Ptx3, Pgrs2
merged_obj2 <- readRDS(file.path(data_dir, "data", "seurat", "merged_obj", "spatial_niches_2.rds"))

#-------------------------------------------------------------------------------

# Spatial plot of cells in SARS-CoV2 infected regions
# Figures 2J & S2F

fovs <- c("Ctrl.uninfected", "R3_KO.infected", "C8_R3_KO.infected", "WT.infected")
sample_ids <- unique(merged_obj$sample)
sample_names_plot <- c("WT uninfected", "R3 KO infected", "C8 R3 KO infected", "WT infected")

# Plot cells within 10 microns of infected region
covid_colours <- c("red2", "lightblue3")

for (i in 1:4) {
  fov <- fovs[i]
  
  plot_covid_colours <- covid_colours
  if (i == 1) plot_covid_colours <- covid_colours[2] # Uninfected sample
  
  p1 <- merged_obj %>%
    subset(subset = !is.na(in_covid_region)) %>%
    ImageDimPlot(fov = fov, 
                 group.by = "in_covid_region", 
                 cols = plot_covid_colours,
                 dark.background = F) +
    ggtitle(sample_names_plot[i]) +
    theme(legend.title=element_blank())
  
  CairoPDF(file.path("output", "seurat", "covid_regions",
                              paste0(fov, "_covid_image_dim_plot.pdf")), 
                     width = 11.69, height = 8.27)
  print(p1)
  dev.off()
}


#-------------------------------------------------------------------------------

# Spatial plots of cell type clusters

ct_name_colours <- c("#D4DFC6", "#DA9161", "#775BDE", "#C484CF", "#7FE85D", 
                     "#E0D285", "#80A3DD", "#A2E095", "#DAB8CD", "#E25984",
                     "#8BD5E3", "#D74ED9", "#6C887C", "#D8DA4A", "#70E3C9")
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
  "Epithelial alveolar cells", 
  "Fosb+ cells", 
  "Ifnb1+ cells",
  "Ptgs2+ cells", 
  "Ptx3+ cells"
)
names(ct_name_colours) <- cluster_labels

for (i in 1:4) {
  fov <- fovs[i]
  p1 <- ImageDimPlot(merged_obj, 
                     fov = fov, 
                     group.by = "harmony_clusters_r0.2_name", 
                     cols = ct_name_colours,
                     dark.background = F) +
    ggtitle(sample_names_plot[i]) +
    theme(legend.title=element_blank())
  
  CairoPDF(file.path("output", "seurat", "clusters", "15_clusters",
                     paste0(fov, "_harmony_cluster_image_dim_plot.pdf")), 
           width = 11.69, height = 8.27)
  print(p1)
  dev.off()
  
}

# Version with 4 genes omitted
ct_name_colours2 <- ct_name_colours[1:11]

for (i in 1:4) {
  fov <- fovs[i]
  p1 <- ImageDimPlot(merged_obj2, 
                     fov = fov, 
                     group.by = "harmony_clusters_r0.2_name", 
                     cols = ct_name_colours2,
                     dark.background = F) +
      ggtitle(sample_names_plot[i]) +
      theme(legend.title=element_blank())
  
  CairoPDF(file.path("output", "seurat", "clusters", "11_clusters",
                              paste0(fov, "_harmony_cluster_image_dim_plot_v2.pdf")), 
                     width = 11.69, height = 8.27)
  print(p1)
  dev.off()
}

#-------------------------------------------------------------------------------

# UMAP plots

# UMAP by cluster

# Figure S2C - both coords are flipped in paper
CairoPDF(file.path("output", "seurat", "clusters", "15_clusters", "umap_harmony_clusters_r0.2.pdf"), 
    width = 11.69, height = 8.27)
DimPlot(merged_obj, reduction = "umap.integrated",
        group.by = c("harmony_clusters_r0.2_name"), 
        cols = ct_name_colours, 
        label = F) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(plot.title = element_blank())
dev.off()

# Version with 4 genes omitted
CairoPDF(file.path("output", "seurat", "clusters", "11_clusters", "umap_harmony_clusters_r0.2_2.pdf"), 
    width = 11.69, height = 8.27)
DimPlot(merged_obj2, reduction = "umap.integrated",
        group.by = c("harmony_clusters_r0.2_name"), 
        cols = ct_name_colours2, 
        label = F) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(plot.title = element_blank())
dev.off()


# UMAP with all 300 genes
# Colour cells by 11 clusters defined with 296 genes (4 genes omitted)
# 4 previous clusters each driven by 1 gene have their cells reallocated to different cell types

all(Cells(merged_obj) == Cells(merged_obj2)) # TRUE
merged_obj$harmony_clusters_r0.2_name_v2 <- merged_obj2$harmony_clusters_r0.2_name

CairoPDF(file.path("output", "seurat", "clusters", "15_clusters", "umap_harmony_clusters_r0.2_realloc.pdf"), 
         width = 11.69, height = 8.27)
DimPlot(merged_obj, reduction = "umap.integrated",
        group.by = c("harmony_clusters_r0.2_name_v2"), 
        cols = ct_name_colours2, 
        label = F) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(plot.title = element_blank())
dev.off()

# UMAP with 4 genes omitted - 11 clusters
# Plot expression of Ifnb1, Fosb, Ptx3, Pgrs2 which drive the 4 previous clusters 
# that have been reallocated.

# Transfer UMAP from merged_obj2 to merged_obj
merged_obj2_umap <- Embeddings(merged_obj2, reduction = "umap.integrated")
colnames(merged_obj2_umap) <- paste0("UMAP2_", 1:2)
merged_obj[["umap_omit_4_genes"]] <- CreateDimReducObject(embeddings = merged_obj2_umap, 
                                                          key = "UMAP2_", 
                                                          assay = DefaultAssay(merged_obj))

ct_short_names <- c("Epi alveolar", "Pericytes", "Myeloid", "B-cells",                                
                    "Endothelial", "T-cells/NK", "DCs", "T-cells/IM",       
                    "Cil/Deut", "Mesenchymal", "Cil/Deut + Sec")
names(ct_short_names) <- c("Epithelial alveolar cells",               "Pericytes",                              
                           "Myeloid cells",                           "B-cells",                                
                           "Endothelial cells",                       "T-cells/natural killer cells",           
                           "Mature dendritic cells",                  "T-cells/interstitial macrophages",       
                           "Ciliated/Deuterosomal cells",             "Mesenchymal cells",                      
                           "Ciliated/Deuterosomal + Secretory cells")
harmony_clusters_r0.2_short_name_v2 <- ct_short_names[match(as.character(merged_obj$harmony_clusters_r0.2_name_v2), names(ct_short_names))]
names(harmony_clusters_r0.2_short_name_v2) <- Cells(merged_obj)
merged_obj$harmony_clusters_r0.2_short_name_v2 <- factor(harmony_clusters_r0.2_short_name_v2)

Idents(merged_obj) <- "harmony_clusters_r0.2_short_name_v2"
genes <- c("Fosb", "Ptgs2", "Ptx3", "Ifnb1")

# Figure S2D
CairoPDF(file.path("output", "seurat", "clusters", "11_clusters",
                   paste0("Four_genes_expression_umap.pdf")), 
         width = 11.69, height = 8.27)
p <- FeaturePlot(merged_obj, 
                 features = genes, 
                 reduction = "umap_omit_4_genes", 
                 slot = "data", 
                 order = T, 
                 label = T, 
                 raster = F)
p[[1]] <- p[[1]] + xlab("UMAP 1") + ylab("UMAP 2")
p[[2]] <- p[[2]] + xlab("UMAP 1") + ylab("UMAP 2")
p[[3]] <- p[[3]] + xlab("UMAP 1") + ylab("UMAP 2")
p[[4]] <- p[[4]] + xlab("UMAP 1") + ylab("UMAP 2")
print(p)
dev.off()

#-------------------------------------------------------------------------------

# Cluster proportions by sample

a1 <- table(merged_obj2$sample, merged_obj2$harmony_clusters_r0.2_name) %>%
  as.data.frame()

# Order cluster cell type names
a1$cluster_name <- factor(a1$Var2, levels = cluster_labels[1:11], ordered = T)

# Change and order sample names
a1$sample_name <- sample_names_plot[match(a1$Var1, sample_ids)]
a1$sample_name <- factor(a1$sample_name, levels = c("WT uninfected", "WT infected", 
                                                    "R3 KO infected", "C8 R3 KO infected"), 
                         ordered = T)

a <- a1 %>%
  ggplot(aes(x = sample_name, y = Freq, fill = cluster_name)) +
  geom_col(position = position_fill()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = ct_name_colours) +
  theme_classic(base_size = 14) +
  labs(x = "", y = "Cluster Proportions") + 
  theme(legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5))

CairoPDF(file.path("output", "seurat", "clusters", "11_clusters",
                   "cluster_proportion_by_sample.pdf"), width = 8, height = 5)
print(a)
dev.off()


# Cluster proportions by sample and distance to SARS-CoV2 infected region (decile)

sum(is.na(merged_obj2$covid_dist_decile))
sum(is.na(merged_obj2$in_covid_region))
table(merged_obj2$covid_dist_decile, merged_obj2$in_covid_region)
table(is.na(merged_obj2$covid_dist_decile), merged_obj2$in_covid_region)


a1 <- table(merged_obj2$sample, merged_obj$harmony_clusters_r0.2_name_v2, merged_obj$covid_dist_decile) %>%
  as.data.frame()

# Order cluster cell type names
a1$cluster_name <- factor(a1$Var2, levels = cluster_labels[1:11], ordered = T)

sample_names_plot <- c("WT uninfected", "R3 KO infected", "C8 R3 KO infected", "WT infected")
names(sample_names_plot) <- c("Ctrl-uninfected", "R3_KO-infected", "C8_R3_KO-infected", "WT-infected")

# Change and order sample names, drop uninfected sample
a1$sample_name <- as.character(sample_names_plot[match(a1$Var1, names(sample_names_plot))])
a1 <- a1[!(a1$sample_name == sample_names_plot[1]), ]
a1$sample_name <- factor(a1$sample_name, levels = sample_names_plot[-1], ordered = T)

a1$sample_name <- sample_names_plot[match(a1$Var1, sample_ids)]
a1 <- a1[!(a1$sample_name == sample_names_plot[1]), ]
a1$sample_name <- factor(a1$sample_name, 
                         levels = c("WT infected", "R3 KO infected", "C8 R3 KO infected"), 
                         ordered = T)

a <- a1 %>%
  ggplot(aes(x = Var3, y = Freq, fill = cluster_name)) +
  geom_col(position = position_fill()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = ct_name_colours) +
  theme_classic(base_size = 14) +
  labs(x = "Distance to infected region (decile)", y = "Cluster Proportions") + 
  theme(legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(cols = vars(sample_name))

# Figure 2G
CairoPDF(file.path("output", "seurat", "clusters", "11_clusters",
                   "cluster_proportion_by_sample_covid_decile.pdf"), width = 11, height = 5)
print(a)
dev.off()


#-------------------------------------------------------------------------------

# Investigate expression of Ifnb1, Fosb, Ptx3, Pgrs2 in clusters determined with 
# these genes omitted

Idents(merged_obj) <- "harmony_clusters_r0.2_name_v2"

my_exp <- FetchData(merged_obj, genes, layer = "scale.data")
my_exp$sample <- merged_obj$sample
my_exp$Cell_type <- merged_obj$harmony_clusters_r0.2_name_v2
my_exp$covid_decile <- merged_obj$covid_dist_decile
my_exp$covid_decile <- as.character(my_exp$covid_decile)
my_exp$covid_decile[my_exp$sample == "Ctrl-uninfected"] <- "-"
my_exp$covid_decile <- factor(my_exp$covid_decile, levels = c(0:10, "-"), ordered = T)
table(my_exp$covid_decile, my_exp$sample)

# Cropped cells for Covid region detection
table(is.na(my_exp$covid_decile), my_exp$sample)

my_exp <- my_exp[!is.na(my_exp$covid_decile), ]
nrow(my_exp) # 170,616 cells

max(my_exp$Ifnb1[my_exp$sample == "Ctrl-uninfected"])
max(my_exp$Ifnb1)
summary(my_exp$Ifnb1)
table(my_exp$sample)
table(my_exp$sample[my_exp$Ifnb1 > 5])
table(my_exp$sample[my_exp$Ifnb1 > 5], my_exp$covid_decile[my_exp$Ifnb1 > 5])
table(my_exp$sample, my_exp$covid_decile)


# Plot percentage of cells with Ifnb1 > 5 by sample and covid decile - stack cell type
my_high_ifnb1_exp <- my_exp[my_exp$Ifnb1 > 5, ]

my_high_ifnb1_table <- data.frame(table(my_high_ifnb1_exp$sample, my_high_ifnb1_exp$Cell_type, my_high_ifnb1_exp$covid_decile))
my_high_ifnb1_table <- my_high_ifnb1_table[my_high_ifnb1_table$Var3 != "-", ]

# Scale by total counts of each sample/decile
my_exp_table <- data.frame(table(my_exp$sample, my_exp$covid_decile))
colnames(my_exp_table)[2:3] <- c("Var3", "totalCells")
my_exp_table <- my_exp_table[my_exp_table$Var3 != "-", ]
my_high_ifnb1_table <- merge(my_high_ifnb1_table, my_exp_table, all = T)
my_high_ifnb1_table$propHighExpr <- my_high_ifnb1_table$Freq / my_high_ifnb1_table$totalCells


a1 <- my_high_ifnb1_table

# Change and order sample names, drop uninfected sample
a1$sample_name <- as.character(sample_names_plot[match(a1$Var1, names(sample_names_plot))])
a1 <- a1[!(a1$sample_name == sample_names_plot[1]), ]
a1$sample_name <- factor(a1$sample_name, levels = sample_names_plot[c(4, 2, 3)], ordered = T)

# Add and order cluster cell type names
a1$cluster_name <- factor(a1$Var2, levels = cluster_labels[1:11], ordered = T)

a <- a1 %>%
  ggplot(aes(x = Var3, y = propHighExpr, fill = cluster_name)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = ct_name_colours) +
  theme_classic(base_size = 14) +
  ylim(c(0, 0.077)) +
  labs(x = "Distance to infected region (decile)", 
       y = "High Ifnb1 cells (proportion of total cells)") + #, title = cluster_name_var) +
  theme(legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        #axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(cols = vars(sample_name))

# Figure 2L
print(a)


# Repeat for all genes
# Count number of cells with scaled expression > 5
# Tabulate by sample and distance to SARS-CoV2 infection
my_exp_long <- pivot_longer(my_exp, cols = Fosb:Ifnb1, names_to = "Gene", values_to = "expression")
my_high_exp_long <- my_exp_long[my_exp_long$expression > 5, ]

my_high_exp_table <- data.frame(table(my_high_exp_long$sample, my_high_exp_long$Cell_type, my_high_exp_long$covid_decile, my_high_exp_long$Gene))

# Scale by total counts of each sample/decile
my_exp_table <- data.frame(table(my_exp$sample, my_exp$covid_decile))
colnames(my_exp_table)[2:3] <- c("Var3", "totalCells")

my_high_exp_table <- merge(my_high_exp_table, my_exp_table, all = T)
my_high_exp_table$propHighExpr <- my_high_exp_table$Freq / my_high_exp_table$totalCells

a1 <- my_high_exp_table

# Change and order sample names, drop uninfected sample
a1$sample_name <- as.character(sample_names_plot[match(a1$Var1, names(sample_names_plot))])
a1$sample_name <- factor(a1$sample_name, levels = sample_names_plot[c(4, 2:3, 1)], ordered = T)

# Add and order cluster cell type names
a1$cluster_name <- factor(a1$Var2, levels = cluster_labels[1:11], ordered = T)

a1 <- a1[a1$totalCells > 0, ]
max_y <- c(0.11, 0.036, 0.026, 0.078) # Max y-axis value for plots of each gene
for (i in 1:4) {
  
  gene <- genes[i]
  
  a <- a1 %>%
    filter(Var4 == gene) %>%
    ggplot(aes(x = Var3, y = propHighExpr, fill = cluster_name)) +
    geom_col() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = ct_name_colours) +
    theme_classic(base_size = 14) +
    ylim(c(0, max_y[i])) + 
    labs(x = "Distance to infected region (decile)", 
         y = paste("High", gene, "cells (proportion of total cells)")) + #, title = cluster_name_var) +
    theme(legend.title = element_blank(),
          axis.text = element_text(color = "black"),
          #axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(hjust = 0.5)) +
    facet_grid(. ~ sample_name, 
               scale = "free_x", space = "free_x")
  
  # Figure 2L
  CairoPDF(file.path("output", "seurat", "clusters", "11_clusters",
                     paste0("high_", gene, "_covid_distance_decile_by_sample.pdf")), 
           width = 12, height = 5)
  print(a)
  dev.off()
}
