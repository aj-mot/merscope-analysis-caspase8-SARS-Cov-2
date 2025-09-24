#-------------------------------------------------------------------------------
# Merscope QC filtering

library(magrittr)
library(tidyverse)
library(Seurat)

set.seed(2025)
options(scipen = 9999)
options(ggrepel.max.overlaps = Inf)

data_dir <- file.path("..", "..", "merscope-data-caspase-SARS-Cov-2")

#-------------------------------------------------------------------------------

# Load unfiltered data
merged_spatial_unfiltered <- readRDS(file.path(data_dir, "data", "seurat", "merged_obj", "spatial_unfiltered.rds"))

# Check distribution of QC metrics
quantile(merged_spatial_unfiltered$nCount_Vizgen, probs=seq(0, 1, 0.1))
mean(merged_spatial_unfiltered$nCount_Vizgen >= 12)

quantile(merged_spatial_unfiltered$nFeature_Vizgen, probs=seq(0, 1, 0.1))
mean(merged_spatial_unfiltered$nFeature_Vizgen >= 10)

quantile(merged_spatial_unfiltered$frac_BlankProbes, probs=seq(0.9, 1, 0.01))
mean(merged_spatial_unfiltered$frac_BlankProbes <= 0.05)

# Violin plots of nCount & nFeature by sample
Idents(merged_spatial_unfiltered) <- "sample"
VlnPlot(merged_spatial_unfiltered, 
        features = c("nCount_Vizgen", "nFeature_Vizgen"), 
        pt.size=0)

# Violin plots of fraction blank probes (for cells passing other QC filters)
merged_spatial_unfiltered[[]] %>%
  filter(nCount_Vizgen >= 12 & nFeature_Vizgen >= 10) %>%
  ggplot(aes(sample, frac_BlankProbes)) +
  geom_violin(aes(fill = sample)) +
  theme_classic() +
  guides(fill="none") +
  labs(y="Fraction of cell count from blank probes")

merged_spatial <- subset(merged_spatial_unfiltered,
                         subset = nCount_Vizgen >= 12 & nFeature_Vizgen >= 10 &
                           frac_BlankProbes <= 0.05)

# Number of cells before and after filtering
bf_cells <- table(merged_spatial_unfiltered$sample)
aft_cells <- table(merged_spatial$sample)
diff_cells <- bf_cells - aft_cells
prop_kept_cells <- round(aft_cells/bf_cells*100, 2)
prop_kept_cells

# Plot the proportion of retained cells by sample
prop_kept_cells <- as.data.frame(prop_kept_cells)
prop_kept_cells$Var1 <- factor(prop_kept_cells$Var1, 
                               levels = c("Ctrl-uninfected", "R3_KO-infected", "C8_R3_KO-infected", "WT-infected"), 
                               ordered = T)
ggplot(prop_kept_cells, aes(Var1, Freq)) +
  geom_col(fill = "royalblue", colour="black") +
  ylab("Percentage of cells passing QC") +
  xlab("") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# Plot the number of cell before and after filtering, by sample
num_kept_cells <- rbind(data.frame(as.data.frame(aft_cells), Group = "Pass QC"), 
                        data.frame(as.data.frame(diff_cells), Group = "Fail QC"))
num_kept_cells$Var1 <- as.character(num_kept_cells$Var1)
num_kept_cells$Var1["Ctrl-uninfected" == num_kept_cells$Var1] <- "WT-uninfected"
num_kept_cells$Var1 <- factor(num_kept_cells$Var1, 
                               levels = c("WT-uninfected", "R3_KO-infected", "C8_R3_KO-infected", "WT-infected"), 
                               ordered = T)

ggplot(num_kept_cells, aes(Var1, Freq, fill = Group)) +
  geom_col(colour="black") +
  ylab("Number of cells") +
  xlab("") +
  theme_classic(base_size = 14) +
  theme(legend.title=element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# Total cells
ncol(merged_spatial_unfiltered) # Unfiltered
ncol(merged_spatial)            # Passing QC
ncol(merged_spatial)/ncol(merged_spatial_unfiltered)*100 # Overall % cells passing QC

# Unfiltered
median(merged_spatial_unfiltered$nCount_Vizgen)
median(merged_spatial_unfiltered$nFeature_Vizgen)

# Cells passing QC
median(merged_spatial$nCount_Vizgen)
median(merged_spatial$nFeature_Vizgen)

# Save filtered data
saveRDS(merged_spatial, file.path(data_dir, "data", "seurat", "merged_obj", "spatial_filtered.rds"))
