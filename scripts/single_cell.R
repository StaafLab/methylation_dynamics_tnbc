#! usr/bin/Rscript

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

#
# Loading data
#

only_tum <- readRDS("/Volumes/Data/Project_3/single_cell_brca/gse161529/GSE161529/seurat_objects/SeuratObject_TNBCTum.rds")
only_tum <- UpdateSeuratObject(only_tum)

only_str <- readRDS("/Volumes/Data/Project_3/single_cell_brca/gse161529/GSE161529/seurat_objects/SeuratObject_TNBCSub.rds")
only_str <- UpdateSeuratObject(only_str)

only_tc <- readRDS("/Volumes/Data/Project_3/single_cell_brca/gse161529/GSE161529/seurat_objects/SeuratObject_TNBCTC.rds")
only_tc <- UpdateSeuratObject(only_tc)


# GENES OF INTEREST IN CANCER CELLS

# Count number of cells per tumor
table(only_tum$group)
table(only_str$group)
table(only_tc$group)


# Count tumor cells per group
cell_counts <- only_tum@meta.data %>%
  group_by(group) %>%
  summarise(n_cells = n())

# Base DotPlot
gene_expression <- DotPlot(
  only_tum,
  features = c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L"),
  group.by = "group",
  assay = "RNA"
) +
  RotatedAxis() +
  xlab("Genes") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) +
  scale_color_gradient2(
    limits = c(-1, 4),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  )

# STROMAL AND IMMUNE CELLS PER SAMPLE
cell_types <- c("0"="T cells", "1"="Macrophagues", "2"="Plasma cells", "3"="Fibroblasts", "4"="T cells","5"="B cells", "6"="Dendritic cells", "7"="Endothelial cells", "8"="Pericytes", "9"="Myeloid cells")

# Get normalized cell counts
cell_type_counts <- table(only_str$group, unname(cell_types[only_str$seurat_clusters]))
tumor_counts <- as.numeric(table(only_tum$group)[rownames(cell_type_counts)])
cell_type_counts <- cbind(cell_type_counts, Tumor = tumor_counts)

total_cells_per_sample <- rowSums(cell_type_counts)
cell_type_fractions <- sweep(cell_type_counts, 1, total_cells_per_sample, FUN = "/")

# Convert to long format for plotting
cell_type_counts_long <- melt(t(cell_type_fractions))

# Add original count values to the long format
cell_type_counts_long$original_counts <- mapply(function(x, y) cell_type_counts[x, y],
                                                cell_type_counts_long$Var2, 
                                                cell_type_counts_long$Var1)

# Plotting with original counts for tile annotation
composition_of_microenvironment <- ggplot(cell_type_counts_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Midpoint set to 0
  labs(x = "Cell types", y = "Sample", fill = "Proportion of cell type") +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_text(aes(label = original_counts),  color = "black", size = 4)  +  # Use raw counts for annotation
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis


# Combine the plots
(gene_expression | composition_of_microenvironment) + 
  plot_layout(widths = c(1, 1.75), guides = "collect") &
  theme(legend.position = "right")


# T CELL SUBTYPES

cell_types <- c("0"="Effector T cells", "1"="Naive/resting T cells", "2"="Regulatory T cells", "3"="Cycling T cells", "4"="TR cells (memory)","5"="NK cells", "6"="Plasma")


cell_type_counts <- t(prop.table(table(only_tc$group, cell_types[only_tc$seurat_clusters]), margin=1))

# Convert to long format for plotting
cell_type_counts_long <- melt(cell_type_counts)

# Plotting
ggplot(cell_type_counts_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Midpoint set to 0
  labs(x = "Group", y = "Cell Type", fill = "Scaled cell counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels if needed
  geom_text(aes(label = round(value, 2)), color = "black", size = 3)  # Add numbers inside the tiles


DotPlot(only_tc, features = c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "MKI67"), group.by = "group", assay = "RNA")
