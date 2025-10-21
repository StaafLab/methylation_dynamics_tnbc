#! usr/bin/Rscript

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(biomaRt)
library(patchwork)
library(ggplot2)



#
# GETTING ENSEMBL IDS FOR GENES OF INTEREST
#

# Defining genes of interest
genes_of_interest <- c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L")


# Getting ensembl ids
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_of_interest_ensembl <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",  
  values = genes_of_interest, 
  mart = ensembl  
)

# Create dictionary of gene ids and ensembl ids
genes_dict <- setNames(genes_of_interest_ensembl$hgnc_symbol, 
                       genes_of_interest_ensembl$ensembl_gene_id)

rm(ensembl)

#
# ANALYSING SINGLE CELL DATA
#

# EPITHELIAL CELLS

#The data comes from: 10.1038/s41588-024-01688-9 
hbca_epithelial <- readRDS("/Volumes/Data/Project_3/normal_breast_single_cell/hbca_epithelial.rds")


# Getting Ensembl IDS from gene names
epi_normal <- DotPlot(hbca_epithelial, 
                      features = names(genes_dict), 
                      assay = "RNA",
                      group.by = "cell_type") + 
  theme_bw(base_size = 14) +
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(labels = NULL) +
  scale_y_discrete(labels = function(y) {
    sapply(y, function(label) {
      paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
    })
  }) +
  scale_color_gradient2(
    limits = c(-2, 3),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  ) +
  scale_size(
    limits = c(0, 40),
    range = c(0, 8)   # Adjust dot size appearance as desired
  )

saveRDS(epi_normal, "PhD/Projects/project_3/analysis/dotplot_epi.rds")
rm(hbca_epithelial)
gc()

# IMMUNE CELLS

#The data comes from: 10.1038/s41588-024-01688-9 
hbca_immune <- readRDS("/Volumes/Data/Project_3/normal_breast_single_cell/hbca_immune.rds")


# Getting Ensembl IDS from gene names
immune_normal <- DotPlot(hbca_immune, 
        features = names(genes_dict), 
        assay = "RNA",
        group.by = "cell_type") + 
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = genes_dict) +
  scale_y_discrete(labels = function(y) {
    sapply(y, function(label) {
      paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
    })
  }) +
  scale_color_gradient2(
    limits = c(-2, 3),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  ) +
  scale_size(
    limits = c(0, 40),
    range = c(0, 8)   # Adjust dot size appearance as desired
  )


saveRDS(immune_normal, "PhD/Projects/project_3/analysis/dotplot_immune.rds")
rm(hbca_immune)
gc()

# STROMAL CELLS

#The data comes from: 10.1038/s41588-024-01688-9 
hbca_stroma <- readRDS("/Volumes/Data/Project_3/normal_breast_single_cell/hbca_stroma.rds")

# Getting Ensembl IDS from gene names
stroma_normal <- DotPlot(hbca_stroma, 
                         features = names(genes_dict), 
                         assay = "RNA",
                         group.by = "cell_type") + 
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = genes_dict) +
  scale_y_discrete(labels = function(y) {
    sapply(y, function(label) {
      paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
    })
  }) +
  scale_color_gradient2(
    limits = c(-2, 3),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  ) +
  scale_size(
    limits = c(0, 40),
    range = c(0, 8)   # Adjust dot size appearance as desired
  )

saveRDS(stroma_normal, "PhD/Projects/project_3/analysis/dotplot_stroma.rds")
rm(hbca_stroma)
gc()

# Epi:Immune = 2:1 in vertical stacking
epi_normal / immune_normal / stroma_normal + plot_layout(heights = c(5, 15), guides = "collect")

