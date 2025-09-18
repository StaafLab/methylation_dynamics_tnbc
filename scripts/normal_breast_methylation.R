#! usr/bin/Rscript

library(ComplexHeatmap)
library(circlize)
library(tidyverse)

#
# LOADING DATA
#

# Loading promoter cassettes
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

x <- x[colnames(betaAdj), ]


# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Loading data (GSE69914)
my_files <- list.files("/Volumes/Data/Project_3/normal_breast_methylation/GSE69914_RAW/", pattern = "Raw.txt$")

# Filtering files. Get only normal tissue
files_to_keep <- paste0("XX0XX0XXXXXXXXXXXX0XXXXXX0XXXXX0XXXXXXXXXXXXX0XXXX",
        "XX0XXX0XXXXXXXXXXX0XXX00XXXXXXXX0XXXXXXXXXXXXXX0XX",
        "XXXXX0XXXXXXXXXXXXXXXX00XXXXXXXXXXXXX0XXXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXX0XXXXX0XXXXXXXXX0XXXXXXXXXXXXX00",
        "XXXXXXX0XXX0XXXX0XXXXXXXXXXXXXX0XXXXXXXXXX0XXXXX0X",
        "XXX0XXXX0XXXXXXX0XXXXXXXXXXXXXX00XXXXXXXXXX00XX0XX",
        "XXXXXXXXXXX0XXX00XXXXXXXXX0XXXXXX0XXXXXXXXXXXX0XXX",
        "XXXX0X0XXXXXXXXX0XXXXXXXXXXX0XXXXXXXXX0XXX0XXXXXX0",
        "XXX0XXX")
files_to_keep <- strsplit(files_to_keep, split="")[[1]]

# samples selection within a series
indices_files <- which(files_to_keep != "X")  # eliminate samples marked as "X
my_files <- my_files[indices_files]

# Readong files and generating matrix

# Define colnames based on a matrix
example <- read.table(paste0("/Volumes/Data/Project_3/normal_breast_methylation/GSE69914_RAW/", "GSM1712369_BCFD3_Raw.txt"))

normal_tissue_methylation <- data.frame(row.names = unname(example$V1))

for(file in my_files) {
  
  my_file <- read.table(paste0("/Volumes/Data/Project_3/normal_breast_methylation/GSE69914_RAW/", file))
  column_name <- strsplit(file, split = "_")[[1]][1]
  
  normal_tissue_methylation[unname(my_file$V1), column_name] <- as.numeric(unname(my_file$V2))             
}


# Loading data (GSE67919)

# Load data (dataframe is called beta and annottaions)
load("/Volumes/Data/Project_3/normal_breast_methylation/GSE67919/GSE67919_Beta.RData")
load("/Volumes/Data/Project_3/normal_breast_methylation/GSE67919/GSE67919_Annotations.RData")

normal_tissue_methylation_96 <- beta

#
# PLOTTING
#


# Prepare long-format data
selected_genes <- c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L")

plot_data <- lapply(selected_genes, function(gene) {
  cpgs <- names(genes)[genes == gene]
  cpgs <- cpgs[cpgs %in% rownames(normal_tissue_methylation_96)]
  data <- normal_tissue_methylation_96[cpgs, , drop = FALSE]
  df <- as.data.frame(t(data))
  df$Sample <- rownames(df)
  df_long <- pivot_longer(df, -Sample, names_to = "CpG", values_to = "Beta")
  df_long$Gene <- gene
  df_long
}) %>% bind_rows()

# Removing cpgs that were not in the detected cassette.
plot_data <- plot_data[plot_data$CpG %in% names(promoter_10$colors)[promoter_10$colors == 10],]

# Plot
ggplot(plot_data, aes(x = CpG, y = Beta)) +
  geom_violin(fill="black") +
  geom_boxplot(fill = "grey90", col= "grey50", outlier.size = 0.5, width=0.2) +
  facet_grid(.~ Gene,scales = "free", space='free') +
  theme_bw(base_size = 12) +
  ylim(0,1) +
  ylab("Beta Value") +
  xlab("CpG Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
