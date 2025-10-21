#! usr/bin/Rscript

library(ggplot2)
library(patchwork)
library(reshape2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(gprofiler2)
library(ggrepel)



# Load data
dist_cassette_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_5.rds")

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

x <- x[colnames(betaAdj), ]

# Epitype annotations
my_annotations <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")

# Loading enhancer data
enhancer_general_annotation <- read.table("/Volumes/Data/GeneHancer/GeneHancer_AnnotSV_elements_v5.24.txt", header = T)
enhancer_general_annotation$chr <- paste0("chr", enhancer_general_annotation$chr)

enhancer_to_genes <- read.table("/Volumes/Data/GeneHancer/GeneHancer_AnnotSV_gene_association_scores_v5.24.txt", header = T)

# Loading cpg density data
cpg_density_5000 <- readRDS("/Users/isasiain/PhD/Projects/project_3/analysis/cpg_density/cpg_density_5000_bp.rds")
cpg_density_10000 <- readRDS("/Users/isasiain/PhD/Projects/project_3/analysis/cpg_density/cpg_density_10000_bp.rds")

# FPKM counts
fpkm_data <- read.table("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")

# Generate obkects for genes linked to CpGs
genes <- sapply(annoObj$nameGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})


names(genes) <- annoObj$illuminaID


#
# OVERLAP OF CPGS WITH ENHANCER REGIONS
#

# DISTAL 1
cpgs_1 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 1]

matches_to_enhancers_1 <- sapply(cpgs_1, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_1 != "") / length(matches_to_enhancers_1)


# DISTAL 2
cpgs_2 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 2]

matches_to_enhancers_2 <- sapply(cpgs_2, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_2 != "") / length(matches_to_enhancers_2)



# DISTAL 3
cpgs_3 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 3]

matches_to_enhancers_3 <- sapply(cpgs_3, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_3 != "") / length(matches_to_enhancers_3)


# DISTAL 4
cpgs_4 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 4]

matches_to_enhancers_4 <- sapply(cpgs_4, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_4 != "") / length(matches_to_enhancers_4)


# DISTAL 5
cpgs_5 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 5]

matches_to_enhancers_5 <- sapply(cpgs_5, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_5 != "") / length(matches_to_enhancers_5)

# DISTAL 6
cpgs_6 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 6]

matches_to_enhancers_6 <- sapply(cpgs_6, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_6 != "") / length(matches_to_enhancers_6)


# DISTAL 7
cpgs_7 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 7]

matches_to_enhancers_7 <- sapply(cpgs_7, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_7 != "") / length(matches_to_enhancers_7)


# Calculate mapping percentages
percentages <- c(
  sum(matches_to_enhancers_1 != "") / length(matches_to_enhancers_1),
  sum(matches_to_enhancers_2 != "") / length(matches_to_enhancers_2),
  sum(matches_to_enhancers_3 != "") / length(matches_to_enhancers_3),
  sum(matches_to_enhancers_4 != "") / length(matches_to_enhancers_4),
  sum(matches_to_enhancers_5 != "") / length(matches_to_enhancers_5),
  sum(matches_to_enhancers_6 != "") / length(matches_to_enhancers_6),
  sum(matches_to_enhancers_7 != "") / length(matches_to_enhancers_7)
) * 100  

# Create a data frame for plotting
df_plot <- data.frame(
  Cassette = paste0(1:7),
  PercentMapped = percentages
)

# Plot
p_enhancers <- ggplot(df_plot, aes(x = Cassette, y = PercentMapped)) +
  geom_bar(stat = "identity", fill = "black") +
  ylab("CpGs overlapping\nwith enhancers (%)") +
  xlab("Distal Cassettes") +
  theme_bw(base_size = 14) +  
  theme(axis.text = element_text(size = 12), 
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())


# 
#  REPETITIVE SEQUENCES PER CASSETTE
#


# Calculate mapping percentages of repeatsöä
percentages_repeats <- c(
  sum(annoObj[cpgs_1, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_1, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_2, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_2, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_3, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_3, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_4, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_4, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_5, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_5, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_6, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_6, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_7, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_7, "hasAnyRepeatOverlap"])
) * 100  

# Create a data frame for plotting
df_plot_repeats <- data.frame(
  Cassette = paste0(1:7),
  PercentMapped = percentages_repeats
)

# Plot
p_repetitive <- ggplot(df_plot_repeats, aes(x = Cassette, y = PercentMapped)) +
  geom_bar(stat = "identity", fill = "black") +
  ylab("CpGs overlapping with \nrepetitive sequences (%)") +
  xlab("Distal Cassettes") +
  theme_bw(base_size = 14)

p_enhancers / p_repetitive


#
# CpG DENSITY PER CASSETTE
#

# Get CpGs to analyse
cpgs_for_density <- names(dist_cassette_10$colors)[dist_cassette_10$colors %in% c(7,6,5,4,3,2,1)]

# Filter density values of interest
rownames(cpg_density_5000) <- cpg_density_5000$cpg_id
density_of_interest <-cpg_density_5000[cpgs_for_density,]

# Adding column for cassette
density_of_interest$cassette <- dist_cassette_10$colors[density_of_interest$cpg_id]
density_of_interest$cassette <- factor(density_of_interest$cassette, levels=c(7,6,5,4,3,2,1))

ggplot(density_of_interest, aes(x = cassette, y = cpg_density)) +
  geom_boxplot(outlier.shape = NA, fill = "grey", color = "black") +
  geom_jitter(width = 0.1, size = 0.001, alpha = 0.5, color = "black") +
  theme_bw(base_size = 14) +
  ylim(0,0.005) +
  labs(x = "Cassette", y = "CpG Density\n(CpGs / bp)") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) + coord_flip()


#
# TF enrichment
#

# Proportions of TFs overlapping

tf_proportions <- rbind(
  "1" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 1],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 1]),
  "2" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 2],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 2]),
  "3" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 3],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 3]),
  "4" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 4],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 4]),
  "5" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 5],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 5]),
  "6" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 6],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 6]),
  "7" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 7],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 7])

)

# Compute max absolute deviation from "All"

# Function to compute max pairwise absolute difference for each TF
max_pairwise_diff <- apply(tf_proportions, 2, function(x) {
  max(abs(outer(x, x, "-")))
})

# Select top 50 most differential TFs based on pairwise difference
top_tfs_pw <- names(sort(max_pairwise_diff, decreasing = TRUE)[1:25])

# Subset and melt the original proportions (excluding "All")
plot_data <- as.data.frame(tf_proportions[1:7, top_tfs_pw])
plot_data$Cassette <- factor(rownames(plot_data))
plot_data_long <- melt(plot_data, id.vars = "Cassette", variable.name = "TF", value.name = "Proportion")

ggplot(plot_data_long, aes(x = TF, y = Cassette)) +
  geom_point(aes(size = Proportion, color = Proportion)) +
  scale_color_gradient(low = "lightgrey", high = "darkblue") +
  scale_size_continuous(range = c(1, 6)) +
  theme_bw(base_size = 14) +
    labs(x = "Transcription Factor", y = "Cassette") +
  coord_flip()
