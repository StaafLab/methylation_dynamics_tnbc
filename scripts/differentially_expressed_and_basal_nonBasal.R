#! usr/bin/Rscript

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
library(energy)
library(mclust)
library(reshape2)


#
# LOADING DATA 
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Create a new grouped feature class
annoObj$CpG_context <- feature_class_grouped <- dplyr::case_when(
  annoObj$featureClass %in% c("distal", "distal body") ~ "Distal",
  annoObj$featureClass %in% c("promoter") ~ "Promoter",
  annoObj$featureClass %in% c("proximal dn", "proximal up") ~ "Proximal",
  TRUE ~ as.character(annoObj$featureClass) 
)

# PROMOTER
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")

summary_prom10 <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_beta_10.csv")
rownames(summary_prom10) <- as.character(summary_prom10$Cassette)
summary_prom10$Cassette <- NULL

# DISTAL
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")

summary_dis10 <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes/summary_beta_10.csv")
rownames(summary_dis10) <- as.character(summary_dis10$Cassette)
summary_dis10$Cassette <- NULL

# PROXIMAL
proximal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_10.rds")

summary_prox10 <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/proximal_cassettes/summary_cassettes/summary_beta_10.csv")
rownames(summary_prox10) <- as.character(summary_prox10$Cassette)
summary_prox10$Cassette <- NULL

# FPKM counts
fpkm_data <- read.table("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Epitype annotations
my_annotations <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")


#
# ASSOCIATION WITH BASAL/NON-BASAL SPLIT. Rand index
#

# Create dataframe to store data
basal_nonbasal_clustering_df <- as.data.frame(matrix(nrow=7, ncol=3))

colnames(basal_nonbasal_clustering_df) <- c("Promoter", "Proximal", "Distal")
rownames(basal_nonbasal_clustering_df) <- c(1,2,3,4,5,6,7)

# Filling the dataframe
for (cassette in c(1,2,3,4,5,6,7)) {
  
  # PROMOTER
  
  # Extract adjusted cpgs per cassette
  betas <- t(betaAdj[names(promoter_10$colors)[promoter_10$colors == cassette], x$PD_ID])
  
  # Compute cludters through hierarchical clustering
  dist_mat <- dist(betas, method = "euclidean")
  hc <- hclust(dist_mat, method = "ward.D2")
  
  # Cut tree into 2 clusters
  clusters <- cutree(hc, k = 2)
  
  # Calculate Adjusted Rand Index between clustering and true basal groups
  basal_nonbasal_clustering_df[cassette, "Promoter"] <- adjustedRandIndex(clusters, as.factor(x$PAM50_Basal_NCN))
  
  # PROXIMAL
  
  # Extract adjusted cpgs per cassette
  betas <- t(betaAdj[names(proximal_10$colors)[proximal_10$colors == cassette], x$PD_ID])
  
  # Compute cludters through hierarchical clustering
  dist_mat <- dist(betas, method = "euclidean")
  hc <- hclust(dist_mat, method = "ward.D2")
  
  # Cut tree into 2 clusters
  clusters <- cutree(hc, k = 2)
  
  # Calculate Adjusted Rand Index between clustering and true basal groups
  basal_nonbasal_clustering_df[cassette, "Proximal"] <- adjustedRandIndex(clusters, as.factor(x$PAM50_Basal_NCN))
  
  # DISTAL
  
  # Extract adjusted cpgs per cassette
  betas <- t(betaAdj[names(distal_10$colors)[distal_10$colors == cassette], x$PD_ID])
  
  # Compute cludters through hierarchical clustering
  dist_mat <- dist(betas, method = "euclidean")
  hc <- hclust(dist_mat, method = "ward.D2")
  
  # Cut tree into 2 clusters
  clusters <- cutree(hc, k = 2)
  
  # Calculate Adjusted Rand Index between clustering and true basal groups
  basal_nonbasal_clustering_df[cassette, "Distal"] <- adjustedRandIndex(clusters, as.factor(x$PAM50_Basal_NCN))
  
}




# Convert matrix to data frame if needed
mat <- as.matrix(basal_nonbasal_clustering_df)
long_df <- melt(mat)

# Rename for clarity
colnames(long_df) <- c("Row", "Column", "ARI")

ggplot(long_df, aes(x = Column, y = factor(Row, levels = c(7,6,5,4,3,2,1)))) +
  geom_point(aes(size = ARI, color = ARI)) +
  scale_size(range = c(1, 10)) +
  scale_color_gradient2(
    low = "blue",      # for negative values
    mid = "grey",     # for zero
    high = "red",      # for positive values (up to 1)
    midpoint = 0,      # center the gradient at zero
    limits = c(-0.2, 1), # optional: fixes the color scale range
    name = "ARI"
  ) +
  ylab("Cassette") +
  xlab("CpG context") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(size = "Adjusted\nRand Index")



#
# CLUSTERING SAMPLES. Basal-NonBasal
#

# Get CpGs of the proximal, distal and promoter cassettes linked to PAM50 Basal/NonBasal
# Using cassettes whose Rand Index is higher than 0.65
my_cpgs_prom <-  c(
  names(promoter_10$colors)[promoter_10$colors == "1"]
)

my_cpgs_dis <-  c(
  names(distal_10$colors)[distal_10$colors == "2"],
  names(distal_10$colors)[distal_10$colors == "3"]
)

my_cpgs_prox <- c(
  names(proximal_10$colors)[proximal_10$colors == "1"],
  names(proximal_10$colors)[proximal_10$colors == "3"]
)

my_cpgs_all <- c(
  names(promoter_10$colors)[promoter_10$colors == "1"],
  names(distal_10$colors)[distal_10$colors == "2"],
  names(distal_10$colors)[distal_10$colors == "3"],
  names(proximal_10$colors)[proximal_10$colors == "1"],
  names(proximal_10$colors)[proximal_10$colors == "3"]
)

# Generate data frame to store groups
groupings_df <- data.frame(matrix(nrow = length(colnames(summary_prox10)), ncol = 5))
rownames(groupings_df) <- colnames(summary_prox10)
colnames(groupings_df) <- c("group_prom", "group_dis", "group_prox", "group_all", "PAM50")        

# CLUSTERING IN TWO GROUPS

# PROMOTER
distance_matrix <- dist(t(betaAdj[my_cpgs_prom,]), method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
groupings_df$group_prom <- cutree(hc, k = 2)

# DISTAL
distance_matrix <- dist(t(betaAdj[my_cpgs_dis,]), method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
groupings_df$group_dis <- cutree(hc, k = 2)

# PROXIMAL
distance_matrix <- dist(t(betaAdj[my_cpgs_prox,]), method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
groupings_df$group_prox <- cutree(hc, k = 2)

# ALL
distance_matrix <- dist(t(betaAdj[my_cpgs_all,]), method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
groupings_df$group_all <- cutree(hc, k = 2)

#PAM50
groupings_df$PAM50 <- ifelse(my_annotations[colnames(betaAdj), "PAM50"] == "Uncl.", 
       "Uncl.", 
       ifelse(my_annotations[colnames(betaAdj), "PAM50"] == "Basal", 
              "Basal", 
              "Non-Basal"))

# Define grouping variables
group_vars <- c("group_prox", "group_prom", "group_dis", "group_all")

# ANALYSISNG DIFFERENCES IN GROUPINGS

# USING PERCENTAGES

# Generate contingency tables and compute percentages
summary_list <- lapply(group_vars, function(var) {
  tbl <- as.data.frame.matrix(table(groupings_df[[var]], groupings_df$PAM50))
  tbl$Grouping_Method <- var
  tbl$Group_Value <- rownames(tbl)
  
  # Convert counts to percentages (row-wise)
  tbl[, 1:(ncol(tbl)-2)] <- tbl[, 1:(ncol(tbl)-2)] / rowSums(tbl[, 1:(ncol(tbl)-2)]) * 100
  
  tbl
})

# Combine all tables into one
summary_table <- do.call(rbind, summary_list)

# Convert to long format for ggplot
summary_long <- summary_table %>%
  pivot_longer(cols = -c(Grouping_Method, Group_Value), 
               names_to = "PAM50", values_to = "Percentage")

# Convert Group_Value to factor for proper ordering
summary_long$Group_Value <- factor(summary_long$Group_Value)

# Plot the data
ggplot(summary_long, aes(x = Group_Value, y = Percentage, fill = PAM50)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  facet_wrap(~Grouping_Method, scales = "free_x") +  # One plot per grouping method
  scale_fill_manual(values = c("Basal" = "indianred1", "Non-Basal" = "blue","Uncl." = "gray")) +  
  labs(x = "Group Value", y = "Percentage", fill = "PAM50") +
  theme_classic()

# USING COUNTS

# Generate contingency tables for each grouping variable against PAM50
summary_list_counts <- lapply(group_vars, function(var) {
  tbl <- as.data.frame.matrix(table(groupings_df[[var]], groupings_df$PAM50))
  tbl$Grouping_Method <- var  # Add column to identify the method
  tbl$Group_Value <- rownames(tbl)  # Store group value
  tbl
})

# Combine all tables into a single data frame
summary_table_counts <- do.call(rbind, summary_list_counts)

# Reorder columns for readability
summary_table_counts <- summary_table_counts[, c("Grouping_Method", "Group_Value", colnames(summary_table_counts)[1:(ncol(summary_table_counts) - 2)])]

# Convert to long format for ggplot
summary_long <- summary_table_counts %>%
  pivot_longer(cols = -c(Grouping_Method, Group_Value), 
               names_to = "PAM50", values_to = "Counts")

# Convert Group_Value to factor for proper ordering
summary_long$Group_Value <- factor(summary_long$Group_Value)

# Plot the data
ggplot(summary_long, aes(x = Group_Value, y = Counts, fill = PAM50)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  facet_wrap(~Grouping_Method, scales = "free_x") +  # One plot per grouping method
  scale_fill_manual(values = c("Basal" = "indianred1", "Non-Basal" = "blue","Uncl." = "gray")) +  
  labs(x = "Group Value", y = "Number of samples", fill = "PAM50") +
  theme_classic()


# Combine all tables into a single data frame
summary_table <- do.call(rbind, summary_list)

# Reorder columns for readability
summary_table <- summary_table[, c("Grouping_Method", "Group_Value", colnames(summary_table)[1:(ncol(summary_table) - 2)])]

# Print the final table
print(summary_table)


#
# PLOTTING DIFFERENTIALLY EXPRESSED GENES
#

# PROMOTER

# Clustering

# Getting data
my_cpgs_prom <-  c(
  names(promoter_10$colors)[promoter_10$colors == "1"]
)

# Clustering
distance_matrix <- dist(t(betaAdj[my_cpgs_prom,]), method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
clusters_to_analyse <- cutree(hc, k = 2)


# Extracting the subset data
fpkm_subset_1 <- fpkm_data[, names(clusters_to_analyse)[clusters_to_analyse== 1]]
fpkm_subset_2 <- fpkm_data[, names(clusters_to_analyse)[clusters_to_analyse == 2]]

# Identifying genes in cassette 1
cpgs_in_cassette <- c(names(promoter_10$colors[promoter_10$colors == 1]))
genes_in_cassette <- unique(genes[cpgs_in_cassette])

# Count the number of CpGs linked to each gene
cpg_counts <- table(genes[cpgs_in_cassette])

# Initialize results dataframe
results <- data.frame(
  Gene = genes_in_cassette,
  Wilcoxon_p = NA,
  Bonferroni_p = NA,
  Log2_fold_change = NA,
  CpG_Count = NA
)

# Compute statistics
for (gene in genes_in_cassette) {
  if (gene %in% rownames(fpkm_data)) {
    expr_group1 <- as.numeric(fpkm_subset_1[gene, ])
    expr_group2 <- as.numeric(fpkm_subset_2[gene, ])
    
    # Wilcoxon test
    test_result <- wilcox.test(expr_group1, expr_group2, exact = FALSE)
    
    # Compute median difference
    log_fold_change <- calculate_logFC(expr_group1, expr_group2)
    
    # Store results
    results[results$Gene == gene, "Wilcoxon_p"] <- test_result$p.value
    results[results$Gene == gene, "Log2_fold_change"] <- log_fold_change
  }
  
  # Store CpG count for the gene
  results[results$Gene == gene, "CpG_Count"] <- cpg_counts[gene]
}

# Apply Bonferroni correction
results$Bonferroni_p <- p.adjust(results$Wilcoxon_p, method = "bonferroni")

# Sort results by Bonferroni p-value
results <- results[order(results$Bonferroni_p), ]


# VOLCANO PLOT

# Compute -log10(Bonferroni p-value)
results$neg_log10_p <- -log10(results$Bonferroni_p)

# Define color based on significance
results$color <- ifelse(results$Bonferroni_p > 0.05, "grey", "blue")

# Volcano plot
ggplot(results, aes(x = Log2_fold_change, y = neg_log10_p, size = CpG_Count, color = color)) +
  geom_point(alpha = 0.7) + 
  geom_text_repel(
    data = subset(results, Bonferroni_p <= 0.00000000001),
    aes(label = Gene),
    size = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("#1f78b4", "grey")) +  
  scale_size(range = c(1, 6)) +                    
  theme_classic() +
  labs(
    x = "Log2 Fold Change in Expression",
    y = "-log10(Bonferroni p-value)",
    size = "CpG Count"
  ) +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  ) +
  theme_bw()


# PROXIMAL

# Getting data
my_cpgs_prox <-  c(
  names(proximal_10$colors)[proximal_10$colors == "1"]
)

# Clustering
distance_matrix <- dist(t(betaAdj[my_cpgs_prox,]), method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
clusters_to_analyse <- cutree(hc, k = 2)

# Extracting the subset data
fpkm_subset_1 <- fpkm_data[, names(clusters_to_analyse)[clusters_to_analyse== 1]]
fpkm_subset_2 <- fpkm_data[, names(clusters_to_analyse)[clusters_to_analyse == 2]]

# Identifying genes in cassette 1
cpgs_in_cassette <- names(proximal_10$colors[proximal_10$colors %in% c(1)])
genes_in_cassette <- unique(genes[cpgs_in_cassette])

# Count the number of CpGs linked to each gene
cpg_counts <- table(genes[cpgs_in_cassette])

# Initialize results dataframe
results <- data.frame(
  Gene = genes_in_cassette,
  Wilcoxon_p = NA,
  Bonferroni_p = NA,
  Log2_fold_change = NA,
  CpG_Count = NA
)

# Compute statistics
for (gene in genes_in_cassette) {
  if (gene %in% rownames(fpkm_data)) {
    expr_group1 <- as.numeric(fpkm_subset_1[gene, ])
    expr_group2 <- as.numeric(fpkm_subset_2[gene, ])
    
    # Wilcoxon test
    test_result <- wilcox.test(expr_group1, expr_group2, exact = FALSE)
    
    # Compute median difference
    log_fold_change <- calculate_logFC(expr_group1, expr_group2)
    
    # Store results
    results[results$Gene == gene, "Wilcoxon_p"] <- test_result$p.value
    results[results$Gene == gene, "Log2_fold_change"] <- log_fold_change
  }
  
  # Store CpG count for the gene
  results[results$Gene == gene, "CpG_Count"] <- cpg_counts[gene]
}

# Apply Bonferroni correction
results$Bonferroni_p <- p.adjust(results$Wilcoxon_p, method = "bonferroni")

# Sort results by Bonferroni p-value
results <- results[order(results$Bonferroni_p), ]


# VOLCANO PLOT

# Compute -log10(Bonferroni p-value)
results$neg_log10_p <- -log10(results$Bonferroni_p)

# Define color based on significance
results$color <- ifelse(results$Bonferroni_p > 0.05, "grey", "blue")

# Volcano plot
ggplot(results, aes(x = Log2_fold_change, y = neg_log10_p, size = CpG_Count, color = color)) +
  geom_point(alpha = 0.7) + 
  geom_text_repel(
    data = subset(results, Bonferroni_p <= 0.0000000000008),
    aes(label = Gene),
    size = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("#1f78b4", "grey")) +  
  scale_size(range = c(1, 6)) +                    
  theme_classic() +
  labs(
    x = "Log2 Fold Change in Expression",
    y = "-log10(Bonferroni p-value)",
    size = "CpG Count"
  ) +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  ) +
  theme_bw()


#
# PLOTTING SPECIFIC GENES
#

# Annotations
pam50_annotations <- my_annotations[colnames(betaAdj), "PAM50"]
pam50_annotations <- ifelse(pam50_annotations == "Uncl.", 
                            "Uncl.", 
                            ifelse(pam50_annotations == "Basal", 
                                   "Basal", 
                                   "Non-Basal"))
tnbc_annotation <- my_annotations[colnames(betaAdj), "TNBC"]
HRD_annotation <- my_annotations[colnames(betaAdj), "HRD"]
epi_annotation <- my_annotations[colnames(betaAdj), "NMF_atacDistal"]
tils_annotation <- as.numeric(x[colnames(betaAdj), "TILs"])


## LDHB
current_gene_id = "LDHB"


# Create top anotation
top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                    TNBC = tnbc_annotation,
                                    HRD = HRD_annotation,
                                    Epitype = epi_annotation,
                                    TILs = anno_points(tils_annotation,
                                                       ylim=c(0,100),
                                                       size=unit(0.75, "mm"),
                                                       axis_param = list(
                                                         side="left",
                                                         at=c(0,25,50,75,100),
                                                         labels=c("0","25","50","75","100")
                                                       )),
                                    col = list(
                                      "PAM50"=c("Basal"="indianred1", "Non-Basal"="darkblue","Uncl."="grey"),
                                      "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                      "IM"=c("Negative"="grey", "Positive"="black"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                      "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                  "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")
                                    )
)
# Generate bottom annotation
bottom_annotation <- HeatmapAnnotation(
  "FPKM" = anno_barplot(as.numeric(fpkm_data[current_gene_id, colnames(betaAdj)]))
)

# Updated left_annotation with color scale
right_annotation <- rowAnnotation(
  "Normal beta" = rowMeans(betaNorm[names(genes)[genes == current_gene_id],]),
  #"ATAC" = annoObj$hasAtacOverlap[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
  col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")),
             "ATAC" = c("0" = "white", "1"= "black"))
)

# CpG annotation. Context and included in cassette
left_annotation <- rowAnnotation("CpG_in_cassette" = names(genes)[genes == current_gene_id] %in% names(promoter_10$colors)[promoter_10$colors == 1],
                                 "Context"= annoObj$CpG_context[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                                 col=list("Context"=c("Distal" = "#f8766d", 
                                                      "Promoter" = "#00ba38", 
                                                      "Proximal" = "#619cff"),
                                          "CpG_in_cassette"=c("TRUE" = "black", 
                                                              "FALSE" = "white")))

# Cluster based on methylation
cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                 annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])

cluster_promoter <- kmeans(t(betaAdj[names(cpgs)[cpgs=="promoter"],]), centers = 2)

# Determine hypo and hypermethylated cluster
promoter_state <- if (mean(betaAdj[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==1]) >
                      mean(betaAdj[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==2])) {
  
  as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
  
} else {
  
  as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
  
}


# Heatmap of genes
Heatmap(
  betaAdj[names(genes)[genes == current_gene_id],],
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  column_split = promoter_state,
  top_annotation = top_annotation,
  bottom_annotation = bottom_annotation,
  right_annotation = right_annotation,
  left_annotation = left_annotation,
  clustering_distance_columns =  "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_rows =  "euclidean",
  clustering_method_rows = "ward.D2",
  name = "Tumor beta"
)


## SPINK 8
current_gene_id = "SPINK8"


# Create top anotation
top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                    TNBC = tnbc_annotation,
                                    HRD = HRD_annotation,
                                    Epitype = epi_annotation,
                                    TILs = anno_points(tils_annotation,
                                                       ylim=c(0,100),
                                                       size=unit(0.75, "mm"),
                                                       axis_param = list(
                                                         side="left",
                                                         at=c(0,25,50,75,100),
                                                         labels=c("0","25","50","75","100")
                                                       )),
                                    col = list(
                                      "PAM50"=c("Basal"="indianred1", "Non-Basal"="darkblue","Uncl."="grey"),
                                      "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                      "IM"=c("Negative"="grey", "Positive"="black"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                      "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                  "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")
                                    )
)
# Generate bottom annotation
bottom_annotation <- HeatmapAnnotation(
  "FPKM" = anno_barplot(as.numeric(fpkm_data[current_gene_id, colnames(betaAdj)]))
)

# Updated left_annotation with color scale
right_annotation <- rowAnnotation(
  "Normal beta" = rowMeans(betaNorm[names(genes)[genes == current_gene_id],]),
  #"ATAC" = annoObj$hasAtacOverlap[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
  col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")),
             "ATAC" = c("0" = "white", "1"= "black"))
)

# CpG annotation. Context and included in cassette
left_annotation <- rowAnnotation("CpG_in_cassette" = names(genes)[genes == current_gene_id] %in% names(proximal_10$colors)[proximal_10$colors == 1],
                                 "Context"= annoObj$CpG_context[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                                 col=list("Context"=c("Distal" = "#f8766d", 
                                                      "Promoter" = "#00ba38", 
                                                      "Proximal" = "#619cff"),
                                          "CpG_in_cassette"=c("TRUE" = "black", 
                                                              "FALSE" = "white")))

# Cluster based on methylation
cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                 annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])

cluster_proximal <- kmeans(betaAdj[names(cpgs)[cpgs=="proximal up" | cpgs=="proximal dn"],], centers = 2)

# Determine hypo and hypermethylated cluster
proximal_state <- if (mean(betaAdj[names(cpgs)[cpgs=="proximal up" | cpgs=="proximal dn"],cluster_proximal$cluster==1]) >
                      mean(betaAdj[names(cpgs)[cpgs=="proximal up" | cpgs=="proximal dn"],cluster_proximal$cluster==2])) {
  
  as.factor(ifelse(cluster_proximal$cluster == 1, "Hypermethylated", "Hypomethylated"))
  
} else {
  
  as.factor(ifelse(cluster_proximal$cluster == 2, "Hypermethylated", "Hypomethylated"))
  
}


# Heatmap of genes
Heatmap(
  betaAdj[names(genes)[genes == current_gene_id],],
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  column_split = proximal_state,
  top_annotation = top_annotation,
  bottom_annotation = bottom_annotation,
  right_annotation = right_annotation,
  left_annotation = left_annotation,
  clustering_distance_columns =  "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_rows =  "euclidean",
  clustering_method_rows = "ward.D2",
  name = "Tumor beta"
)

