#! usr/bin/Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

#
# LOADING DATA
#

# CosMs data
oas2_per_cell_12 <-  read.csv("/Users/isasiain/PhD/Projects/project_3/data/block12_OAS2_in_tumour.csv")

# Pdid to cores mapping
core_to_pdid <- read.table("/Volumes/Data/CosMx/mapping.txt", sep = "\t", header = T)
core_to_pdid <- unique(core_to_pdid[,c("tmaID", "PDid")])
rownames(core_to_pdid) <- core_to_pdid$tmaID

# Loading betas
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

#
# PLOT MEAN EXPRESSION OF ALL GENES (OAS2 + CONTROLS)
#

# Genes to plot
my_genes <- c("OAS2_Count", "BRCA1_Count", "SMO_Count", "GSTP1_Count" , "AR_Count")

# Reshape to long format
df_long <- oas2_per_cell_12 %>%
  filter(GMM_Label == "Tumor") %>%
  pivot_longer(cols = all_of(my_genes), names_to = "Gene", values_to = "Expression")

# Calculate mean expression per gene
mean_expr <- df_long %>%
  group_by(Gene) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE))

# Plot
ggplot(mean_expr, aes(x = Gene, y = Mean_Expression, fill = Gene)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Mean Expression per Tumour Cell",
       y = "Mean Expression", x = "Gene") +
  theme_minimal() +
  theme(legend.position = "none")


#
# DEFINING METHYLATION STATE
#

current_gene_id = "OAS2"

# Cluster based on methylation
cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                 annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])

cluster_promoter <- kmeans(t(betaAdj[names(cpgs),]), centers = 2)

# Determine hypo and hypermethylated cluster
promoter_state <- if (mean(betaAdj[names(cpgs),cluster_promoter$cluster==1]) >
                           mean(betaAdj[names(cpgs),cluster_promoter$cluster==2])) {
  
  as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
  
} else {
  
  as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
  
}


#
# FILTERING AND SUMMARIZING TUMOR CELLS
#

### BLOCK 1 and 2

par(mfrow=c(1,1))

# Plotting panCK intensity vs cell type. Group 1 corresponds to tumour
boxplot(oas2_per_cell_12$Mean.PanCK ~ oas2_per_cell_12$GMM_Label)

# Remove non-tumour cells from df
oas2_per_cell_12 <- oas2_per_cell_12[oas2_per_cell_12$GMM_Label == "Tumour", ]

# Adding pdid
oas2_per_cell_12$"PDid" <- sapply(
  oas2_per_cell_12$Tissue_ID, function(tissue_name) {
    core_to_pdid[strsplit(tissue_name, "-")[[1]][1], "PDid"]
  }
)

# Methylation state
oas2_per_cell_12$"Methylation" <- sapply(
  oas2_per_cell_12$PDid, function(pdid) {
    promoter_state[pdid]
  }
)


# SUMMARIZING

# Summarizing per Core (Tissue ID)
summary_per_tissue_12 <- oas2_per_cell_12 %>%
  group_by(PDid) %>%
  summarise(
    PDid = dplyr::first(PDid),
    Methylation = dplyr::first(Methylation),
    
    mean_OAS2 = mean(OAS2_Count, na.rm = TRUE),
    prop_OAS2 = mean(OAS2_Count > 0, na.rm = TRUE),
    
    #mean_BRCA1 = mean(BRCA1_Count, na.rm = TRUE),
    #prop_BRCA1 = mean(BRCA1_Count > 0, na.rm = TRUE),
    
    #mean_SMO = mean(SMO_Count, na.rm = TRUE),
    #prop_SMO = mean(SMO_Count > 0, na.rm = TRUE),
    
    #mean_GSTP1 = mean(GSTP1_Count, na.rm = TRUE),
    #prop_GSTP1 = mean(GSTP1_Count > 0, na.rm = TRUE),
    
    #mean_AR = mean(AR_Count, na.rm = TRUE),
    #prop_AR = mean(AR_Count > 0, na.rm = TRUE)
  )


#
# ADDING COHORT COMPOSITION INFORMATION
#

# Distinct pdIDS
summary_per_tissue_12_distinct_pdid <- summary_per_tissue_12 %>%
  distinct(PDid, .keep_all = TRUE)

# Adding annotations
summary_per_tissue_12$PAM50 <- x[summary_per_tissue_12$PDid, "PAM50_Basal_NCN"]
summary_per_tissue_12$Lehmann_4 <- x[summary_per_tissue_12$PDid, "TNBCtype4_n235_notPreCentered"]
summary_per_tissue_12$IM <- x[summary_per_tissue_12$PDid, "TNBCtype_IMpositive"]
summary_per_tissue_12$LN <- x[summary_per_tissue_12$PDid, "LNbinary"]
summary_per_tissue_12$grade <- x[summary_per_tissue_12$PDid, "Grade"]
summary_per_tissue_12$hrd <- x[summary_per_tissue_12$PDid, "HRD.2.status"]


#
# PLOTTING
#

## BLOCK 1 and 2

# Mean OAS2 expression

# Count number of data points per class
counts_cores <- table(summary_per_tissue_12$Methylation)
counts_samples <- c("Hypermethylated" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$Methylation == "Hypermethylated","PDid"]))),
                    "Hypomethylated" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$Methylation == "Hypomethylated","PDid"]))))


labels_with_counts <- paste0(
  "n=", counts_cores[names(counts_cores)]
)

# Prepare the labels as a dataframe
label_df <- data.frame(
  Methylation = c("Hypermethylated", "Hypomethylated"),
  label = labels_with_counts)


ex_plot <- ggplot(filter(summary_per_tissue_12, !is.na(summary_per_tissue_12$Methylation)), aes(x=Methylation, y=mean_OAS2)) +
  geom_violin(fill="black") +
  geom_boxplot(width=0.1, size=0.2) +
  theme_bw(base_size = 14) +
  xlab(NULL) +
  ylim(0,1.52) +
  ylab("Mean OAS2 expression in tumor cells") + 
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     comparisons = list(c("Hypomethylated", "Hypermethylated")), 
                     label.x = c("Hypomethylated", "Hypermethylated"),
                     label.y = 1.4, size = 5) +
  geom_text(data = label_df, aes(x = Methylation, y = 1.2, label = label), 
            inherit.aes = FALSE, size = 5, vjust = 0)

# Proportion of cells with detected expressin of OAS2
prop_plot <- ggplot(filter(summary_per_tissue_12, !is.na(summary_per_tissue_12$Methylation)), aes(x=Methylation, y=prop_OAS2 * 100)) +
  geom_violin(fill="black") +
  geom_boxplot(width=0.1, size=0.2) +
  theme_bw(base_size = 14) +
  xlab(NULL) +
  ylim(0,38) +
  ylab("% of tumor cells expressing OAS2") + 
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     comparisons = list(c("Hypomethylated", "Hypermethylated")), 
                     label.x = c("Hypomethylated", "Hypermethylated"),
                     label.y = 35, size = 5) +
  geom_text(data = label_df, aes(x = Methylation, y = 30, label = label), 
            inherit.aes = TRUE, size = 5, vjust = 0)

ex_plot | prop_plot


##geom_violin()### CONTROL. Plotting AR vs Basal/NonBasal

# Mean AR expression

# Count number of data points per class
counts_cores <- table(summary_per_tissue_12$PAM50)
counts_samples <- c("Basal" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$PAM50 == "Basal","PDid"]))),
                    "nonBasal" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$PAM50 == "nonBasal","PDid"]))))


labels_with_counts <- paste0(
  "n=", counts_cores[names(counts_cores)], 
  "\ns=", counts_samples[names(counts_cores)]
)


# Draw boxplot with custom x-axis labels
boxplot(summary_per_tissue_12$mean_AR ~ summary_per_tissue_12$PAM50,
        ylab = "Mean Expression in Tumour cells",
        xlab = "Promoter Methylation",,
        ylim = c(0, 2),
        frame = FALSE)

# Add the annotation
text(x = 1:2, y = 1.7, labels = labels_with_counts, xpd = TRUE, cex = 0.8)


# Add jittered points
stripchart(summary_per_tissue_12$mean_AR ~ summary_per_tissue_12$PAM50,
           method = "jitter", 
           pch = 16,
           cex = 0.6,
           col = rgb(0, 0, 0, 0.5),
           vertical = TRUE,
           add = TRUE)

# Perform Wilcoxon test
wilcox_res <- wilcox.test(summary_per_tissue_12$mean_AR ~ summary_per_tissue_12$PAM50)

# Add p-value to plot
pval <- wilcox_res$p.value
text(x = 1.3, 
     y = 1.2, 
     labels = paste0("Wilcoxon's p = ", signif(pval, 3)),
     pos = 3, cex = 0.9)


#
# SAVING DATA
#

write.csv(oas2_per_cell_12[c("Tissue_ID", "GMM_Label", "Mean.PanCK", "OAS2_Count", "batch", "PDid", "Methylation")],
          "/Users/isasiain/PhD/Projects/project_3/data/supp_data/cosmx_supp_data.csv")



