#! usr/bin/Rscript

library(Boruta)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(ggrepel)
library(energy)


#
# LOADING DATA 
#

set.seed(123)

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Defining gene-cpg dictionary
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Load cassettes
prom_cassettes <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_beta_10.csv")
rownames(prom_cassettes) <- prom_cassettes$Cassette
prom_cassettes$Cassette <- NULL
prom_cassettes <- prom_cassettes[, x$PD_ID]

# Epitype annotations
my_annotations <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")

#
# ASSOCIATION OF CASSETTES WITH TILS
#

# Getting TILs
main_var <- as.numeric(x[colnames(prom_cassettes),"TILs"])

# PROMOTER

# Define the variables to test
variables <- rownames(prom_cassettes)

# Compute correlations and wilocoxon p values. Kendall + Distance
results <- data.frame(
  Cassette = variables,
  Wilcoxon_P_value = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    if (all(is.na(y))) return(NA)
    
    # Keep only non-NA pairs
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    
    if (length(unique(y)) < 2) return(NA)
    
    # K-means clustering into 2 groups
    cluster <- kmeans(y, centers = 2)$cluster
    
    # Wilcoxon test between groups for TIL levels
    wilcox.test(x ~ cluster)$p.value
  }),
  
  Kendall_Tau = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    if (all(is.na(y))) return(NA)
    
    # Keep only non-NA pairs
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    
    # Kendall's Tau correlation
    cor.test(x, y, method = "kendall")$estimate
  }),
  
  Distance_Correlation = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    if (all(is.na(y))) return(NA)
    
    # Keep only non-NA pairs
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    
    # Distance Correlation
    dcor(x, y)
  }),
  
  Median_TIL_Hypo = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    
    cluster <- kmeans(y, centers = 2)$cluster
    
    # Selecting hypermethylated cluster
    if (mean(y[cluster == 1]) > mean(y[cluster == 2])) {
      mean(x[cluster == 1])
    } else {
      median(x[cluster == 2])
    }
    
    
  }),
  
  Median_TIL_Hyper = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    if (length(unique(y)) < 2) return(NA)
    cluster <- kmeans(y, centers = 2)$cluster
    
    # Selecting hypermethylated cluster
    if (mean(y[cluster == 2]) > mean(y[cluster == 1])) {
      mean(x[cluster == 1])
    } else {
      median(x[cluster == 2])
    }
  }),
  
  gene = sapply(variables, function(v) {
    data <- promoter_10$colors
    cpgs <- names(data)[data == v]
    paste0(unique(unname(genes[cpgs])), collapse = ",")
  }) 
  
)

# Apply multiple testing correction
results$P_adj_Bonf <- p.adjust(results$Wilcoxon_P_value, method = "bonferroni")
results$P_adj_FDR <- p.adjust(results$Wilcoxon_P_value, method = "fdr")

# Convert p-values to -log10 scale
results$logP_FDR <- -log10(results$P_adj_FDR)
results$logP_bonferroni <- -log10(results$P_adj_Bonf)

#Change rownames
rownames(results) <- results$Cassette


# Set threshold and number of top hits to label
threshold <- 0.05
top_n <- 30

# Subset the top significant hits to label
top_labels <- results %>%
  filter(P_adj_FDR < threshold) %>%
  arrange(desc(logP_FDR)) %>%
  slice_head(n = top_n)

# Volcano plot
ggplot(results, aes(x = Kendall_Tau, y = logP_FDR)) +
  geom_point(aes(color = P_adj_FDR < threshold), size = 3, alpha = 0.8) + 
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(
    data = top_labels,
    aes(label = Cassette),
    size = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", color = "blue") +
  theme_bw() +
  labs(
    x = "Kendall’s Tau",
    y = "-log10(FDR. P-value)"
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )

# Sorting results
results_sorted <- results[order(-abs(results$Kendall_Tau)), ]


#
# DEFINING IMPORTANT CASSETTES USING BORUTA ALGORITHM
#

# Getting data to run boruta
X <- t(prom_cassettes)
tils <- x$TILs

# Remove NAs (from TILs)
to_remove <- is.na(tils)

X <- X[!to_remove,]
tils <- tils[!to_remove]

# Combine into data frame
df_boruta <- data.frame(X)
df_boruta$TILs <- tils

# Run Boruta for feture selection
boruta_result <- Boruta(TILs ~ ., data = df_boruta, doTrace = 2, maxRuns = 500)

# Check results
print(boruta_result)

# Plot boruta results
plot(boruta_result, las = 2, cex.axis = 0.7)

# Get confirmed attributes
confirmed_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
cat("Confirmed important features:\n")
print(confirmed_features)

# Get feature importance scores
feature_importance <- attStats(boruta_result)


# Plot importance of each selected variable

# 1. Extract variable importance and convert to a tidy data frame
importance_df <- data.frame("Median_Importance" = feature_importance[confirmed_features, "medianImp"],
                            "Cassette" = rownames(feature_importance[confirmed_features, ]))

importance_df$Cassette <- substring(importance_df$Cassette, 2)

importance_df <- tibble(
  variable = importance_df$Cassette,
  importance = as.numeric(importance_df$Median_Importance)
)

# 2. Extract Kendall Tau named vector and convert to data frame
kendall_tau_vec <- results_sorted[, "Kendall_Tau"] 
names(kendall_tau_vec) <- rownames(results_sorted)
kendall_df <- tibble(
  variable = names(kendall_tau_vec),
  kendall_tau = as.numeric(kendall_tau_vec)
)

# 3. Join importance and Kendall Tau
merged_df <- importance_df %>%
  left_join(kendall_df, by = "variable")

# 4. Sort op importance
top_vars <- merged_df %>%
  arrange(desc(importance)) %>%
  mutate(variable = factor(variable, levels = rev(variable)))

# 5. Dot plot (left)
p1 <- ggplot(top_vars, aes(x = kendall_tau, y = variable)) +
  geom_point(aes(
    size = abs(kendall_tau),
    color = "red"
  )) +
  scale_size(range = c(0.5, 3)) +
  scale_color_identity() +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(-0.4, 0.4) +
  xlab("Kendall Tau")

# 6. Bar plot (right)
p2 <- ggplot(top_vars, aes(x = importance, y = variable)) +
  geom_col(fill = "steelblue") +
  theme_classic() +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab("Boruta Median Importance")

# 7. Combine plots
p1 + p2 + plot_layout(ncol = 2, widths = c(1, 2))


#
# PLOTTING CPGS AFFECTING 
#

data <- promoter_10$colors
cpgs <- names(data)[data == 10]
genes[cpgs]

pam50_annotations <- my_annotations[colnames(betaAdj), "PAM50"]
tnbc_annotation <- my_annotations[colnames(betaAdj), "TNBC"]
HRD_annotation <- my_annotations[colnames(betaAdj), "HRD"]
epi_annotation <- my_annotations[colnames(betaAdj), "NMF_atacDistal"]
im_annotation <- my_annotations[colnames(betaAdj), "IM"]
tils_annotation <- as.numeric(x[colnames(betaAdj), "TILs"])

# Create top anotation
top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                    TNBC = tnbc_annotation,
                                    HRD = HRD_annotation,
                                    IM = im_annotation,
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
                                      "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                      "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                      "IM"=c("Negative"="grey", "Positive"="black"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                      "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                  "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")
                                    )
)


# Updated left_annotation with color scale
right_annotation <- rowAnnotation(
  "Normal beta" = rowMeans(betaNorm[cpgs,]),
  col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")),
             "ATAC" = c("0" = "white", "1"= "black"))
)

# CpG context annotation
left_annotation <- rowAnnotation("Context"= annoObj$featureClass[annoObj$illuminaID %in% cpgs]
)

# Cluster based on methylation
cluster_promoter <- kmeans(t(betaAdj[cpgs,]), centers = 2)

# Determine hypo and hypermethylated cluster
promoter_state <- if (mean(betaAdj[cpgs,cluster_promoter$cluster==1]) >
                      mean(betaAdj[cpgs,cluster_promoter$cluster==2])) {
  
  as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
  
} else {
  
  as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
  
}


# Heatmap of genes
Heatmap(
  betaAdj[cpgs,],
  cluster_rows = FALSE,
  row_order = order(annoObj$start[annoObj$illuminaID %in% cpgs]),
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  column_split = promoter_state,
  top_annotation = top_annotation,
  right_annotation = right_annotation,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  name = "Tumor beta"
)
