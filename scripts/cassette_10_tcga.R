#! usr/bin/Rscript

library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(dplyr)
library(grid)

#
# LOADING DATA
#

set.seed(07102000)

# TCGA-BRCA data
tcga <- new.env()
load("/Users/isasiain/PhD/Projects/project_3/data/tcga_brca_withAnnotations.RData", envir = tcga)
load("/Users/isasiain/PhD/Projects/project_3/data/jvc_PAM50_NCN_subtype.RData", envir = tcga)
tcga$annoObj[, "featureClass"] <- sapply(tcga$annoObj[, "featureClass"], function(val) {
  if (val == "distal" | val == "distal body") {"Distal"}
  else if (val == "proximal dn" | val == "proximal up") {"Proximal"}
  else {"Promoter"}
})


# SCANB TNBC discovery cohort
scanb <- new.env()
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData", envir=scanb)
scanb$annoObj[, "featureClass"] <- sapply(scanb$annoObj[, "featureClass"], function(val) {
  if (val == "distal" | val == "distal body") {"Distal"}
  else if (val == "proximal dn" | val == "proximal up") {"Proximal"}
  else {"Promoter"}
})

# Distal cassettes
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")

# Reading in silico TILs
tcga_tils <- read.csv2("/Users/isasiain/PhD/Projects/project_3/data/tcga_tils_from_mpath.csv", skip = 1)
tcga_tils$sample_id <- sub("^([A-Z0-9-]+?-[A-Z0-9]+-[A-Z0-9]+)-.*$", "\\1", tcga_tils$Biospecimen_barcode)

# LINKING TCGA SAMPLES TO IN SILICO TILS
rownames(tcga_tils) <- tcga_tils$sample_id

tcga$sampleAnno[, "TILs"] <- tcga_tils[tcga$sampleAnno[, "bcr_patient_barcode"],"TIL_score"]

# Cassette 10
cpgs_10 <- names(promoter_10$colors)[promoter_10$colors == 10]
cpgs_10 <- cpgs_10[cpgs_10 %in% rownames(tcga$betaAdj)]

# PAM50 annotations (my.pam50.subtype)
load("/Users/isasiain/PhD/Projects/project_3/data/jvc_PAM50_NCN_subtype_TCGA.RData")


#
# TILs VS TCGA SUBTYPE
#

# Prepare the data frame for ggplot
df <- data.frame(
  TILs = as.numeric(tcga$sampleAnno[, "TILs"]),
  Subtype = as.factor(my.pam50.subtype[rownames(tcga$sampleAnno)])
)

# Drop NAs if any
df <- na.omit(df)

# Drop unclassified samples
df <- df[!df$Subtype == "unclassified",]

# Custom colors
subtype_colors <- c(
  "Her2" = "purple",
  "LumA" = "darkblue",
  "Basal" = "indianred1",
  "LumB" = "lightblue",
  "Normal" = "green"
)

# Plot
ggplot(df, aes(x = Subtype, y = TILs, fill = Subtype)) +
  geom_boxplot(width = 0.7, outlier.shape = 21, alpha = 0.9, ) +
  scale_fill_manual(values = subtype_colors) +
  theme_bw(base_size = 11) +
  labs(
    x = "PAM50 Subtype",
    y = "TILs (%)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) + 
  stat_compare_means(
    method = "kruskal.test",
    label.x = 4.5,
    label.y = max(df$TILs, na.rm = TRUE) * 0.9,   
    label = "p.format"
  )


#
# ANALYSES IN TNBC
#


# Subset the beta matrix
beta_mat <- tcga$betaAdj[cpgs_10, ]
beta_mat <- beta_mat[,tcga$sampleAnno$TNBC]

# Cluster methylation state of cassette
clusters_cassette_10 <- as.factor(kmeans(t(beta_mat), 2)$cluster)
cluster_means <- tapply(
  colMeans(beta_mat), 
  clusters_cassette_10,   
  mean                   
)
cluster_label <- ifelse(cluster_means[clusters_cassette_10] == min(cluster_means), "Hypo.", "Hyper.")

# COMPARE TILS BETWEEN CLUSTERS

# Prepare data for ggplot
tils_df <- data.frame(
  Sample = colnames(beta_mat),
  TILs = tcga$sampleAnno[colnames(beta_mat), "TILs"],
  Cluster = cluster_label
)


ggplot(tils_df, aes(x = Cluster, y = TILs)) +
  geom_violin(, fill="black") +
  geom_boxplot(width=0.1) +
  theme_bw(base_size = 14)+
  stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Hypo.", "Hyper.")),
      label = "p.format",
      label.y = max(tils_df$TILs, na.rm = TRUE) * 1.05, 
      size = 5) +
  labs(,
    x = "Cluster",
    y = "TILs (%)"
  ) +
  ylim(0,max(tils_df$TILs, na.rm = TRUE) * 1.2) +
  theme(axis.title.y = ggtext::element_markdown())



# ANALYSE METHYLATION PATTERNS

# Generate annotations
subtypes <- factor(
  my.pam50.subtype[colnames(beta_mat)],
  levels = c("unclassified","Normal","LumA","LumB","Her2","Basal")
)

top_annotation <- HeatmapAnnotation(
  PAM50 = subtypes,
  TILs = anno_points(tcga$sampleAnno[colnames(beta_mat), "TILs"]), # dot size
  #ER = tcga$sampleAnno[colnames(beta_mat),"ER"],
  #PR = tcga$sampleAnno[colnames(beta_mat),"PR"],
  #HER2 = tcga$sampleAnno[colnames(beta_mat),"HER2"],
  col = list(
    ER = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey"),
    PR = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey", "Indeterminate"="grey"),
    HER2 = c("Negative"="white", "Positive"="black", "NA"="grey", "Indeterminate"="grey", "Equivocal"="grey"),
    PAM50 = c("Her2"="purple","LumA"="darkblue","Basal"="indianred1","LumB"="lightblue","Normal"="green","unclassified"="grey")
  )
)


# Plotting heatmap
Heatmap(
  beta_mat,
  column_split = cluster_label,
  show_column_names = FALSE,
  show_row_names = FALSE,
  name = "Tumor\nBeta",
  top_annotation = top_annotation
)


# CHECK PER GENE
genes <- c("GBP4", "OAS2", "ZBP1", "CARD16")

list_of_tils_plots <- list()
list_of_gex_plots <- list()
list_of_heatmaps <- list()

for (gene in genes) {
  
  my_cpgs <- tcga$betaAdj[tcga$annoObj[,"nameUCSCknownGeneOverlap"] == gene,]
  
  # Subset the beta matrix
  beta_mat <- my_cpgs[,tcga$sampleAnno$TNBC,drop=FALSE]
  
  # Cluster methylation state of cassette
  clusters_cassette_10 <- as.factor(kmeans(t(beta_mat[rownames(beta_mat),, drop=FALSE]), 2)$cluster)
  cluster_means <- tapply(
    colMeans(beta_mat), 
    clusters_cassette_10,   
    mean                   
  )
  cluster_label <- ifelse(cluster_means[clusters_cassette_10] == min(cluster_means), "Hypo.", "Hyper.")
  
  # Getting gene expression (fpkm)
  gene_ensembl <- bitr(
    gene,
    fromType = "SYMBOL",
    toType = "ENSEMBL",
    OrgDb = org.Hs.eg.db
  )[1,2]
  
  gex_data <- tcga$gexFpkm[grep(gene_ensembl,rownames(tcga$gexFpkm)), colnames(beta_mat)]
  
  # Prepare data for ggplot
  tils_df <- data.frame(
    Sample = colnames(beta_mat),
    TILs = tcga$sampleAnno[colnames(beta_mat), "TILs"],
    Cluster = cluster_label,
    GEX = gex_data
  )
  
  # Plotting TILS
  list_of_tils_plots[[gene]] <- ggplot(tils_df, aes(x = Cluster, y = TILs)) +
    geom_violin(, fill="black") +
    geom_boxplot(width=0.1) +
    theme_bw(base_size = 14)+
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Hypo.", "Hyper.")),
      label = "p.format",
      label.y = max(tils_df$TILs, na.rm = TRUE) * 1.05, # position above violins
      size = 5) +
    labs(,
         x = "Cluster",
         y = "TILs (%)"
    ) +
    ylim(0, max(tils_df$TILs, na.rm = TRUE) * 1.25) +
    theme(axis.title.y = ggtext::element_markdown())
  
  # Plotting Expression
  list_of_gex_plots[[gene]] <- ggplot(tils_df, aes(x = Cluster, y = GEX)) +
    geom_violin(, fill="black") +
    geom_boxplot(width=0.1) +
    theme_bw(base_size = 14)+
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Hypo.", "Hyper.")),
      label = "p.format",
      label.y = max(tils_df$GEX, na.rm = TRUE) * 1.2, # position above violins
      size = 5) +
    labs(,
         x = "Cluster",
         y = paste0(gene, " FPKM")
    ) +
    ylim(0, max(tils_df$GEX, na.rm = TRUE) * 1.35) +
    theme(axis.title.y = ggtext::element_markdown())
  
  
  # Generate annotations
  subtypes <- factor(
    my.pam50.subtype[colnames(beta_mat)],
    levels = c("unclassified","Normal","LumA","LumB","Her2","Basal")
  )
  
  #Getting betas to plot
  cpgs_overlapping <- tcga$annoObj[tcga$annoObj$nameUCSCknownGeneOverlap == gene, "illuminaID"]
  betas_to_plot <- tcga$betaAdj[cpgs_overlapping, colnames(beta_mat)]
  
  top_annotation <- HeatmapAnnotation(
    PAM50 = subtypes,
    TILs = anno_points(tcga$sampleAnno[colnames(beta_mat), "TILs"]), # dot size
    #ER = tcga$sampleAnno[colnames(beta_mat),"ER"],
    #PR = tcga$sampleAnno[colnames(beta_mat),"PR"],
    #HER2 = tcga$sampleAnno[colnames(beta_mat),"HER2"],
    col = list(
      ER = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey"),
      PR = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey", "Indeterminate"="grey"),
      HER2 = c("Negative"="white", "Positive"="black", "NA"="grey", "Indeterminate"="grey", "Equivocal"="grey"),
      PAM50 = c("Her2"="purple","LumA"="darkblue","Basal"="indianred1","LumB"="lightblue","Normal"="green","unclassified"="grey")
    )
  )
  
  bottom_annotation <- HeatmapAnnotation(TILs = anno_barplot(tils_df[colnames(beta_mat), "GEX"]))
  

  
  right_annotation = rowAnnotation("In Cassette" = cpgs_overlapping %in% cpgs_10,
                                        col = list(
                                          "In Cassette" = c("TRUE" = "black", 
                                                            "FALSE" = "white")
                                          )
                                        )
  
  # Plotting heatmap
  list_of_heatmaps[[gene]] <- Heatmap(
    betas_to_plot,
    column_split = cluster_label,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    name = "Tumor\nBeta",
    top_annotation = top_annotation,
    bottom_annotation = bottom_annotation,
    left_annotation = right_annotation
  )
  
}

#Plotting tils vs methylation state
list_of_tils_plots[[1]] | list_of_tils_plots[[2]] | list_of_tils_plots[[3]] | list_of_tils_plots[[4]]

# Plotting gex vs methylation state
list_of_gex_plots[[1]] | list_of_gex_plots[[2]] | list_of_gex_plots[[3]] | list_of_gex_plots[[4]]

#Heatmaps of genes
heatmap_grobs <- lapply(list_of_heatmaps, function(ht) {
  grid.grabExpr(draw(ht,
                     show_heatmap_legend = FALSE,
                     show_annotation_legend = FALSE))
})

wrap_plots(heatmap_grobs, ncol = 4)



#
# ANALYSES IN ALL BRCA
#


# Subset the beta matrix
beta_mat <- tcga$betaAdj[cpgs_10, ,drop=FALSE]
#beta_mat <- beta_mat[, my.pam50.subtype[colnames(beta_mat)] == "LumA",drop=FALSE]

# Cluster methylation state of cassette
clusters_cassette_10 <- as.factor(kmeans(t(beta_mat), 2)$cluster)
cluster_means <- tapply(
  colMeans(beta_mat), 
  clusters_cassette_10,   
  mean                   
)
cluster_label <- ifelse(cluster_means[clusters_cassette_10] == min(cluster_means), "Hypo.", "Hyper.")

# COMPARE TILS BETWEEN CLUSTERS

# Prepare data for ggplot
tils_df <- data.frame(
  Sample = colnames(beta_mat),
  TILs = tcga$sampleAnno[colnames(beta_mat), "TILs"],
  Cluster = cluster_label
)


ggplot(tils_df, aes(x = Cluster, y = TILs)) +
  geom_boxplot(fill="lightgrey") +
  theme_bw(base_size = 14)+
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Hypo.", "Hyper.")),
    label = "p.format",
    label.y = max(tils_df$TILs, na.rm = TRUE) * 1.05, # position above violins
    size = 5) +
  labs(,
       x = "Cluster",
       y = "TILs (%)"
  ) +
  ylim(0, max(tils_df$TILs, na.rm = TRUE) * 1.2) +
  theme(axis.title.y = ggtext::element_markdown())


# ANALYSE METHYLATION PATTERNS

# Generate annotations
subtypes <- factor(
  my.pam50.subtype[colnames(beta_mat)],
  levels = c("unclassified","Normal","LumA","LumB","Her2","Basal")
)

top_annotation <- HeatmapAnnotation(
  PAM50 = subtypes,
  TILs = anno_points(tcga$sampleAnno[colnames(beta_mat), "TILs"]), # dot size
  ER = tcga$sampleAnno[colnames(beta_mat),"ER"],
  PR = tcga$sampleAnno[colnames(beta_mat),"PR"],
  HER2 = tcga$sampleAnno[colnames(beta_mat),"HER2"],
  col = list(
    ER = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey"),
    PR = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey", "Indeterminate"="grey"),
    HER2 = c("Negative"="white", "Positive"="black", "NA"="grey", "Indeterminate"="grey", "Equivocal"="grey"),
    PAM50 = c("Her2"="purple","LumA"="darkblue","Basal"="indianred1","LumB"="lightblue","Normal"="green","unclassified"="grey")
  )
)


# Plotting heatmap
Heatmap(
  beta_mat,
  column_split = cluster_label,
  show_column_names = FALSE,
  show_row_names = FALSE,
  name = "Tumor\nBeta",
  top_annotation = top_annotation
)


# CLUSTERING BASED ON ALL CPGS

genes <- c("GBP4", "OAS2", "ZBP1", "CARD16")
subtypes <- c("LumA", "LumB", "Her2", "Basal", "Normal")

# To store plots
list_of_tils_plots <- list()
list_of_gex_plots <- list()
list_of_heatmaps <- list()
list_of_barplots <- list()

# To store proportion of hypomethylation
hypo_prop_summary <- data.frame()
all_data_long_df <- data.frame()

for (gene in genes) {
  
  # --- select all promoter and proximal CpGs for this gene ---
  my_cpgs <- rownames(tcga$annoObj)[
    tcga$annoObj$nameUCSCknownGeneOverlap == gene &
      tcga$annoObj$featureClass %in% c("Promoter", "Proximal", "Distal")
  ]
  
  # Subset beta matrix
  beta_mat <- tcga$betaAdj[my_cpgs, , drop = FALSE]
  
  # Cluster methylation pattern
  clusters_gene <- as.factor(kmeans(t(beta_mat), 2)$cluster)
  cluster_means <- tapply(colMeans(beta_mat), clusters_gene, mean)
  cluster_label <- ifelse(
    cluster_means[clusters_gene] == min(cluster_means),
    "Hypo.",
    "Hyper."
  )
  names(cluster_label) <- names(clusters_gene)
  
  # Get expression data
  gene_ensembl <- bitr(
    gene,
    fromType = "SYMBOL",
    toType = "ENSEMBL",
    OrgDb = org.Hs.eg.db
  )[1, 2]
  
  gex_data <- tcga$gexFpkm[grep(gene_ensembl, rownames(tcga$gexFpkm)), colnames(beta_mat)]
  
  # Prepare heatmap
  cpgs_overlapping <- my_cpgs
  betas_to_plot <- tcga$betaAdj[cpgs_overlapping, colnames(beta_mat), drop = FALSE]
  
  top_annotation <- HeatmapAnnotation(
    TILs = anno_points(tcga$sampleAnno[colnames(beta_mat), "TILs"], 
                       size = unit(1, "mm")),
    PAM50 = my.pam50.subtype[colnames(beta_mat)],
    ER = tcga$sampleAnno[colnames(beta_mat), "ER"],
    PR = tcga$sampleAnno[colnames(beta_mat), "PR"],
    HER2 = tcga$sampleAnno[colnames(beta_mat), "HER2"],
    col = list(
      ER = c("Negative" = "white", "Positive" = "black", "[Not Evaluated]" = "grey"),
      PR = c("Negative" = "white", "Positive" = "black", "[Not Evaluated]" = "grey", "Indeterminate" = "grey"),
      HER2 = c("Negative" = "white", "Positive" = "black", "NA" = "grey", "Indeterminate" = "grey", "Equivocal" = "grey"),
      PAM50 = c("Her2"="purple","LumA"="darkblue","Basal"="indianred1","LumB"="lightblue","Normal"="green","unclassified"="grey")
    )
  )
  
  bottom_annotation <- HeatmapAnnotation(GEX = anno_barplot(gex_data[colnames(tcga$betaAdj)]))
  
  row_annotation <- rowAnnotation(
    "CpG in Cassette" = rownames(beta_mat) %in% cpgs_10,
    col = list(
      "CpG in Cassette" = c("TRUE" = "black", "FALSE" = "white")
    )
  )
  
  list_of_heatmaps[[gene]]<- Heatmap(
    betas_to_plot,
    column_split = cluster_label[colnames(beta_mat)],
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    name = "Tumor\nBeta",
    top_annotation = top_annotation,
    left_annotation = row_annotation,
    bottom_annotation = bottom_annotation
  )
  
  # Iterate over subtypes
  for (subtype in subtypes) {
    samples_in_subtype <- names(my.pam50.subtype[colnames(beta_mat)][my.pam50.subtype[colnames(beta_mat)] == subtype])
    
    # Prepare dataframe
    tils_df <- data.frame(
      Gene = gene,
      Subtype = subtype,
      Sample = samples_in_subtype,
      TILs = tcga$sampleAnno[samples_in_subtype, "TILs"],
      Cluster = cluster_label[samples_in_subtype],
      GEX = gex_data[samples_in_subtype]
    )
    
    all_data_long_df <- rbind(all_data_long_df, tils_df)
    
    # TILs plots
    list_of_tils_plots[[gene]][[subtype]] <- ggplot(tils_df, aes(x = Cluster, y = TILs)) +
      geom_violin(fill = "black", size = 0.1) +
      geom_boxplot(width = 0.1, size = 0.1) +
      theme_bw(base_size = 14) +
      stat_compare_means(
        method = "wilcox.test",
        comparisons = list(c("Hypo.", "Hyper.")),
        label = "p.format",
        label.y = max(tils_df$TILs, na.rm = TRUE) * 1.05,
        size = 5
      ) +
      labs(x = "Cluster", y = "TILs (%)") +
      ylim(0, max(tils_df$TILs, na.rm = TRUE) * 1.2) +
      theme(axis.title.y = ggtext::element_markdown())
    
    # GEX plots
    list_of_gex_plots[[gene]][[subtype]] <- ggplot(tils_df, aes(x = Cluster, y = GEX)) +
      geom_violin(fill = "black", size = 0.1) +
      geom_boxplot(width = 0.1, size = 0.1) +
      theme_bw(base_size = 14) +
      stat_compare_means(
        method = "wilcox.test",
        comparisons = list(c("Hypo.", "Hyper.")),
        label = "p.format",
        label.y = max(tils_df$GEX, na.rm = TRUE) * 1.05,
        size = 5
      ) +
      labs(x = "Cluster", y = paste0(gene, " FPKM")) +
      ylim(0, max(tils_df$GEX, na.rm = TRUE) * 1.2) +
      theme(axis.title.y = ggtext::element_markdown())
    
    
    # Proportion of hypomethylation per subtype
    n_sub <- length(samples_in_subtype)
    n_hypo <- sum(cluster_label[samples_in_subtype] == "Hypo.", na.rm = TRUE)
    prop_df <- data.frame(
      Gene = gene,
      PAM50 = subtype,
      n = n_sub,
      n_Hypo. = n_hypo,
      Hypo. = if (n_sub > 0) n_hypo / n_sub else NA,
      stringsAsFactors = FALSE
    )
    hypo_prop_summary <- rbind(hypo_prop_summary, prop_df)
  }
}

# Plotting TILs vs subtype and gene
(list_of_tils_plots$CARD16$Basal | list_of_tils_plots$CARD16$Her2 | list_of_tils_plots$CARD16$LumA  | list_of_tils_plots$CARD16$LumB | list_of_tils_plots$CARD16$Normal)/
  (list_of_tils_plots$GBP4$Basal | list_of_tils_plots$GBP4$Her2 | list_of_tils_plots$GBP4$LumA  | list_of_tils_plots$GBP4$LumB | list_of_tils_plots$GBP4$Normal)/
  (list_of_tils_plots$OAS2$Basal | list_of_tils_plots$OAS2$Her2 | list_of_tils_plots$OAS2$LumA  | list_of_tils_plots$OAS2$LumB | list_of_tils_plots$OAS2$Normal)/
  (list_of_tils_plots$ZBP1$Basal | list_of_tils_plots$ZBP1$Her2 | list_of_tils_plots$ZBP1$LumA  | list_of_tils_plots$ZBP1$LumB | list_of_tils_plots$ZBP1$Normal)

# Plotting GEX vs subtype and gene methylation
(list_of_gex_plots$CARD16$Basal | list_of_gex_plots$CARD16$Her2 | list_of_gex_plots$CARD16$LumA  | list_of_gex_plots$CARD16$LumB | list_of_gex_plots$CARD16$Normal)/
  (list_of_gex_plots$GBP4$Basal | list_of_gex_plots$GBP4$Her2 | list_of_gex_plots$GBP4$LumA  | list_of_gex_plots$GBP4$LumB | list_of_gex_plots$GBP4$Normal)/
  (list_of_gex_plots$OAS2$Basal | list_of_gex_plots$OAS2$Her2 | list_of_gex_plots$OAS2$LumA  | list_of_gex_plots$OAS2$LumB | list_of_gex_plots$OAS2$Normal)/
  (list_of_gex_plots$ZBP1$Basal | list_of_gex_plots$ZBP1$Her2 | list_of_gex_plots$ZBP1$LumA  | list_of_gex_plots$ZBP1$LumB | list_of_gex_plots$ZBP1$Normal)


# Plotting piecharts
pie_data <- hypo_prop_summary %>%
  mutate(
    Hyper. = 1 - Hypo.
  ) %>%
  tidyr::pivot_longer(cols = c("Hypo.", "Hyper."), 
                      names_to = "State", values_to = "Proportion")


ggplot(pie_data, aes(x = "", y = Proportion, fill = State)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_grid(Gene ~ PAM50) +
  scale_fill_manual(
    values = c("Hypo." = "steelblue", "Hyper." = "indianred1"),
    name = NULL
  ) +
  theme_void(base_size = 14) +
  labs(
  ) +
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) 


# Analyse statistical significance
counts_df <- all_data_long_df %>%
  count(Gene, PAM50 = Subtype, Cluster) %>%
  tidyr::pivot_wider(names_from = PAM50, values_from = n, values_fill = 0)

chi_results <- all_data_long_df %>%
  group_by(Gene) %>%
  summarise(
    p_value = tryCatch({
      tbl <- table(Cluster, Subtype)
      chisq.test(tbl)$p.value
    }, error = function(e) NA)
  )

print(chi_results)

#Heatmaps
heatmap_grobs <- lapply(list_of_heatmaps, function(ht) {
  grid.grabExpr(draw(ht,
                show_heatmap_legend = FALSE,
                show_annotation_legend = FALSE))
})

wrap_plots(heatmap_grobs, ncol = 4)


# Plotting co-ocurrence of hypomethylation
gene_sets <- all_data_long_df %>%
  filter(Cluster == "Hyper.") %>%
  group_by(Gene) %>%
  summarize(Samples = list(Sample)) %>%
  deframe()  # produces a named list: Gene -> vector of Samples

venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,          # NULL to plot directly in R
  fill = c("steelblue", "indianred1", "gold", "forestgreen"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  margin = 0.1
)

grid.draw(venn.plot)


gene_sets <- all_data_long_df %>%
  filter(Cluster == "Hypo.") %>%
  group_by(Gene) %>%
  summarize(Samples = list(Sample)) %>%
  deframe()  # produces a named list: Gene -> vector of Samples

venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,          # NULL to plot directly in R
  fill = c("steelblue", "indianred1", "gold", "forestgreen"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  margin = 0.1
)

grid.draw(venn.plot)
