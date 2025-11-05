#! usr/bin/Rscript

library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

#
# LOADING DATA
#

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

# Distal cassettes
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")

# Reading in silico TILs
tcga_tils <- read.csv2("/Users/isasiain/PhD/Projects/project_3/data/tcga_tils_from_mpath.csv", skip = 1)
tcga_tils$sample_id <- sub("^([A-Z0-9-]+?-[A-Z0-9]+-[A-Z0-9]+)-.*$", "\\1", tcga_tils$Biospecimen_barcode)

# LINKING TCGA SAMPLES TO IN SILICO TILS
rownames(tcga_tils) <- tcga_tils$sample_id

tcga$sampleAnno[, "TILs"] <- tcga_tils[tcga$sampleAnno[, "bcr_patient_barcode"],"TIL_score"]

# Cassette 10
cpgs_10 <- names(distal_10$colors)[distal_10$colors == 10]
cpgs_10 <- cpgs_10[cpgs_10 %in% rownames(tcga$betaAdj)]

# PAM50 annotations
load("/Users/isasiain/PhD/Projects/project_3/data/jvc_PAM50_NCN_subtype_TCGA.RData")


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
cluster_label <- ifelse(cluster_means[clusters_cassette_10] == min(cluster_means), "Hypomethylated", "Hypermethylated")

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
      comparisons = list(c("Hypomethylated", "Hypermethylated")),
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
  tcga$my.pam50.subtype[colnames(beta_mat)],
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
  
  my_cpgs <- cpgs_10[tcga$annoObj[cpgs_10,"nameUCSCknownGeneOverlap"] == gene]
  
  # Subset the beta matrix
  beta_mat <- tcga$betaAdj[my_cpgs, ,drop=FALSE]
  beta_mat <- beta_mat[,tcga$sampleAnno$TNBC,drop=FALSE]
  
  # Cluster methylation state of cassette
  clusters_cassette_10 <- as.factor(kmeans(t(beta_mat), 2)$cluster)
  cluster_means <- tapply(
    colMeans(beta_mat), 
    clusters_cassette_10,   
    mean                   
  )
  cluster_label <- ifelse(cluster_means[clusters_cassette_10] == min(cluster_means), "Hypomethylated", "Hypermethylated")
  
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
      comparisons = list(c("Hypomethylated", "Hypermethylated")),
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
      comparisons = list(c("Hypomethylated", "Hypermethylated")),
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
    tcga$my.pam50.subtype[colnames(beta_mat)],
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
  
  row_annotation <- rowAnnotation(
    "CpG context" = tcga$annoObj[cpgs_overlapping, "featureClass"],
    col = list(
      "CpG context" = c("Distal" = "#f8766d", 
        "Promoter" = "#00ba38", 
        "Proximal" = "#619cff")
    )
  )
  
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
    name = "Tumor\nBeta",
    top_annotation = top_annotation,
    bottom_annotation = bottom_annotation,
    left_annotation = row_annotation,
    right_annotation = right_annotation
  )
  
}

list_of_tils_plots[[1]] | list_of_tils_plots[[2]] | list_of_tils_plots[[3]] | list_of_tils_plots[[4]]
list_of_gex_plots[[1]] | list_of_gex_plots[[2]] | list_of_gex_plots[[3]] | list_of_gex_plots[[4]]

list_of_heatmaps[[1]]
list_of_heatmaps[[2]]
list_of_heatmaps[[3]]
list_of_heatmaps[[4]]


#
# ANALYSES IN ALL BRCA
#


# Subset the beta matrix
beta_mat <- tcga$betaAdj[cpgs_10, ,drop=FALSE]
#beta_mat <- beta_mat[, tcga$my.pam50.subtype[colnames(beta_mat)] == "LumA",drop=FALSE]

# Cluster methylation state of cassette
clusters_cassette_10 <- as.factor(kmeans(t(beta_mat), 2)$cluster)
cluster_means <- tapply(
  colMeans(beta_mat), 
  clusters_cassette_10,   
  mean                   
)
cluster_label <- ifelse(cluster_means[clusters_cassette_10] == min(cluster_means), "Hypomethylated", "Hypermethylated")

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
    comparisons = list(c("Hypomethylated", "Hypermethylated")),
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
  tcga$my.pam50.subtype[colnames(beta_mat)],
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


# CHECK PER GENE
genes <- c("GBP4", "OAS2", "ZBP1", "CARD16")

list_of_tils_plots <- list()
list_of_gex_plots <- list()
list_of_heatmaps <- list()
list_of_barplots <- list()

for (gene in genes) {
  
  my_cpgs <- cpgs_10[tcga$annoObj[cpgs_10,"nameUCSCknownGeneOverlap"] == gene]
  
  # Subset the beta matrix
  beta_mat <- tcga$betaAdj[my_cpgs, ,drop=FALSE]
  #beta_mat <- beta_mat[, tcga$my.pam50.subtype[colnames(beta_mat)] == "LumB",drop=FALSE]
  
  # Cluster methylation state of cassette
  clusters_cassette_10 <- as.factor(kmeans(t(beta_mat), 2)$cluster)
  cluster_means <- tapply(
    colMeans(beta_mat), 
    clusters_cassette_10,   
    mean                   
  )
  cluster_label <- ifelse(cluster_means[clusters_cassette_10] == min(cluster_means), "Hypomethylated", "Hypermethylated")
  
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
    GEX = gex_data,
    PAM50 = as.factor(tcga$my.pam50.subtype[colnames(beta_mat)])
  )
  
  # Plotting TILS
  list_of_tils_plots[[gene]] <- ggplot(tils_df, aes(x = Cluster, y = TILs)) +
    geom_violin(, fill="black") +
    geom_boxplot(width=0.1) +
    theme_bw(base_size = 14)+
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Hypomethylated", "Hypermethylated")),
      label = "p.format",
      label.y = max(tils_df$TILs, na.rm = TRUE) * 1.05, # position above violins
      size = 5) +
    labs(,
         x = "Cluster",
         y = "TILs (%)"
    ) +
    ylim(0, max(tils_df$TILs, na.rm = TRUE) * 1.2) +
    theme(axis.title.y = ggtext::element_markdown())
  
  # Plotting Expression
  list_of_gex_plots[[gene]] <- ggplot(tils_df, aes(x = Cluster, y = GEX)) +
    geom_violin(, fill="black") +
    geom_boxplot(width=0.1) +
    theme_bw(base_size = 14)+
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Hypomethylated", "Hypermethylated")),
      label = "p.format",
      label.y = max(tils_df$GEX, na.rm = TRUE) * 1.05, # position above violins
      size = 5) +
    labs(,
         x = "Cluster",
         y = paste0(gene, " FPKM")
    ) +
    ylim(0, max(tils_df$GEX, na.rm = TRUE) * 1.2) +
    theme(axis.title.y = ggtext::element_markdown())
  
  
  # Generate annotations
  subtypes <- factor(
    tcga$my.pam50.subtype[colnames(beta_mat)],
    levels = c("unclassified","Normal","LumA","LumB","Her2","Basal")
  )
  
  #Getting betas to plot
  cpgs_overlapping <- tcga$annoObj[tcga$annoObj$nameUCSCknownGeneOverlap == gene, "illuminaID"]
  betas_to_plot <- tcga$betaAdj[cpgs_overlapping, colnames(beta_mat)]
  
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
  
  bottom_annotation <- HeatmapAnnotation(TILs = anno_barplot(tils_df[colnames(beta_mat), "GEX"]))
  
  row_annotation <- rowAnnotation(
    "CpG context" = tcga$annoObj[cpgs_overlapping, "featureClass"],
    col = list(
      "CpG context" = c("Distal" = "#f8766d", 
                        "Promoter" = "#00ba38", 
                        "Proximal" = "#619cff"),
      "ATAC" = c("0"="white", "1"="black")
    )
  )
  
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
    name = "Tumor\nBeta",
    top_annotation = top_annotation,
    bottom_annotation = bottom_annotation,
    left_annotation = row_annotation,
    right_annotation = right_annotation
  )
  
  # BARPLOT OF METHYLATION PER SUBTYPE
  list_of_barplots[[gene]] <- ggplot(data=tils_df, aes(x=PAM50, fill=Cluster)) +
    geom_bar(position = "fill") +
    theme_bw(base_size = 14)
  
}

list_of_tils_plots[[1]] | list_of_tils_plots[[2]] | list_of_tils_plots[[3]] | list_of_tils_plots[[4]]
list_of_gex_plots[[1]] | list_of_gex_plots[[2]] | list_of_gex_plots[[3]] | list_of_gex_plots[[4]]

list_of_heatmaps[[1]]
list_of_heatmaps[[2]]
list_of_heatmaps[[3]]
list_of_heatmaps[[4]]

list_of_barplots[[1]] / list_of_barplots[[2]] / list_of_barplots[[3]] / list_of_barplots[[4]]



#
# ANALYSES PER PAM50 SUBTYPE
#


# List of PAM50 subtypes to analyse
pam50_subtypes <- c("LumA", "LumB", "Her2", "Basal", "Normal")

# Genes of interest
genes <- c("GBP4", "OAS2", "ZBP1", "CARD16")

# Storage for results
results <- list()

for (subtype in pam50_subtypes) {
  message("Processing subtype: ", subtype)
  
  # Subset the beta matrix
  beta_mat <- tcga$betaAdj[cpgs_10, names(tcga$my.pam50.subtype)[tcga$my.pam50.subtype == subtype], drop = FALSE]
  
  # Cluster methylation state of cassette
  clusters_cassette_10 <- as.factor(kmeans(t(beta_mat), 2)$cluster)
  cluster_means <- tapply(colMeans(beta_mat), clusters_cassette_10, mean)
  cluster_label <- ifelse(cluster_means[clusters_cassette_10] == min(cluster_means),
                          "Hypomethylated", "Hypermethylated")
  
  # Compare TILs between clusters
  tils_df <- data.frame(
    Sample = colnames(beta_mat),
    TILs = tcga$sampleAnno[colnames(beta_mat), "TILs"],
    Cluster = cluster_label
  )
  
  p_tils <- ggplot(tils_df, aes(x = Cluster, y = TILs)) +
    geom_violin(fill = "black") +
    geom_boxplot(width = 0.1) +
    theme_bw(base_size = 14) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Hypomethylated", "Hypermethylated")),
      label = "p.format",
      label.y = max(tils_df$TILs, na.rm = TRUE) * 1.05,
      size = 5
    ) +
    labs(title = paste0(subtype, " subtype"),
         x = "Cluster",
         y = "TILs (%)") +
    ylim(0, max(tils_df$TILs, na.rm = TRUE) * 1.25)
  
  # Save results for each gene in this subtype
  list_of_tils_plots <- list()
  list_of_gex_plots <- list()
  list_of_heatmaps <- list()
  
  for (gene in genes) {
    my_cpgs <- cpgs_10[tcga$annoObj[cpgs_10, "nameUCSCknownGeneOverlap"] == gene]
    beta_gene <- tcga$betaAdj[my_cpgs, names(tcga$my.pam50.subtype)[tcga$my.pam50.subtype == subtype], drop = FALSE]
    
    # Skip genes with too few CpGs
    if (nrow(beta_gene) == 0 || ncol(beta_gene) < 5) next
    
    # Cluster methylation
    clusters_gene <- as.factor(kmeans(t(beta_gene), 2)$cluster)
    cluster_means <- tapply(colMeans(beta_gene), clusters_gene, mean)
    cluster_label <- ifelse(cluster_means[clusters_gene] == min(cluster_means),
                            "Hypomethylated", "Hypermethylated")
    
    # Get gene expression
    gene_ensembl <- bitr(gene, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)[1, 2]
    gex_data <- tcga$gexFpkm[grep(gene_ensembl, rownames(tcga$gexFpkm)), colnames(beta_gene)]
    
    tils_gene <- data.frame(
      Sample = colnames(beta_gene),
      TILs = tcga$sampleAnno[colnames(beta_gene), "TILs"],
      Cluster = cluster_label,
      GEX = gex_data
    )
    
    # TIL violin
    list_of_tils_plots[[gene]] <- ggplot(tils_gene, aes(x = Cluster, y = TILs)) +
      geom_violin(fill = "black") +
      geom_boxplot(width = 0.1) +
      theme_bw(base_size = 14) +
      stat_compare_means(
        method = "wilcox.test",
        comparisons = list(c("Hypomethylated", "Hypermethylated")),
        label = "p.format",
        label.y = max(tils_gene$TILs, na.rm = TRUE) * 1.05,
        size = 5
      ) +
      labs(title = paste0(gene, " (", subtype, ")"),
           x = "Cluster", y = "TILs (%)") +
      ylim(0, max(tils_gene$TILs, na.rm = TRUE) * 1.25)
    
    # GEX violin
    list_of_gex_plots[[gene]] <- ggplot(tils_gene, aes(x = Cluster, y = GEX)) +
      geom_violin(fill = "black") +
      geom_boxplot(width = 0.1) +
      theme_bw(base_size = 14) +
      stat_compare_means(
        method = "wilcox.test",
        comparisons = list(c("Hypomethylated", "Hypermethylated")),
        label = "p.format",
        label.y = max(tils_gene$GEX, na.rm = TRUE) * 1.2,
        size = 5
      ) +
      labs(title = paste0(gene, " (", subtype, ")"),
           x = "Cluster", y = paste0(gene, " FPKM")) +
      ylim(0, max(tils_gene$GEX, na.rm = TRUE) * 1.35)
    
    # Heatmap
    cpgs_overlapping <- tcga$annoObj[tcga$annoObj$nameUCSCknownGeneOverlap == gene, "illuminaID"]
    betas_to_plot <- tcga$betaAdj[cpgs_overlapping, colnames(beta_gene)]
    
    row_annotation <- rowAnnotation(
      "CpG context" = as.factor(as.vector(tcga$annoObj[cpgs_overlapping, "featureClass"])),
      col = list("CpG context" = c("Distal" = "#f8766d",
                                   "Promoter" = "#00ba38",
                                   "Proximal" = "#619cff"))
    )
    
    list_of_heatmaps[[gene]] <- Heatmap(
      betas_to_plot,
      column_split = cluster_label,
      show_column_names = FALSE,
      show_row_names = FALSE,
      name = paste0("Beta (", subtype, ")"),
      left_annotation = row_annotation
    )
  }
  
  # Store results
  results[[subtype]] <- list(
    TILs = p_tils,
    TILs_genes = list_of_tils_plots,
    GEX_genes = list_of_gex_plots,
    Heatmaps = list_of_heatmaps
  )
}

results$LumA$TILs_genes$GBP4 | results$LumA$TILs_genes$OAS2 | results$LumA$TILs_genes$ZBP1 | results$LumA$TILs_genes$CARD16
results$LumB$TILs_genes$GBP4 | results$LumB$TILs_genes$OAS2 | results$LumB$TILs_genes$ZBP1 | results$LumB$TILs_genes$CARD16
results$Her2$TILs_genes$GBP4 | results$Her2$TILs_genes$OAS2 | results$Her2$TILs_genes$ZBP1 | results$Her2$TILs_genes$CARD16
results$Normal$TILs_genes$GBP4 | results$Normal$TILs_genes$OAS2 | results$Normal$TILs_genes$ZBP1 | results$Normal$TILs_genes$CARD16
results$Basal$TILs_genes$GBP4 | results$Basal$TILs_genes$OAS2 | results$Basal$TILs_genes$ZBP1 | results$Basal$TILs_genes$CARD16
