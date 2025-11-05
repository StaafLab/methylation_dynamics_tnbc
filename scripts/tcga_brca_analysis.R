#! usr/bin/Rscript

library(ComplexHeatmap)
library(ggplot2)

#
# LOADING DATA
#

# TCGA-BRCA data
tcga <- new.env()
load("/Users/isasiain/PhD/Projects/project_3/data/tcga_brca_withAnnotations.RData", envir = tcga)
load("/Users/isasiain/PhD/Projects/project_3/data/jvc_PAM50_NCN_subtype.RData", envir = tcga)

# SCANB TNBC discovery cohort
scanb <- new.env()
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData", envir=scanb)

# Distal cassettes
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")


# Plotting heatmap of distal cassettes 1 and 2 in TCGA-BRCA
# Using CpG annotations form EPIC

# Cassette 1
cpgs_1 <- names(distal_10$colors)[distal_10$colors == 1]
cpgs_1 <- cpgs_1[cpgs_1 %in% rownames(tcga$betaAdj)]

cpgs_2 <- names(distal_10$colors)[distal_10$colors == 2]
cpgs_2 <- cpgs_2[cpgs_2 %in% rownames(tcga$betaAdj)]

cpgs_3 <- names(distal_10$colors)[distal_10$colors == 3]
cpgs_3 <- cpgs_3[cpgs_3 %in% rownames(tcga$betaAdj)]

cpgs_4 <- names(distal_10$colors)[distal_10$colors == 4]
cpgs_4 <- cpgs_4[cpgs_4 %in% rownames(tcga$betaAdj)]

cpgs_5 <- names(distal_10$colors)[distal_10$colors == 5]
cpgs_5 <- cpgs_5[cpgs_5 %in% rownames(tcga$betaAdj)]

cpgs_6 <- names(distal_10$colors)[distal_10$colors == 6]
cpgs_6 <- cpgs_6[cpgs_6 %in% rownames(tcga$betaAdj)]

cpgs_7 <- names(distal_10$colors)[distal_10$colors == 7]
cpgs_7 <- cpgs_7[cpgs_7 %in% rownames(tcga$betaAdj)]

# Combine CpGs from all three cassettes
all_cpgs <- c(cpgs_1, cpgs_2, cpgs_3, cpgs_4, cpgs_5, cpgs_6, cpgs_7)

# Subset the beta matrix
beta_mat <- tcga$betaAdj[all_cpgs, ]

cassette_group <- factor(
  rep(c("1","2","3", "4","5","6", "7"),
      times = c(length(cpgs_1), length(cpgs_2), length(cpgs_3), length(cpgs_4), length(cpgs_5), length(cpgs_6), length(cpgs_7))),
  levels = c("1","2","3", "4","5","6", "7")
)

row_annotation <- rowAnnotation(
  "CpG context" = scanb$annoObj[all_cpgs, "featureClass"],
  "ATAC" = as.factor(tcga$annoObj[all_cpgs, "hasAtacOverlap"]),
  col = list(
    "CpG context" = c("distal"="#FFB3B3", "distal body"="darkred"),
    "ATAC" = c("0"="white", "1"="black")
  )
)

subtypes <- factor(
  tcga$my.pam50.subtype[colnames(beta_mat)],
  levels = c("unclassified","Normal","LumA","LumB","Her2","Basal")
)

top_annotation <- HeatmapAnnotation(
  PAM50 = subtypes,
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

# Define PNG file
png("/Users/isasiain/PhD/Manuscripts/methylation_dynamics/plots_for_figures/fig_3/tcga_heatmap.png",
    width = 6000,   # width in pixels
    height = 8000,  # height in pixels
    res = 800)      # resolution in dpi

# Draw the heatmap
Heatmap(
  beta_mat,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,                
  column_split = subtypes,
  column_title = NULL,
  row_split = cassette_group,   # use cassette group to split rows, not all_cpgs
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2", 
  top_annotation = top_annotation,
  left_annotation = row_annotation,
  name = "Beta"
)

dev.off()


