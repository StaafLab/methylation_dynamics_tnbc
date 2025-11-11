#! usr/bin/Rscript

library(ComplexHeatmap)
library(ggplot2)
library(gprofiler2)

#
# LOADING DATA
#

# Loading dicovery betas
#load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

load("/Volumes/Data/Project_3/validation_cohort/Combined_annotations_rel4SCANB_deNovoMainNMF_distalAtac5000_supervisedSubNMF.RData")
load("/Volumes/Data/Project_3/validation_cohort/FPKM_rel4TNBC_validationCohort_n136.RData")
load("/Volumes/Data/Project_3/validation_cohort/PurBeta_adjustedTumor_betaMatrix_V1_V2_reduced_717459commonCpGs_TNBCs_n136.RData")

# Defining gene-cpg dictionary
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID


# Object name: SCANBrel4_rdata
load("/Users/isasiain/PhD/Projects/project_3/data/SCANBrel4valcohort_annotations.RData") 
rownames(SCANBrel4_rdata) <- SCANBrel4_rdata$id

# Reading cassettes detected in discovery cohort

promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")
proximal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_10.rds")


#
# CONVERT ENSEMBL IDS INTO GENE IDS
#

# Extract Ensembl IDs without versioning (i.e., without ".1", ".2", etc.)
ensembl_ids <- sapply(rownames(gex.data), function(id) {strsplit(id, "\\.")[[1]][1]})

# Convert Ensembl IDs to gene names (e.g., ENTREZGENE)
new_names <- gconvert(ensembl_ids, organism="hsapiens", target="ENTREZGENE", filter_na = F)


combined_names <- sapply(sort(unique(new_names$input_number)), function(num) {

  paste(new_names$name[grep(paste0("^", num, "\\."), new_names$target_number)], collapse = "_")
  
})

# Renaming gex 
rownames(gex.data) <- combined_names

#
# ANALYSING GENES OF INTEREST
#

# Defining genes of interest
gene_ids = c("GBP4", "ZBP1", "OAS2", "CARD16", "SAMD9L")
list_of_heatmaps <- list()

for (current_gene_id in gene_ids) {
  
  
  # Cluster based on methylation
  cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                   annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])
  
  cpgs <- cpgs[names(cpgs) %in% rownames(beta.adjusted)]
  
  cluster_promoter <- kmeans(t(beta.adjusted[names(cpgs),]), centers = 2)
  
  # Determine hypo and hypermethylated cluster
  promoter_state <- if (mean(beta.adjusted[names(cpgs),cluster_promoter$cluster==1]) >
                        mean(beta.adjusted[names(cpgs),cluster_promoter$cluster==2])) {
    
    as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
    
  } else {
    
    as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
    
  }
  
  # Generating annotation for heatmap
  tnbc_annotation <- annotations[,"TNBCtype4"]
  epi_annotation <- annotations[,"NMF_ATAC_finalSubClusters"]
  pam50_annotation <- annotations[,"NCN_PAM50"]
  pam50_annotation <- ifelse(pam50_annotation == "Uncl.", 
                             "Uncl.", 
                             ifelse(pam50_annotation == "Basal", 
                                    "Basal", 
                                    "Non-Basal"))
  
  
  # Generate bottom annotation
  bottom_annotation <- HeatmapAnnotation(
    "FPKM" = anno_barplot(as.numeric(gex.data[current_gene_id,]))
  )
  
  # Generate column annotation
  column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotation,
                                         TNBC = tnbc_annotation,
                                         Epitype = epi_annotation,
                                         col = list(
                                           "PAM50"=c("Basal"="indianred1", "Non-Basal"="darkblue","Uncl."="grey"),
                                           "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey", "UNS"="black"),
                                           "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                       "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
  )
  
  # Generate row annotation
  cpgs_in_cassette <- names(promoter_10$colors)[promoter_10$colors == 10]
  
  row_annotation <- rowAnnotation("CpG_in_cassette" = names(cpgs)[names(cpgs) %in% rownames(beta.adjusted)] %in% cpgs_in_cassette,
                                  col=list("CpG_in_cassette"=c("TRUE" = "black", "FALSE" = "white")))
  
  # Plotting heatmaps
  list_of_heatmaps[[current_gene_id]] <- Heatmap(beta.adjusted[names(cpgs)[names(cpgs) %in% rownames(beta.adjusted)],],
          column_split = promoter_state,
          show_column_names = FALSE,
          show_row_names = FALSE,
          show_row_dend =  FALSE,
          bottom_annotation = bottom_annotation, 
          top_annotation = column_annotation,
          left_annotation = row_annotation,
          clustering_distance_columns =  "euclidean",
          clustering_method_columns = "ward.D2",
          clustering_distance_rows =  "euclidean",
          clustering_method_rows = "ward.D2",
          name = "Beta")
  
  # Testing significance in expression between groups
  print(current_gene_id)
  print(wilcox.test(as.numeric(gex.data[current_gene_id,]) ~ promoter_state))
  
}

list_of_heatmaps[[5]]

#
# PLOTTING TILEPLOT OF PROMOTER HYPERMETHYLATION OF SELECTED GENES
#

gene_ids = c("GBP4", "ZBP1", "OAS2", "CARD16", "SAMD9L")

clusters_methylation <- data.frame(matrix(ncol = ncol(beta.adjusted), nrow = length(gene_ids)))
colnames(clusters_methylation) <- colnames(beta.adjusted)
rownames(clusters_methylation) <- gene_ids

for(current_gene_id in gene_ids) {
  
  # Cluster based on methylation
  cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                   annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])
  
  cpgs <- cpgs[names(cpgs) %in% rownames(beta.adjusted)]
  
  cluster_promoter <- kmeans(t(beta.adjusted[names(cpgs),]), centers = 2)
  
  # Determine hypo and hypermethylated cluster
  promoter_state <- if (mean(beta.adjusted[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==1]) >
                        mean(beta.adjusted[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==2])) {
    
    as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
    
  } else {
    
    as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
    
  }
  
  # Adding data to dataframe
  clusters_methylation[current_gene_id,] <- promoter_state
}

# Convert the character values to numeric for the entire matrix
# "hypermethylated" = 1, "hypomethylated" = 0
numeric_matrix <- apply(clusters_methylation, c(1, 2), function(x) {
  if (x == "Hypermethylated") {
    return(1)
  } else if (x == "Hypomethylated") {
    return(0)
  } else {
    return(NA)
  }
})


# Getting number of hypermethylated samples
agreement_hyper <- colSums(numeric_matrix == 1)

# Defining annotations
top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotation,
                                    TNBC = tnbc_annotation,
                                    Epitype = epi_annotation,
                                    col = list(
                                      "PAM50"=c("Basal"="indianred1", "Non-Basal"="darkblue","Uncl."="grey"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey", "UNS"="white"),
                                      "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                  "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")
                                    )
)

agreement_hyper_annotation <- HeatmapAnnotation("Consenus" = anno_barplot(unname(agreement_hyper)))

# Function to convert pvals to stars
pval_to_stars <- function(pval) {
  if (pval <= 0.0001) return("****")
  else if (pval <= 0.001) return("***")
  else if (pval <= 0.01) return("**")
  else if (pval <= 0.05) return("*")
  else return("ns") 
}

# Basal/nonBasal

pval_stars <- sapply(1:5, function(idx) {
  pval <- chisq.test(clusters_methylation[idx, ], pam50_annotation)$p.value
  pval_to_stars(pval)
})

chi_p_annotation <- rowAnnotation(
  pvalue = anno_text(pval_stars, gp = gpar(fontsize = 10))
)

Heatmap(clusters_methylation,
        col = c("red", "blue"),
        top_annotation = top_annotation,
        left_annotation = chi_p_annotation,
        bottom_annotation = agreement_hyper_annotation,
        show_column_names = FALSE,
        show_row_names = TRUE,
        name = "Methylation",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = as.factor(ifelse(pam50_annotation == "Basal", "Basal", "NonBasal")))

# Lehmann

keep <- !is.na(tnbc_annotation) & tnbc_annotation != "UNS"

pval_stars <- sapply(1:5, function(idx) {
  pval <- chisq.test(clusters_methylation[idx, keep], tnbc_annotation[keep])$p.value
  pval_to_stars(pval)
})


chi_p_annotation <- rowAnnotation(
  pvalue = anno_text(pval_stars, gp = gpar(fontsize = 10))
)

# Setting order
tnbc_annotation <- factor(tnbc_annotation, levels = c("BL1", "BL2", "M", "LAR"))

Heatmap(clusters_methylation[,keep],
        col = c("red", "blue"),
        top_annotation = top_annotation[keep],
        bottom_annotation = agreement_hyper_annotation[keep],
        left_annotation = chi_p_annotation,
        show_column_names = FALSE,
        show_row_names = TRUE,
        name = "Methylation",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = tnbc_annotation[keep])



#
# ANALYSING CASSETTES OF INTEREST. 10
#

# Defining cpgs of interest
cpgs_of_interest <- names(promoter_10$colors)[promoter_10$colors == 10]
cpgs_of_interest <- cpgs_of_interest[cpgs_of_interest %in% rownames(beta.adjusted)]

# Cluster based on methylation
cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% cpgs_of_interest],
                 annoObj$illuminaID[annoObj$illuminaID %in% cpgs_of_interest])

cluster_promoter <- kmeans(t(beta.adjusted[names(cpgs),]), centers = 2)


# Determine hypo and hypermethylated cluster
promoter_state <- if (mean(beta.adjusted[names(cpgs),cluster_promoter$cluster==1]) >
                      mean(beta.adjusted[names(cpgs),cluster_promoter$cluster==2])) {
  
  as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
  
} else {
  
  as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
  
}


# Reorder
promoter_state <- promoter_state[annotations$PD_ID]


column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotation,
                                       TNBC = tnbc_annotation,
                                       Epitype = epi_annotation,
                                       col = list(
                                         "PAM50"=c("Basal"="indianred1", "Non-Basal"="darkblue","Uncl."="grey"),
                                         "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey", "UNS"="black"),
                                         "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                     "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
)


# Plotting heatmap
Heatmap(beta.adjusted[cpgs_of_interest,],
        column_split = promoter_state,
        show_column_names = FALSE,
        show_row_names = FALSE,
        show_row_dend =  FALSE,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2",
        top_annotation = column_annotation,
        use_raster = FALSE,
        name = "Beta")

# Promoter state vs subtypes
pam50_vs_methyl <- as.data.frame(table(pam50_annotation, as.factor(unname(promoter_state))))

ggplot(pam50_vs_methyl, aes(x = Var2, y = Freq, fill = pam50_annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  labs(x = "Promoter state", y = "Count", fill = "PAM50 subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

lehmann_vs_methyl <- as.data.frame(table(tnbc_annotation, as.factor(unname(promoter_state))))

lehmann_vs_methyl <- lehmann_vs_methyl[lehmann_vs_methyl$tnbc_annotation != "UNS",]

ggplot(lehmann_vs_methyl, aes(x = Var2, y = Freq, fill = tnbc_annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  labs(x = "Promoter state", y = "Count", fill = "PAM50 subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Fredlund immue response vs cluster
sample_ids <- sapply(colnames(beta.adjusted), function(id) strsplit(id, ".l")[[1]][1])
immune_scores <- SCANBrel4_rdata[sample_ids, "Fredlund.ImmuneResponse"]
promoter_clusters <- promoter_state[colnames(beta.adjusted)]

df <- data.frame(
  ImmuneScore = immune_scores,
  PromoterCluster = factor(promoter_clusters)
)

# Compute group sizes
df_counts <- df %>%
  group_by(PromoterCluster) %>%
  summarise(n = n(), .groups = "drop")

# Plot with group size annotation
ggplot(df, aes(x = PromoterCluster, y = ImmuneScore)) +
  geom_violin(fill = "black") +
  geom_boxplot(outlier.shape = NA, fill = "gray90", col = "gray60", width = 0.25) +
  stat_compare_means(method = "wilcox.test", label.y = 1050000) +
  geom_text(
    data = df_counts,
    aes(x = PromoterCluster, y = 1150000, label = paste0("n = ", n)),  # position text
    inherit.aes = FALSE
  ) +
  labs(x = "Promoter cluster", y = "Fredlund immune score") +
  theme_bw(base_size = 14)



# Kaplan meier plot of based on methylation state

IDFS <- SCANBrel4_rdata[sapply(colnames(beta.adjusted), function (id) {strsplit(id, ".l")[[1]][1]}), "IDFS"]
IDFS_bin <- SCANBrel4_rdata[sapply(colnames(beta.adjusted), function (id) {strsplit(id, ".l")[[1]][1]}), "IDFSbin"]


chemo <- SCANBrel4_rdata[sapply(colnames(beta.adjusted), function (id) {strsplit(id, ".l")[[1]][1]}), "Chemotherapy"]

combined <- paste(promoter_state, chemo, sep = "_")

# Create the survival object
surv_obj <- Surv(time = IDFS, event = IDFS_bin)

# Fit Kaplan-Meier survival curves stratified by promoter_state
fit <- survfit(surv_obj ~ combined)

# Plot Kaplan-Meier curves
ggsurvplot(
  fit,
  data = data.frame(IDFS = IDFS, IDFS_bin = IDFS_bin, promoter_state = promoter_state),
  pval = TRUE,                 
  conf.int = FALSE,
  risk.table = TRUE,      
  legend.title = "Promoter State",
  xlab = "Time",
  ylab = "Survival Probability (IDFS)",
  palette = "Dark2"             
)

