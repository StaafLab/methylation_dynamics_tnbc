#! usr/bin/Rscript

#
# LIBARIES AND USER DEFINED FUNCTIONS
#
library(httr)
library(KEGGREST)

#Function:  Merging duplicated metabolites from lipids. Taking the mean
merge_duplicate_columns <- function(df) {
  
  cols <- colnames(df)
  unique_cols <- unique(cols)
  
  # Initialize merged_df with the correct number of rows
  merged_df <- data.frame(matrix(nrow = nrow(df), ncol = 0))
  rownames(merged_df) <- rownames(df)
  
  for (name in unique_cols) {
    
    col_indices <- which(cols == name)
    
    if (length(col_indices) == 1) {
      merged_df[, name] <- df[, col_indices]
    } else {
      merged_df[, name] <- rowMeans(df[, col_indices], na.rm = TRUE)
    }
  }
  
  return(merged_df)
}

#
# LOADING DATA
#

lipids <- read.csv2("/Users/isasiain/PhD/Projects/metabolic_subtyping/data/Chinese_cohort_TNBC_PRJNA486023/metabolomics_lipid.csv")

#
# PREPROCESSING LIPIDS
#

# Rotate and using names as colnames
lipids <- t(lipids)
colnames(lipids) <- lipids["Name",]
lipids <- lipids[!rownames(lipids) %in% c("Name", "Peak"), ]

# Convert to numeric
names_of_rows <- rownames(lipids)
lipids <- apply(lipids, 2, function(x) as.numeric(as.character(x)))
rownames(lipids) <- names_of_rows

# Convert to original scale. 2^(value)
lipids <- 2 ^ lipids

# Remove duplicates in lipids. Take the mean
lipids <- merge_duplicate_columns(lipids)

#
# GETTING LIPID CLASSIFICATION
#

# Initialize data frame
lipid_df <- data.frame(
  Name = colnames(lipids),
  Molecule_name = sub(" .*", "", colnames(lipids)),
  Category = NA,
  Family = NA,
  Subfamily = NA,
  stringsAsFactors = FALSE
)


# Glycerophospholipids
lipid_df$Category[grepl("^PC|^PE|^PS|^PI|^PG|^PA|^LPC|^LPE|^LPI|^LPS|^PEtOH|^PMeOH|^Ox", lipid_df$Molecule_name)] <- "Glycerophospholipids"
  lipid_df$Family[grepl("^PC|^OxPC", lipid_df$Molecule_name)] <- "Phosphatidylcholines"
    lipid_df$Subfamily[grepl("^PC", lipid_df$Molecule_name)] <- "Phosphatidylcholines"
    lipid_df$Subfamily[grepl("^OxPC", lipid_df$Molecule_name)] <- "Oxidized phosphatidylcholine"
  
  lipid_df$Family[grepl("^PE|^OxPE", lipid_df$Molecule_name)] <- "Phosphatidylethanolamines"
    lipid_df$Subfamily[grepl("^PE", lipid_df$Molecule_name)] <- "Phosphatidylethanolamine"
    lipid_df$Subfamily[grepl("^OxPE", lipid_df$Molecule_name)] <- "Oxidized Phosphatidylethanolamine"
    
  lipid_df$Family[grepl("^PEtOH", lipid_df$Molecule_name)] <- "Phosphatidylethanol"
    lipid_df$Subfamily[grepl("^PEtOH", lipid_df$Molecule_name)] <- "Phosphatidylethanol"
    
  lipid_df$Family[grepl("^PMeOH", lipid_df$Molecule_name)] <- "Phosphatidylmethanol"
    lipid_df$Subfamily[grepl("^PMeOH", lipid_df$Molecule_name)] <- "Phosphatidylmethanol"
    
  lipid_df$Family[grepl("^PS|^OxPS", lipid_df$Molecule_name)] <- "Phosphatidylserines"
    lipid_df$Subfamily[grepl("^PS", lipid_df$Molecule_name)] <- "Phosphatidylserine"
    lipid_df$Subfamily[grepl("^OxPS", lipid_df$Molecule_name)] <- "Oxidized Phosphatidylserine"
    
  lipid_df$Family[grepl("^PI|^LPI|^OxPI", lipid_df$Molecule_name)] <- "Phosphatidylinositols"
    lipid_df$Subfamily[grepl("^PI", lipid_df$Molecule_name)] <- "Phosphatidylinositol"
    lipid_df$Subfamily[grepl("^LPI", lipid_df$Molecule_name)] <- "Lysophosphatidylinositol"
    lipid_df$Subfamily[grepl("^OxPI", lipid_df$Molecule_name)] <- "Oxidized Phosphatidylinositol"
    
  lipid_df$Family[grepl("^PG|^LPG|^OxPG", lipid_df$Molecule_name)] <- "Phosphatidylglycerols"
    lipid_df$Subfamily[grepl("^PG", lipid_df$Molecule_name)] <- "Phosphatidylglycerol"
    lipid_df$Subfamily[grepl("^LPG", lipid_df$Molecule_name)] <- "Lysophosphatidylglycerol"
    lipid_df$Subfamily[grepl("^OxPG", lipid_df$Molecule_name)] <- "Oxidized Phosphatidylglycerol"
    
  lipid_df$Family[grepl("^PA", lipid_df$Molecule_name)] <- "Phosphatidic acids"
    lipid_df$Subfamily[grepl("^PA", lipid_df$Molecule_name)] <- "Phosphatidic acid"
    
  lipid_df$Family[grepl("^LPC", lipid_df$Molecule_name)] <- "Lysophosphatidylcholines"
    lipid_df$Subfamily[grepl("^LPC", lipid_df$Molecule_name)] <- "Lysophosphatidylcholine"
    
  lipid_df$Family[grepl("^LPE", lipid_df$Molecule_name)] <- "Lysophosphatidylethanolamines"
    lipid_df$Subfamily[grepl("^LPE", lipid_df$Molecule_name)] <- "Lysophosphatidylethanolamine"
    
  lipid_df$Family[grepl("^LPS", lipid_df$Molecule_name)] <- "Lysophosphatidylserines"
    lipid_df$Subfamily[grepl("^LPS", lipid_df$Molecule_name)] <- "Lysophosphatidylserine"

# Sphingolipids
lipid_df$Category[grepl("^Cer|^HexCer|^SHexCer|^SM|^ACar|^GM3", lipid_df$Molecule_name)] <- "Sphingolipids"
  lipid_df$Family[grepl("^Cer-|^HexCer|^SHexCer", lipid_df$Molecule_name)] <- "Ceramides"
    lipid_df$Subfamily[grepl("^Cer-", lipid_df$Molecule_name)] <- "Ceramide"
    lipid_df$Subfamily[grepl("^HexCer", lipid_df$Molecule_name)] <- "Hexosylceramide"
    lipid_df$Subfamily[grepl("^SHexCer", lipid_df$Molecule_name)] <- "Sulfatide"

  lipid_df$Family[grepl("^SM", lipid_df$Molecule_name)] <- "Phosphosphingolipids"
    lipid_df$Subfamily[grepl("^SM", lipid_df$Molecule_name)] <- "Sphingomyelin"
  
  lipid_df$Family[grepl("^GM3", lipid_df$Molecule_name)] <- "Glycosphingolipids"
    lipid_df$Subfamily[grepl("^GM3", lipid_df$Molecule_name)] <- "Ganglioside"

# Fatty acyls
lipid_df$Category[grepl("^FA|^FAHFA|^ACar", lipid_df$Molecule_name)] <- "Fatty acyls"
    lipid_df$Family[grepl("^FAHFA", lipid_df$Molecule_name)] <- "Fatty acid estolides"
      lipid_df$Subfamily[grepl("^FAHFA", lipid_df$Molecule_name)] <- "Fatty Acyl Ester of Hydroxy Fatty Acid"
    
    lipid_df$Family[grepl("^FA$", lipid_df$Molecule_name)] <- "Fatty Acids"
      lipid_df$Subfamily[grepl("^FA$", lipid_df$Molecule_name)] <- "Fatty Acid"
    
    lipid_df$Family[grepl("^ACar", lipid_df$Molecule_name)] <- "Fatty esters"
      lipid_df$Subfamily[grepl("^ACar", lipid_df$Molecule_name)] <- "Acylcarnitine"

# Sterol lipids
lipid_df$Category[grepl("^CE", lipid_df$Molecule_name)] <- "Sterol lipids"
    lipid_df$Family[grepl("^CE", lipid_df$Molecule_name)] <- "Sterols"
      lipid_df$Subfamily[grepl("^CE", lipid_df$Molecule_name)] <- "Cholesteryl ester"
    
      
    # Glycerolipids
lipid_df$Category[grepl("^DAG|^MAG|^TAG|^DGTS|^MGDG|^AcylGlcADG|^BMP|^HBMP|^SQDG|^GlcADG", lipid_df$Molecule_name)] <- "Glycerolipids"
    lipid_df$Family[grepl("^MAG|^DAG|^TAG", lipid_df$Molecule_name)] <- "Neutral glycerolipids"
      lipid_df$Subfamily[grepl("^MAG", lipid_df$Molecule_name)] <- "Monoacylglycerol"
      lipid_df$Subfamily[grepl("^DAG", lipid_df$Molecule_name)] <- "Diacylglycerol"
      lipid_df$Subfamily[grepl("^TAG", lipid_df$Molecule_name)] <- "Triacylglycerol"
    
 
    lipid_df$Family[grepl("^MGDG|^SQDG|^GlcADG|^AcylGlcADG", lipid_df$Molecule_name)] <- "Glycosyldiacylglycerols"
      lipid_df$Subfamily[grepl("^MGDG", lipid_df$Molecule_name)] <- "Monogalactosyldiacylglycerol"
      lipid_df$Subfamily[grepl("^SQDG", lipid_df$Molecule_name)] <- "Sulfoquinovosyl diacylglycerol"
      lipid_df$Subfamily[grepl("^GlcADG", lipid_df$Molecule_name)] <- "Glucosyl-diacylglycerol"
      lipid_df$Subfamily[grepl("^AcylGlcADG", lipid_df$Molecule_name)] <- "Acyl-Glucosyl-diacylglycerol"
    
    
    lipid_df$Family[grepl("^BMP|^HBMP", lipid_df$Molecule_name)] <- "Bis(monoacylglycero)phosphate"
      lipid_df$Subfamily[grepl("^BMP", lipid_df$Molecule_name)] <- "Bis(monoacylglycero)phosphate"
      lipid_df$Subfamily[grepl("^HBMP", lipid_df$Molecule_name)] <- "Hydroxilated bis(monoacylglycero)phosphate"
    
    lipid_df$Family[grepl("^DGTS", lipid_df$Molecule_name)] <- "Betaine lipids"
      lipid_df$Subfamily[grepl("^DGTS", lipid_df$Molecule_name)] <- "Diacylglyceryl-N,N,N-trimethylhomoserine"




# View result
rownames(lipid_df) <- lipid_df$Name
lipid_df[lipid_df$Subfamily == "Cholesteryl ester",]



query <- "2',4'-Dihydroxyacetophenone"
url <- paste0("https://rest.kegg.jp/find/compound/", URLencode(query))
res <- GET(url)
content(res, "text")

#
# PREPROCESSING POLAR
#

# Getting named vector of all compounds annotated in KEGG
all_compounds <- keggList("compound")

# Split each entry by "; " and keep the KEGG ID
split_compounds <- unlist(
  lapply(seq_along(all_compounds), function(i) {
    id <- names(all_compounds)[i]
    names_split <- strsplit(all_compounds[i], ";\\s*")[[1]]
    # assign the same KEGG ID to each split name
    setNames(names_split, rep(id, length(names_split)))
  })
)

# Locating polar metabolites
metabolites <- polar$Name[1:40]

# 1. Remove anything in parentheses
clean_names <- sapply(metabolites, function(x) {
  # check if it contains parentheses
  if (grepl("\\(", x)) {
    # extract the text inside parentheses
    sub(".*\\(([^)]*)\\).*", "\\1", x)
  } else {
    # otherwise keep as is
    x
  }
})

# 2. Replace apostrophes and other tricky characters with nothing
clean_names <- gsub("['`]", "", clean_names)

# 3. Trim leading/trailing whitespace
clean_names <- trimws(clean_names)

# 4. Optionally, replace multiple spaces with single space
clean_names <- gsub("\\s+", " ", clean_names)

clean_names

test_matches <- lapply(clean_names, function(metabolite) {
  
  # Step 1: exact match
  exact_matches <- sapply(split_compounds, function(x) metabolite %in% x)
  
  if (any(exact_matches)) {
    matches <- split_compounds[exact_matches]
    
  } else {
    matches <- split_compounds[grep(metabolite, unname(split_compounds))]
  }
  
  return(matches)
  
})


test_matches

