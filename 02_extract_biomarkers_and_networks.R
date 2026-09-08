# ==============================================================================
# Script 02: Multi-Omics Biomarker Signatures & Relevance Network Matrices
# ==============================================================================

library(mixOmics)
library(igraph)
library(dplyr)
library(tidyr)

# Helper function to locate input files in current or parent directory
find_file <- function(filename) {
  if (file.exists(filename)) {
    return(filename)
  } else if (file.exists(file.path("..", filename))) {
    return(file.path("..", filename))
  } else {
    stop(sprintf("File '%s' not found in current directory or parent directory.", filename))
  }
}

cat("==============================================================================\n")
cat("1. Loading Datasets and Fitting 5-Block Model...\n")
cat("==============================================================================\n")

meta_df       <- read.csv(find_file("Metadata.csv"), row.names = 1, check.names = FALSE)
gcms_np_raw   <- read.csv(find_file("GCMS non-polar.csv"), row.names = 1, check.names = FALSE)
gcms_pol_raw  <- read.csv(find_file("GCSM polar.csv"), row.names = 1, check.names = FALSE)
lcms_raw      <- read.csv(find_file("LCMS.csv"), row.names = 1, check.names = FALSE)
micro_raw     <- read.csv(find_file("microbiome.csv"), row.names = 1, check.names = FALSE)
chg_raw       <- read.csv(find_file("CHG.csv"), row.names = 1, check.names = FALSE)

sample_ids <- rownames(meta_df)

# Microbiome
taxa_label <- ifelse(!is.na(micro_raw$Genus) & micro_raw$Genus != "" & !grepl("uncultured", micro_raw$Genus, ignore.case = TRUE), 
                     micro_raw$Genus, 
                     ifelse(!is.na(micro_raw$Family) & micro_raw$Family != "" & !grepl("uncultured", micro_raw$Family, ignore.case = TRUE), 
                            micro_raw$Family, "Taxon"))
clean_taxa_label <- gsub("[\\(\\)>\\*\"',;: ]", "", taxa_label)
micro_feature_names <- make.unique(paste0(clean_taxa_label, "_", rownames(micro_raw)), sep = "")
micro_counts <- micro_raw[, sample_ids]
rownames(micro_counts) <- micro_feature_names
micro_mat <- t(micro_counts)
mode(micro_mat) <- "numeric"

# Metabolites
clean_metabolite_matrix <- function(raw_df, samples, block_label = "Block", max_len = 24) {
  df_sub <- raw_df[, samples, drop = FALSE]
  feat_names <- rownames(df_sub)
  
  extracted <- sapply(feat_names, function(nm) {
    if (grepl("^[0-9\\.]+\\([0-9\\.]+,\\s*(.*)\\)$", nm)) {
      sub("^[0-9\\.]+\\([0-9\\.]+,\\s*(.*)\\)$", "\\1", nm)
    } else if (grepl("^[0-9\\.]+\\([0-9\\.]+,\\s*(.*)", nm)) {
      sub("^[0-9\\.]+\\([0-9\\.]+,\\s*(.*)", "\\1", nm)
    } else {
      nm
    }
  }, USE.NAMES = FALSE)
  extracted <- gsub("^[>\\?\\s]+", "", extracted)
  extracted <- gsub("^NF[0-9]+_", "", extracted)

  is_unknown <- grepl("unknown", extracted, ignore.case = TRUE)
  is_num_only <- grepl("^[0-9\\.\\-_: /]+$", extracted)
  is_rt_mz <- grepl("^[0-9]+(\\.[0-9]+)?[-_][0-9]+(\\.[0-9]+)?$", extracted)
  
  valid_idx <- which(!(is_unknown | is_num_only | is_rt_mz))
  df_sub <- df_sub[valid_idx, , drop = FALSE]
  
  cleaned <- extracted[valid_idx]
  cleaned <- gsub(",?O,O,O-TMS", "-TMS", cleaned)
  cleaned <- gsub(",O-TMS", "-TMS", cleaned)
  cleaned <- gsub("-EI-B;-MS", "", cleaned)
  cleaned <- gsub("; EI-B; MS", "", cleaned)
  cleaned <- gsub("-EI-TOF; MS; 3 TMS; BP", "-TMS", cleaned)
  cleaned <- gsub("phosphocholine", "PC", cleaned, ignore.case = TRUE)
  cleaned <- gsub("phosphatidylcholine", "PC", cleaned, ignore.case = TRUE)
  cleaned <- gsub("phosphatidylethanolamine", "PE", cleaned, ignore.case = TRUE)
  cleaned <- gsub("galactopyranoside", "Gal", cleaned, ignore.case = TRUE)
  cleaned <- gsub("octamethyl-", "octamethyl", cleaned)
  cleaned <- gsub("OCTAMETHYLCYCLOTETRASILOXANE", "Octamethylsilox", cleaned)
  cleaned <- gsub("dihydrooxazole", "oxazole", cleaned)
  cleaned <- gsub("[\\(\\)>\\*\"',;: ]", "", cleaned)
  cleaned <- substr(cleaned, 1, max_len)
  
  rownames(df_sub) <- make.unique(cleaned, sep = "")
  mat <- t(df_sub)
  mode(mat) <- "numeric"
  return(mat)
}

gcms_np_mat  <- clean_metabolite_matrix(gcms_np_raw, sample_ids, "GCMS Non-Polar")
gcms_pol_mat <- clean_metabolite_matrix(gcms_pol_raw, sample_ids, "GCMS Polar")
lcms_mat     <- clean_metabolite_matrix(lcms_raw, sample_ids, "LCMS")

# CHG (RNA-seq)
chg_sub <- chg_raw[, sample_ids, drop = FALSE]
clean_chg_names <- substr(gsub("[\\(\\)>\\*\"',;: ]", "_", rownames(chg_sub)), 1, 20)
rownames(chg_sub) <- make.unique(clean_chg_names, sep = "")
chg_mat <- t(chg_sub)
mode(chg_mat) <- "numeric"

# Common Samples
common_samples <- intersect(rownames(meta_df), rownames(gcms_np_mat))
common_samples <- intersect(common_samples, rownames(gcms_pol_mat))
common_samples <- intersect(common_samples, rownames(lcms_mat))
common_samples <- intersect(common_samples, rownames(micro_mat))
common_samples <- intersect(common_samples, rownames(chg_mat))

X <- list(
  Microbiome     = as.matrix(micro_mat[common_samples, , drop = FALSE]),
  CHG_RNAseq     = as.matrix(chg_mat[common_samples, , drop = FALSE]),
  GCMS_Polar     = as.matrix(gcms_pol_mat[common_samples, , drop = FALSE]),
  GCMS_NonPolar  = as.matrix(gcms_np_mat[common_samples, , drop = FALSE]),
  LCMS           = as.matrix(lcms_mat[common_samples, , drop = FALSE])
)

for (b in names(X)) {
  vars <- apply(X[[b]], 2, var, na.rm = TRUE)
  zero_idx <- which(is.na(vars) | vars == 0)
  if (length(zero_idx) > 0) {
    X[[b]] <- X[[b]][, -zero_idx, drop = FALSE]
  }
}

Y <- as.factor(meta_df[common_samples, "Group"])
names(Y) <- common_samples

list.keepX <- list(
  Microbiome     = c(15, 15),
  CHG_RNAseq     = c(20, 20),
  GCMS_Polar     = c(15, 15),
  GCMS_NonPolar  = c(15, 15),
  LCMS           = c(15, 15)
)

diablo_model <- block.splsda(X = X, Y = Y, ncomp = 2, keepX = list.keepX, near.zero.var = TRUE)

cat("==============================================================================\n")
cat("2. Extracting Component 1 Inter-Block Correlation Matrix...\n")
cat("==============================================================================\n")

variates_comp1 <- data.frame(
  Microbiome    = diablo_model$variates$Microbiome[, 1],
  RNAseq        = diablo_model$variates$CHG_RNAseq[, 1],
  GCMS_Polar    = diablo_model$variates$GCMS_Polar[, 1],
  GCMS_NonPolar = diablo_model$variates$GCMS_NonPolar[, 1],
  LCMS          = diablo_model$variates$LCMS[, 1],
  Group         = diablo_model$Y
)

cor_matrix <- cor(variates_comp1[, 1:5])
print(round(cor_matrix, 3))
write.csv(round(cor_matrix, 4), "MultiOmics_Component1_Correlation_Matrix.csv")

cat("\n==============================================================================\n")
cat("3. Extracting Discriminative Multi-Omics Biomarkers (Component 1 & 2)...\n")
cat("==============================================================================\n")

get_comp_biomarkers <- function(model, comp_num = 1) {
  all_biomarkers <- data.frame()
  
  for(b in names(model$X)) {
    loadings_vec <- model$loadings[[b]][, comp_num]
    sel <- names(loadings_vec)[loadings_vec != 0]
    
    if(length(sel) > 0) {
      sub_data <- model$X[[b]][, sel, drop = FALSE]
      medians <- aggregate(sub_data, by = list(Group = model$Y), FUN = median)
      rownames(medians) <- medians$Group
      medians$Group <- NULL
      
      max_grp <- apply(medians, 2, function(x) rownames(medians)[which.max(x)])
      
      df_tmp <- data.frame(
        Block = b,
        Component = comp_num,
        Feature = sel,
        Loading = round(loadings_vec[sel], 4),
        Enriched_In = max_grp,
        stringsAsFactors = FALSE
      )
      all_biomarkers <- rbind(all_biomarkers, df_tmp)
    }
  }
  return(all_biomarkers)
}

biomarker_table_comp1 <- get_comp_biomarkers(diablo_model, comp_num = 1)
biomarker_table_comp2 <- get_comp_biomarkers(diablo_model, comp_num = 2)
biomarker_table_all   <- rbind(biomarker_table_comp1, biomarker_table_comp2)

write.csv(biomarker_table_comp1, "MultiOmics_Component1_Biomarkers.csv", row.names = FALSE)
write.csv(biomarker_table_all, "MultiOmics_AllComponents_Biomarkers.csv", row.names = FALSE)
cat(sprintf(" [OK] Exported %d Component 1 biomarkers to MultiOmics_Component1_Biomarkers.csv\n", nrow(biomarker_table_comp1)))
cat(sprintf(" [OK] Exported %d total biomarkers to MultiOmics_AllComponents_Biomarkers.csv\n", nrow(biomarker_table_all)))

cat("\n==============================================================================\n")
cat("4. Extracting Relevance Network Adjacency Matrix & Edge List (|r| >= 0.80)...\n")
cat("==============================================================================\n")

net_obj <- network(diablo_model, blocks = 1:5, cutoff = 0.80, plot.graph = FALSE)
g <- net_obj$gR
adj_mat <- as.matrix(as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
write.csv(adj_mat, "MultiOmics_Relevance_Adjacency_Matrix.csv")
cat(sprintf(" [OK] Exported %d x %d Adjacency Matrix to MultiOmics_Relevance_Adjacency_Matrix.csv\n", nrow(adj_mat), ncol(adj_mat)))

edge_df <- igraph::as_data_frame(g, what = "edges")
write.csv(edge_df, "MultiOmics_Relevance_Edge_List.csv", row.names = FALSE)
cat(sprintf(" [OK] Exported %d relevance edges to MultiOmics_Relevance_Edge_List.csv\n", nrow(edge_df)))

unlink(c("Rplots.pdf", "Rplots1.pdf", "Rplots2.pdf", "Rplots3.pdf"))
cat("\n==============================================================================\n")
cat("Script 02 Completed Successfully!\n")
cat("==============================================================================\n")
