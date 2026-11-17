# ==============================================================================
# Script 05: DIABLO Cross-Validation & Model Performance Evaluation
# Blocks: Microbiome, CHG (RNA-seq), GCMS Polar, GCMS Non-Polar, LCMS
# Performs Leave-One-Out (LOO-CV), Repeated M-Fold CV, and Performance Plotting
# ==============================================================================

library(mixOmics)
library(ggplot2)

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
cat("1. Loading Multi-Omics Datasets...\n")
cat("==============================================================================\n")

meta_df       <- read.csv(find_file("Metadata.csv"), row.names = 1, check.names = FALSE)
gcms_np_raw   <- read.csv(find_file("GCMS non-polar.csv"), row.names = 1, check.names = FALSE)
gcms_pol_raw  <- read.csv(find_file("GCSM polar.csv"), row.names = 1, check.names = FALSE)
lcms_raw      <- read.csv(find_file("LCMS.csv"), row.names = 1, check.names = FALSE)
micro_raw     <- read.csv(find_file("microbiome.csv"), row.names = 1, check.names = FALSE)
chg_raw       <- read.csv(find_file("CHG.csv"), row.names = 1, check.names = FALSE)

sample_ids <- rownames(meta_df)

# Microbiome (Taxonomic Genus/Family + SV ID)
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

# Metabolite Cleaning (Extract Names, Drop Unknowns & Pure Numerics, Standardize)
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

gcms_np_mat  <- clean_metabolite_matrix(gcms_np_raw, sample_ids, "GCMS Non-Polar", max_len = 24)
gcms_pol_mat <- clean_metabolite_matrix(gcms_pol_raw, sample_ids, "GCMS Polar", max_len = 24)
lcms_mat     <- clean_metabolite_matrix(lcms_raw, sample_ids, "LCMS", max_len = 24)

# Format CHG (RNA-seq / Transcripts)
chg_sub <- chg_raw[, sample_ids, drop = FALSE]
clean_chg_names <- substr(gsub("[\\(\\)>\\*\"',;: ]", "_", rownames(chg_sub)), 1, 20)
rownames(chg_sub) <- make.unique(clean_chg_names, sep = "")
chg_mat <- t(chg_sub)
mode(chg_mat) <- "numeric"

# Match Samples Across All 5 Omics Blocks
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

# Zero variance check
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
  Microbiome     = c(min(15, ncol(X$Microbiome)), min(15, ncol(X$Microbiome))),
  CHG_RNAseq     = c(min(20, ncol(X$CHG_RNAseq)), min(20, ncol(X$CHG_RNAseq))),
  GCMS_Polar     = c(min(15, ncol(X$GCMS_Polar)), min(15, ncol(X$GCMS_Polar))),
  GCMS_NonPolar  = c(min(15, ncol(X$GCMS_NonPolar)), min(15, ncol(X$GCMS_NonPolar))),
  LCMS           = c(min(15, ncol(X$LCMS)), min(15, ncol(X$LCMS)))
)

cat("\n==============================================================================\n")
cat("2. Fitting 5-Block DIABLO Model...\n")
cat("==============================================================================\n")

diablo_model <- block.splsda(X = X, Y = Y, ncomp = 2, keepX = list.keepX, near.zero.var = TRUE)
cat("DIABLO model fitted successfully!\n\n")

# Setup output directory
output_dir <- if (dir.exists("diablo_output_plots")) "diablo_output_plots" else if (dir.exists("../diablo_output_plots")) "../diablo_output_plots" else "diablo_output_plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

cat("==============================================================================\n")
cat("3. Evaluating Cross-Validation: Leave-One-Out (LOO-CV)...\n")
cat("==============================================================================\n")

set.seed(42)
perf_loo <- perf(
  diablo_model,
  validation = "loo",
  dist = c("max.dist", "centroids.dist", "mahalanobis.dist"),
  progressBar = FALSE
)

cat("\n--- Majority Vote Error Rates (LOO-CV) ---\n")
print(perf_loo$MajorityVote.error.rate)

cat("\n--- Weighted Vote Error Rates (LOO-CV) ---\n")
print(perf_loo$WeightedVote.error.rate)

# Export LOO performance tables
loo_maj_df <- do.call(rbind, lapply(names(perf_loo$MajorityVote.error.rate), function(d) {
  df <- as.data.frame(perf_loo$MajorityVote.error.rate[[d]])
  df$Metric <- rownames(df)
  df$Distance <- d
  return(df)
}))
write.csv(loo_maj_df, file.path(output_dir, "Performance_LOO_MajorityVote_ErrorRates.csv"), row.names = FALSE)

cat("\n==============================================================================\n")
cat("4. Evaluating Repeated 3-Fold Cross-Validation (10 Repeats)...\n")
cat("==============================================================================\n")

set.seed(42)
perf_mfold <- perf(
  diablo_model,
  validation = "Mfold",
  folds = 3,
  nrepeat = 10,
  dist = c("max.dist", "centroids.dist", "mahalanobis.dist"),
  progressBar = FALSE
)

cat("\n--- Majority Vote Error Rates (Repeated 3-Fold CV) ---\n")
print(perf_mfold$MajorityVote.error.rate)

mfold_maj_df <- do.call(rbind, lapply(names(perf_mfold$MajorityVote.error.rate), function(d) {
  df <- as.data.frame(perf_mfold$MajorityVote.error.rate[[d]])
  df$Metric <- rownames(df)
  df$Distance <- d
  return(df)
}))
write.csv(mfold_maj_df, file.path(output_dir, "Performance_Repeated_Mfold_ErrorRates.csv"), row.names = FALSE)

cat("\n==============================================================================\n")
cat("5. Exporting Cross-Validation Performance Plots...\n")
cat("==============================================================================\n")

pdf(file.path(output_dir, "Panel_Performance_ErrorRates.pdf"), width = 10, height = 7)
plot(perf_loo)
dev.off()

png(file.path(output_dir, "Panel_Performance_ErrorRates.png"), width = 10, height = 7, units = "in", res = 300)
plot(perf_loo)
dev.off()

unlink(c("Rplots.pdf", "Rplots1.pdf"))
cat("\n==============================================================================\n")
cat("Script 05 Completed: Cross-Validation & Performance Evaluation Finished!\n")
cat("==============================================================================\n")
