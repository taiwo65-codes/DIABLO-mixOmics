# ==============================================================================
# Script 03: Figure 5 - Component 1 Large-Text DIABLO Loadings Multi-Panel
# Blocks: Microbiome, CHG (RNA-seq), GCMS Polar, GCMS Non-Polar, LCMS
# ==============================================================================

library(mixOmics)
library(ggplot2)
library(dplyr)
library(patchwork)

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
cat("1. Loading Datasets and Fitting Model for Figure 5 Loadings...\n")
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
  Microbiome     = c(min(15, ncol(X$Microbiome)), min(15, ncol(X$Microbiome))),
  CHG_RNAseq     = c(min(20, ncol(X$CHG_RNAseq)), min(20, ncol(X$CHG_RNAseq))),
  GCMS_Polar     = c(min(15, ncol(X$GCMS_Polar)), min(15, ncol(X$GCMS_Polar))),
  GCMS_NonPolar  = c(min(15, ncol(X$GCMS_NonPolar)), min(15, ncol(X$GCMS_NonPolar))),
  LCMS           = c(min(15, ncol(X$LCMS)), min(15, ncol(X$LCMS)))
)

diablo_model <- block.splsda(X = X, Y = Y, ncomp = 2, keepX = list.keepX, near.zero.var = TRUE)

cat("==============================================================================\n")
cat("2. Extracting Component 1 Loadings and Median Group Associations...\n")
cat("==============================================================================\n")

extract_diablo_loadings <- function(model, block_name, comp_num = 1) {
  loadings_vec <- model$loadings[[block_name]][, comp_num]
  selected_feats <- names(loadings_vec)[loadings_vec != 0]
  if (length(selected_feats) == 0) return(NULL)
  
  data_mat <- model$X[[block_name]][, selected_feats, drop = FALSE]
  groups <- model$Y
  group_medians <- aggregate(data_mat, by = list(Group = groups), FUN = median)
  rownames(group_medians) <- group_medians$Group
  group_medians$Group <- NULL
  
  max_group <- apply(group_medians, 2, function(x) rownames(group_medians)[which.max(x)])
  
  df <- data.frame(
    Feature = selected_feats,
    Importance = loadings_vec[selected_feats],
    Group = max_group,
    Block = block_name,
    stringsAsFactors = FALSE
  )
  df <- df[order(df$Importance), ]
  df$Feature <- factor(df$Feature, levels = df$Feature)
  return(df)
}

plot_single_block_loading <- function(df, title_name) {
  ggplot(df, aes(x = Importance, y = Feature, fill = Group)) +
    geom_col(width = 0.72, color = "#1A1C1E", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "#1A1C1E", linewidth = 0.6) +
    scale_fill_manual(
      values = c("HFD" = "#3A86FF", "HFD-GRN" = "#F77F00"),
      name = "Outcome",
      drop = FALSE
    ) +
    labs(
      title = paste0("Contribution on comp 1\nBlock '", title_name, "'"),
      x = "Loading Value",
      y = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 11, face = "bold", color = "#1A1C1E", hjust = 1),
      axis.text.x = element_text(size = 10, color = "#1A1C1E"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, lineheight = 1.2),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(t = 15, r = 20, b = 15, l = 25)
    )
}

df_micro <- extract_diablo_loadings(diablo_model, "Microbiome", 1)
df_chg   <- extract_diablo_loadings(diablo_model, "CHG_RNAseq", 1)
df_pol   <- extract_diablo_loadings(diablo_model, "GCMS_Polar", 1)
df_np    <- extract_diablo_loadings(diablo_model, "GCMS_NonPolar", 1)
df_lcms  <- extract_diablo_loadings(diablo_model, "LCMS", 1)

p1 <- plot_single_block_loading(df_micro, "Microbiome")
p2 <- plot_single_block_loading(df_chg, "RNAseq")
p3 <- plot_single_block_loading(df_pol, "GCMS_Polar")
p4 <- plot_single_block_loading(df_np, "GCMS_NonPolar")
p5 <- plot_single_block_loading(df_lcms, "LCMS")

cat("==============================================================================\n")
cat("3. Assembling and Exporting Figure 5 (Large Text)...\n")
cat("==============================================================================\n")

combined_figure <- (p1 + p2) / (p3 + p4) / (p5 + plot_spacer()) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    plot.margin = margin(15, 20, 15, 20)
  )

output_dir <- if (dir.exists("diablo_output_plots")) "diablo_output_plots" else if (dir.exists("../diablo_output_plots")) "../diablo_output_plots" else "diablo_output_plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

ggsave(
  filename = file.path(output_dir, "Figure_5_DIABLO_Loadings_Comp1_LargeText.pdf"),
  plot = combined_figure,
  width = 18,
  height = 20,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Figure_5_DIABLO_Loadings_Comp1_LargeText.png"),
  plot = combined_figure,
  width = 18,
  height = 20,
  dpi = 300,
  units = "in",
  limitsize = FALSE
)

cat(" [OK] Exported Figure_5_DIABLO_Loadings_Comp1_LargeText.pdf (18x20 in)\n")
cat(" [OK] Exported Figure_5_DIABLO_Loadings_Comp1_LargeText.png (300 DPI)\n")

unlink(c("Rplots.pdf", "Rplots1.pdf"))
cat("\n==============================================================================\n")
cat("Script 03 Completed Successfully!\n")
cat("==============================================================================\n")
