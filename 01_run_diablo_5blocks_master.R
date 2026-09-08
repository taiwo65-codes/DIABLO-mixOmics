# ==============================================================================
# Script 01: Master Multi-Omics 5-Block DIABLO Analysis & Figure Generation
# Blocks: Microbiome, CHG (RNA-seq), GCMS Polar, GCMS Non-Polar, LCMS
# ==============================================================================

library(mixOmics)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(ggrepel)

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
  message(sprintf(" [%s] Retained %d of %d annotated features.", block_label, ncol(mat), nrow(raw_df)))
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
message(sprintf(" [CHG / RNAseq] Retained %d of %d genes.", ncol(chg_mat), nrow(chg_raw)))

# ------------------------------------------------------------------
# 2. Match Samples Across All 5 Omics Blocks
# ------------------------------------------------------------------
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
    message(sprintf(" Removing %d zero-variance variables from block %s", length(zero_idx), b))
    X[[b]] <- X[[b]][, -zero_idx, drop = FALSE]
  }
}

Y <- as.factor(meta_df[common_samples, "Group"])
names(Y) <- common_samples

cat("\n==============================================================================\n")
cat("2. Fitting 5-Block DIABLO (Multi-Block sPLS-DA) Model...\n")
cat("==============================================================================\n")

list.keepX <- list(
  Microbiome     = c(min(15, ncol(X$Microbiome)), min(15, ncol(X$Microbiome))),
  CHG_RNAseq     = c(min(20, ncol(X$CHG_RNAseq)), min(20, ncol(X$CHG_RNAseq))),
  GCMS_Polar     = c(min(15, ncol(X$GCMS_Polar)), min(15, ncol(X$GCMS_Polar))),
  GCMS_NonPolar  = c(min(15, ncol(X$GCMS_NonPolar)), min(15, ncol(X$GCMS_NonPolar))),
  LCMS           = c(min(15, ncol(X$LCMS)), min(15, ncol(X$LCMS)))
)

diablo_model <- block.splsda(X = X, Y = Y, ncomp = 2, keepX = list.keepX, near.zero.var = TRUE)
cat("DIABLO model fitted successfully!\n\n")

# Setup output directory
output_dir <- if (dir.exists("diablo_output_plots")) "diablo_output_plots" else if (dir.exists("../diablo_output_plots")) "../diablo_output_plots" else "diablo_output_plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 5-Block Color Palette & Group Colors
block_colors <- c(
  "Microbiome"    = "#F4A261",
  "CHG_RNAseq"    = "#6B5B95",
  "GCMS_Polar"    = "#2A9D8F",
  "GCMS_NonPolar" = "#264653",
  "LCMS"          = "#E76F51"
)
group_colors <- c("HFD" = "#3A86FF", "HFD-GRN" = "#F77F00")

cat("==============================================================================\n")
cat("3. Exporting Publication Panels (A through F)...\n")
cat("==============================================================================\n")

# Panel A: Block Scores
cat("--> Exporting Panel A: Sample Scores (5 Blocks)...\n")
pdf(file.path(output_dir, "Panel_A_Scores_5Blocks.pdf"), width = 12, height = 8)
plotIndiv(diablo_model, legend = TRUE, col = group_colors[levels(Y)],
          title = "A: Block Scores", ellipse = TRUE)
dev.off()

png(file.path(output_dir, "Panel_A_Scores_5Blocks.png"), width = 12, height = 8, units = "in", res = 300)
plotIndiv(diablo_model, legend = TRUE, col = group_colors[levels(Y)],
          title = "A: Block Scores", ellipse = TRUE)
dev.off()

# Panel B: Correlation Matrix
cat("--> Exporting Panel B: 5x5 Inter-Block Correlation Matrix...\n")
pdf(file.path(output_dir, "Panel_B_CorrelationMatrix_5Blocks.pdf"), width = 9, height = 9)
plotDiablo(diablo_model, ncomp = 1, col.per.group = group_colors[levels(Y)])
dev.off()

png(file.path(output_dir, "Panel_B_CorrelationMatrix_5Blocks.png"), width = 9, height = 9, units = "in", res = 300)
plotDiablo(diablo_model, ncomp = 1, col.per.group = group_colors[levels(Y)])
dev.off()

# Panel C: Clustered Heatmap
cat("--> Exporting Panel C: 5-Block Clustered Heatmap (CIM)...\n")
pdf(file.path(output_dir, "Panel_C_Heatmap_5Blocks.pdf"), width = 12, height = 10)
cimDiablo(
  diablo_model,
  color.Y = group_colors[levels(Y)],
  color.blocks = unname(block_colors),
  margins = c(12, 10),
  row.cex = 0.85,
  col.cex = 0.65,
  size.legend = 0.8,
  legend.position = "topright"
)
dev.off()

png(file.path(output_dir, "Panel_C_Heatmap_5Blocks.png"), width = 12, height = 10, units = "in", res = 300)
cimDiablo(
  diablo_model,
  color.Y = group_colors[levels(Y)],
  color.blocks = unname(block_colors),
  margins = c(12, 10),
  row.cex = 0.85,
  col.cex = 0.65,
  size.legend = 0.8,
  legend.position = "topright"
)
dev.off()

# Panel D: Circos Plot
cat("--> Exporting Panel D: 5-Block Circos Ribbon Plot...\n")
sel_vars_c1 <- selectVar(diablo_model, comp = 1)
X_circos <- list(
  Microbiome     = X$Microbiome[, sel_vars_c1$Microbiome$name, drop = FALSE],
  CHG_RNAseq     = X$CHG_RNAseq[, sel_vars_c1$CHG_RNAseq$name, drop = FALSE],
  GCMS_Polar     = X$GCMS_Polar[, sel_vars_c1$GCMS_Polar$name, drop = FALSE],
  GCMS_NonPolar  = X$GCMS_NonPolar[, sel_vars_c1$GCMS_NonPolar$name, drop = FALSE],
  LCMS           = X$LCMS[, sel_vars_c1$LCMS$name, drop = FALSE]
)
diablo_circos <- block.splsda(X = X_circos, Y = Y, ncomp = 1)

pdf(file.path(output_dir, "Panel_D_Circos_5Blocks.pdf"), width = 14, height = 14)
circosPlot(
  diablo_circos, 
  comp = 1,
  cutoff = 0.70, 
  color.blocks = unname(block_colors), 
  color.cor = c("#1A1C1E", "#F4D03F"), 
  color.Y = c("#3A86FF", "#F77F00"), 
  line = FALSE, 
  size.variables = 0.85, 
  size.labels = 1.5,
  size.legend = 1.1,
  var.adj = 0.5
)
dev.off()

png(file.path(output_dir, "Panel_D_Circos_5Blocks.png"), width = 14, height = 14, units = "in", res = 300)
circosPlot(
  diablo_circos, 
  comp = 1,
  cutoff = 0.70, 
  color.blocks = unname(block_colors), 
  color.cor = c("#1A1C1E", "#F4D03F"), 
  color.Y = c("#3A86FF", "#F77F00"), 
  line = FALSE, 
  size.variables = 0.85, 
  size.labels = 1.5,
  size.legend = 1.1,
  var.adj = 0.5
)
dev.off()

# Panel E: Repelled Correlation Circle Plot
cat("--> Exporting Panel E: Repelled Correlation Circle Plot...\n")
df_coords_raw <- plotVar(diablo_model, plot = FALSE)
df_coords <- data.frame(
  x = df_coords_raw$x,
  y = df_coords_raw$y,
  Feature = df_coords_raw$names,
  Block = df_coords_raw$Block,
  stringsAsFactors = FALSE
)
df_coords$r <- sqrt(df_coords$x^2 + df_coords$y^2)
df_clean <- subset(df_coords, r >= 0.65)
df_clean$Feature <- substr(gsub("[\\(\\)>\\*_\",;: ]", "", df_clean$Feature), 1, 14)

circle_fun <- function(r = 1, npoints = 200) {
  tt <- seq(0, 2 * pi, length.out = npoints)
  data.frame(x = r * cos(tt), y = r * sin(tt))
}
circle_inner <- circle_fun(r = 0.5)
circle_outer <- circle_fun(r = 1.0)

p_var <- ggplot() +
  geom_path(data = circle_outer, aes(x = x, y = y), color = "#1A1C1E", linewidth = 0.8) +
  geom_path(data = circle_inner, aes(x = x, y = y), color = "#6C7278", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "#6C7278") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "#6C7278") +
  geom_point(data = df_clean, aes(x = x, y = y, color = Block), size = 2.8, alpha = 0.9) +
  geom_text_repel(
    data = df_clean,
    aes(x = x, y = y, label = Feature, color = Block),
    size = 3.2,
    fontface = "bold",
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "#6C7278",
    segment.size = 0.35,
    show.legend = FALSE
  ) +
  scale_color_manual(values = block_colors) +
  scale_x_continuous(limits = c(-1.15, 1.15), breaks = seq(-1, 1, by = 0.5)) +
  scale_y_continuous(limits = c(-1.15, 1.15), breaks = seq(-1, 1, by = 0.5)) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Correlation Circle Plot (5-Block Multi-Omics)",
    subtitle = "Selected multi-omics variables (r >= 0.65)",
    x = "Component 1",
    y = "Component 2",
    color = "Omics Block"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#6C7278", margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "#1A1C1E", size = 10)
  )

ggsave(file.path(output_dir, "Panel_E_CorrelationCircle_5Blocks.pdf"), plot = p_var, width = 9, height = 8, dpi = 300)
ggsave(file.path(output_dir, "Panel_E_CorrelationCircle_5Blocks.png"), plot = p_var, width = 9, height = 8, dpi = 300)

# Panel F: Clean 5-Block ggraph Relevance Network
cat("--> Exporting Panel F: Clean Relevance Network (|r| >= 0.85)...\n")
net <- network(diablo_model, blocks = 1:5, cutoff = 0.85, plot.graph = FALSE)
g <- net$gR

node_raw_names <- V(g)$name
node_blocks <- ifelse(grepl("_GCMS_NonPolar$", node_raw_names), "GCMS_NonPolar",
               ifelse(grepl("_GCMS_Polar$", node_raw_names), "GCMS_Polar",
               ifelse(grepl("_LCMS$", node_raw_names), "LCMS",
               ifelse(grepl("_Microbiome$", node_raw_names), "Microbiome",
               ifelse(grepl("_CHG_RNAseq$", node_raw_names), "CHG_RNAseq", "Other")))))

for(i in seq_along(node_raw_names)) {
  if(node_blocks[i] == "Other") {
    nm <- node_raw_names[i]
    if (nm %in% colnames(X$GCMS_NonPolar)) node_blocks[i] <- "GCMS_NonPolar"
    else if (nm %in% colnames(X$GCMS_Polar)) node_blocks[i] <- "GCMS_Polar"
    else if (nm %in% colnames(X$LCMS)) node_blocks[i] <- "LCMS"
    else if (nm %in% colnames(X$Microbiome)) node_blocks[i] <- "Microbiome"
    else if (nm %in% colnames(X$CHG_RNAseq)) node_blocks[i] <- "CHG_RNAseq"
  }
}

clean_labels <- gsub("_(GCMS_NonPolar|GCMS_Polar|LCMS|Microbiome|CHG_RNAseq)$", "", node_raw_names)
V(g)$Block <- factor(node_blocks, levels = names(block_colors))
V(g)$label <- clean_labels

tbl_g <- as_tbl_graph(g)

set.seed(42)
p_net <- ggraph(tbl_g, layout = "graphopt", charge = 0.05, spring.length = 50) +
  geom_edge_link(aes(edge_alpha = abs(weight)), color = "#6C7278", width = 0.8, show.legend = FALSE) +
  geom_node_point(aes(fill = Block), shape = 21, size = 11, color = "#1A1C1E", stroke = 0.9) +
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = 3.2,
    fontface = "bold",
    color = "#1A1C1E",
    box.padding = 0.35,
    point.padding = 0.5,
    max.overlaps = 50,
    segment.color = "#6C7278",
    segment.size = 0.3
  ) +
  scale_fill_manual(
    name = "Omics Block",
    values = block_colors,
    drop = FALSE
  ) +
  theme_void() +
  labs(title = "5-Block Multi-Omics Relevance Network (|r| >= 0.85)") +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 15)),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(file.path(output_dir, "Panel_F_Network_Clean.pdf"), plot = p_net, width = 11, height = 11, dpi = 300)
ggsave(file.path(output_dir, "Panel_F_Network_Clean.png"), plot = p_net, width = 11, height = 11, dpi = 300)

unlink(c("Rplots.pdf", "Rplots1.pdf", "Rplots2.pdf", "Rplots3.pdf"))
cat("\n==============================================================================\n")
cat("Script 01 Completed: All 5-Block DIABLO publication panels exported!\n")
cat("==============================================================================\n")
