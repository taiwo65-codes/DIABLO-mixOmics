# ==============================================================================
# Master Pipeline Runner: Run All Multi-Omics Analyses End-to-End
# ==============================================================================

cat("##############################################################################\n")
cat("Starting Complete Multi-Omics DIABLO & Pathway Enrichment Pipeline\n")
cat("##############################################################################\n\n")

# Detect script directory across Rscript CLI, RStudio, and interactive sessions
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(normalizePath(dirname(rstudioapi::getActiveDocumentContext()$path)))
  }
  if (file.exists("00_install_dependencies.R")) {
    return(normalizePath("."))
  } else if (file.exists(file.path("R script", "00_install_dependencies.R"))) {
    return(normalizePath("R script"))
  }
  return(normalizePath("."))
}

script_dir <- get_script_dir()
setwd(script_dir)
cat(sprintf("Working directory set to: %s\n\n", script_dir))

# 1. Check and Install Dependencies
cat("\n>>> Step 0: Checking Environment and Dependencies...\n")
source("00_install_dependencies.R")

# 2. Fit Master 5-Block DIABLO Model & Export Panels A to F
cat("\n>>> Step 1: Running Master 5-Block DIABLO Model & Generating Panels A-F...\n")
source("01_run_diablo_5blocks_master.R")

# 3. Extract Correlation Matrices, Biomarker Tables, and Relevance Networks
cat("\n>>> Step 2: Extracting Biomarker Tables & Relevance Network Matrices...\n")
source("02_extract_biomarkers_and_networks.R")

# 4. Generate Figure 5 (Component 1 Large-Text Loadings)
cat("\n>>> Step 3: Generating Publication Figure 5 (Component 1 Loadings)...\n")
source("03_generate_figure5_loadings_comp1.R")

# 5. Generate Figure 6 (Component 2 Clean Loadings)
cat("\n>>> Step 4: Generating Publication Figure 6 (Component 2 Loadings)...\n")
source("04_generate_figure6_loadings_comp2.R")

cat("\n##############################################################################\n")
cat("All Multi-Omics Pipeline Analyses and Figure Exports Finished Successfully!\n")
cat("##############################################################################\n")
