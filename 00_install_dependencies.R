# ==============================================================================
# Script 00: Automated Dependency Installer
# Multi-Omics DIABLO & Pathway Enrichment Analysis
# ==============================================================================

cat("==============================================================================\n")
cat("Checking and Installing Required R Packages...\n")
cat("==============================================================================\n\n")

# 1. Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager from CRAN...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# 2. List of CRAN Packages
cran_packages <- c(
  "ggplot2", "dplyr", "tidyr", "patchwork", "ggrepel", 
  "igraph", "ggraph", "tidygraph", "stringr"
)

# 3. List of Bioconductor Packages
bioc_packages <- c(
  "mixOmics", "clusterProfiler", "ReactomePA", "org.Mm.eg.db", "enrichplot"
)

# Function to check and install
install_if_missing <- function(pkgs, is_bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("--> Installing missing package: %s ...\n", pkg))
      if (is_bioc) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    } else {
      cat(sprintf(" [OK] %s is already installed.\n", pkg))
    }
  }
}

cat("--- Checking CRAN Packages ---\n")
install_if_missing(cran_packages, is_bioc = FALSE)

cat("\n--- Checking Bioconductor Packages ---\n")
install_if_missing(bioc_packages, is_bioc = TRUE)

cat("\n==============================================================================\n")
cat("All required dependencies are installed and ready!\n")
cat("==============================================================================\n")
