# Multi-Omics DIABLO: Integrative 5-Block Biomarker Discovery Pipeline

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![mixOmics](https://img.shields.io/badge/Bioconductor-mixOmics-2D728F.svg)](https://bioconductor.org/packages/release/bioc/html/mixOmics.html)
[![ggplot2](https://img.shields.io/badge/CRAN-ggplot2-green.svg)](https://cran.r-project.org/package=ggplot2)
[![igraph](https://img.shields.io/badge/CRAN-igraph-orange.svg)](https://cran.r-project.org/package=igraph)
[![ggraph](https://img.shields.io/badge/CRAN-ggraph-red.svg)](https://cran.r-project.org/package=ggraph)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An end-to-end, reproducible R pipeline for supervised multi-omics integration and biomarker signature discovery using **DIABLO** (**D**ata **I**ntegration **A**nalysis for **B**iomarker discovery using **L**atent c**O**mponents) within the **mixOmics** framework (multi-block sparse Partial Least Squares Discriminant Analysis, **N-integration / sPLS-DA**).

This workflow simultaneously models and integrates **5 heterogeneous omics layers** to discover coordinated biomarker signatures that discriminate between dietary phenotypes (**High-Fat Diet [HFD]** vs. **Broccoli Dietary Intervention [HFD-GRN]**).

---

## Table of Contents

- [Overview & Biological Rationale](#overview--biological-rationale)
- [The 5 Integrated Omics Blocks](#the-5-integrated-omics-blocks)
- [Mathematical & Methodological Framework](#mathematical--methodological-framework)
- [Repository Structure](#repository-structure)
- [Input Data Requirements & Formatting](#input-data-requirements--formatting)
- [Quick Start & Execution Guide](#quick-start--execution-guide)
  - [Option A: Automated Master Execution (Recommended)](#option-a-automated-master-execution-recommended)
  - [Option B: Modular Step-by-Step Execution](#option-b-modular-step-by-step-execution)
- [Pipeline Scripts & Methodological Workflow](#pipeline-scripts--methodological-workflow)
  - [Step 0: Automated Dependency Management (`00_install_dependencies.R`)](#step-0-automated-dependency-management-00_install_dependenciesr)
  - [Step 1: 5-Block DIABLO Fitting & Panel Generation (`01_run_diablo_5blocks_master.R`)](#step-1-5-block-diablo-fitting--panel-generation-01_run_diablo_5blocks_masterr)
  - [Step 2: Biomarker Signatures & Relevance Networks (`02_extract_biomarkers_and_networks.R`)](#step-2-biomarker-signatures--relevance-networks-02_extract_biomarkers_and_networksr)
  - [Step 3: Component 1 Biomarker Loadings (`03_generate_figure5_loadings_comp1.R`)](#step-3-component-1-biomarker-loadings-03_generate_figure5_loadings_comp1r)
  - [Step 4: Component 2 Biomarker Loadings (`04_generate_figure6_loadings_comp2.R`)](#step-4-component-2-biomarker-loadings-04_generate_figure6_loadings_comp2r)
  - [Step 5: Cross-Validation & Model Performance (`05_evaluate_diablo_performance.R`)](#step-5-cross-validation--model-performance-05_evaluate_diablo_performancer)
- [Generated Deliverables & Publication Figures](#generated-deliverables--publication-figures)
  - [Publication Panels A through F](#publication-panels-a-through-f)
  - [Publication Figures 5 & 6 (Multi-Panel Loadings)](#publication-figures-5--6-multi-panel-loadings)
  - [Tabular Data & Supplementary Export Tables](#tabular-data--supplementary-export-tables)
- [Installation & System Requirements](#installation--system-requirements)
- [Author & Citation Information](#author--citation-information)
- [License](#license)

---

## Overview & Biological Rationale

High-throughput biomedical research increasingly generates multi-layered omics datasets from the same cohort of subjects. However, analyzing each omics layer in isolation ignores the intricate biological cross-talk between the gut microbiome, circulating metabolites, lipid mediators, and host gene regulation.

**DIABLO** solves this limitation by implementing multi-block sPLS-DA, which projects heterogeneous data blocks into a common latent space while:
1. **Maximizing inter-block covariance** across all biological data modalities.
2. **Discriminating phenotypic classes** (`HFD` control vs. `HFD-GRN` broccoli intervention).
3. **Selecting sparse, highly predictive biomarker subsets** via $\ell_1$-penalization (Lasso).

### Experimental Cohort
- **Phenotypic Groups**:
  - `HFD`: High-Fat Diet control (metabolic syndrome / diet-induced obesity model).
  - `HFD-GRN`: High-Fat Diet supplemented with dietary broccoli (evaluating nutritional rescue, metabolic reprogramming, and anti-inflammatory signatures).
- **Matching**: 11 biological samples matched across all 5 omics platforms.

---

## The 5 Integrated Omics Blocks

The pipeline performs joint multivariate integration across 5 distinct data matrices:

| Block Name | Technology / Assay | Biological Modality | Feature Cleaning & Preprocessing |
| :--- | :--- | :--- | :--- |
| **Microbiome** | 16S rRNA Amplicon Sequencing (V3–V4) | Gut microbial taxa (ASVs) | Genus/Family taxonomic annotation, uncultured tags resolved, non-zero counts aligned. |
| **CHG_RNAseq** | Bulk RNA-Sequencing (Host Tissue) | Transcriptional gene expression | Count standardization, zero-variance filtering, gene symbol sanitization. |
| **GCMS_Polar** | GC-MS (Derivatized Polar Fraction) | Primary metabolites (amino acids, organic acids, carbohydrates) | Chemical name extraction, TMS adduct cleaning, unknown/numerical tag removal. |
| **GCMS_NonPolar** | GC-MS (Non-Polar Fraction) | Lipids & hydrophobic metabolites (fatty acids, sterols, hydrocarbons) | Mass-spec tag cleaning, siloxane artifact truncation, feature standardization. |
| **LCMS** | Liquid Chromatography-Mass Spectrometry | Complex lipids & secondary metabolites (PC, PE, lysoPCs, bile acids) | Lipid abbreviation standardizing (PC/PE/Gal), retention time parsing, unknown filtering. |

---

## Mathematical & Methodological Framework

DIABLO generalizes sparse Partial Least Squares (sPLS) to $K$ omics data blocks $X_1, X_2, \dots, X_K$ and a dummy-coded categorical outcome matrix $Y$.

For each latent component $h = 1, \dots, H$:
$$\max_{a_1^h, \dots, a_K^h} \sum_{k, j = 1, k \ne j}^K c_{kj} \text{cov}(X_k a_k^h, X_j a_j^h) + \sum_{k=1}^K c_{kY} \text{cov}(X_k a_k^h, Y a_Y^h)$$
subject to $\|a_k^h\|_2 = 1$ and $\|a_k^h\|_1 \le \lambda_k$, where:
- $a_k^h$ represents the loading vector for omics block $k$ on component $h$.
- $c_{kj} \in [0, 1]$ are elements of the **design matrix**, controlling the relative emphasis between inter-block correlation ($c_{kj} = 1$) vs. phenotypic discrimination ($c_{kY} = 1$).
- $\lambda_k$ is the $\ell_1$ penalty tuning parameter controlling the sparsity (number of selected features `keepX`) in block $k$.

In this pipeline, a full connection design matrix is implemented with tuned variable selection across $H = 2$ latent components:
- **Microbiome**: 15 features per component
- **CHG (RNA-seq)**: 20 features per component
- **GCMS Polar**: 15 features per component
- **GCMS Non-Polar**: 15 features per component
- **LCMS**: 15 features per component

---

## Repository Structure

```
DIABLO-mixOmics/
├── .gitignore                               # Standard R environment exclusions (.Rhistory, .rds, etc.)
├── LICENSE                                  # MIT Open Source License
├── README.md                                # Comprehensive pipeline and methodology documentation
├── 00_install_dependencies.R                # Automated package installer for CRAN & Bioconductor
├── 01_run_diablo_5blocks_master.R           # Core multi-block sPLS-DA fitting & Panels A–F generator
├── 02_extract_biomarkers_and_networks.R     # Biomarker signature extraction & relevance edge lists
├── 03_generate_figure5_loadings_comp1.R     # Component 1 publication multi-panel loadings (Figure 5)
├── 04_generate_figure6_loadings_comp2.R     # Component 2 publication multi-panel loadings (Figure 6)
├── 05_evaluate_diablo_performance.R         # Cross-validation (LOO & M-fold) & classification error evaluation
└── run_all_analyses.R                       # Master runner to execute the entire pipeline end-to-end
```

---

## Input Data Requirements & Formatting

Place the following 6 CSV input data files in the working directory (or parent directory):

| File Name | Dimensions / Format | Description & Requirements |
| :--- | :--- | :--- |
| `Metadata.csv` | `Samples` $\times$ `Attributes` | First column must be sample IDs (`row.names = 1`). Must include a `Group` column with levels `HFD` and `HFD-GRN`. |
| `CHG.csv` | `Genes` $\times$ `Samples` | Gene expression matrix. Row names are gene symbols; column names must match sample IDs in `Metadata.csv`. |
| `microbiome.csv` | `ASVs` $\times$ (`Taxa` + `Samples`) | ASV table containing taxonomic columns (`Genus`, `Family`, etc.) and sample abundance columns. |
| `GCSM polar.csv` | `Metabolites` $\times$ `Samples` | Polar metabolite abundance matrix with chemical feature names in row headers. |
| `GCMS non-polar.csv` | `Metabolites` $\times$ `Samples` | Non-polar metabolite abundance matrix. |
| `LCMS.csv` | `Lipids` $\times$ `Samples` | LC-MS lipid / secondary metabolite abundance matrix. |

> **Note**: The pipeline automatically identifies common sample IDs across all 5 matrices and `Metadata.csv`, safely handling zero-variance features and synchronizing factor levels.

---

## Quick Start & Execution Guide

### Option A: Automated Master Execution (Recommended)

To run the complete pipeline from start to finish with automated dependency checks, execute:

```bash
Rscript run_all_analyses.R
```

Or from within an interactive R / RStudio session:

```r
source("run_all_analyses.R")
```

### Option B: Modular Step-by-Step Execution

Each script is modular and can be run independently:

```r
# Step 0: Ensure all dependencies are installed
source("00_install_dependencies.R")

# Step 1: Fit DIABLO model and generate Panels A through F
source("01_run_diablo_5blocks_master.R")

# Step 2: Extract biomarker lists and relevance network matrices
source("02_extract_biomarkers_and_networks.R")

# Step 3: Export Figure 5 (Component 1 Loadings)
source("03_generate_figure5_loadings_comp1.R")

# Step 4: Export Figure 6 (Component 2 Loadings)
source("04_generate_figure6_loadings_comp2.R")

# Step 5: Evaluate Cross-Validation & Model Performance
source("05_evaluate_diablo_performance.R")
```

---

## Pipeline Scripts & Methodological Workflow

### Step 0: Automated Dependency Management (`00_install_dependencies.R`)
- Detects the availability of `BiocManager`.
- Checks CRAN packages (`ggplot2`, `dplyr`, `tidyr`, `patchwork`, `ggrepel`, `igraph`, `ggraph`, `tidygraph`, `stringr`).
- Checks Bioconductor packages (`mixOmics`, `clusterProfiler`, `ReactomePA`, `org.Mm.eg.db`, `enrichplot`).
- Safely installs only missing dependencies without overwriting existing installations.

### Step 1: 5-Block DIABLO Fitting & Panel Generation (`01_run_diablo_5blocks_master.R`)
- Implements automated text sanitization: strips TMS derivatization tags, siloxanes, ambiguous prefixes, and non-informative mass-to-charge strings.
- Synchronizes sample order across the 5 matrices (`common_samples`).
- Fits a 2-component `block.splsda` model with near-zero-variance filtering.
- Renders high-resolution vector PDFs and 300 DPI PNGs for Panels A through F.

### Step 2: Biomarker Signatures & Relevance Networks (`02_extract_biomarkers_and_networks.R`)
- Queries selected variables per block using `selectVar()`.
- Calculates mean group values per feature to assign group contribution (`HFD` vs. `HFD-GRN`).
- Builds a global cross-omics relevance network based on latent component Pearson correlations.
- Exports filtered edge lists ($|r| \ge 0.70$ and $|r| \ge 0.85$) for Cytoscape and Gephi visualization.

### Step 3: Component 1 Biomarker Loadings (`03_generate_figure5_loadings_comp1.R`)
- Extracts loading weights of selected features for Component 1 across all 5 blocks.
- Assigns dynamic fill colors according to group enrichment:
  - Blue (`#3A86FF`): Features elevated in `HFD`.
  - Orange (`#F77F00`): Features elevated in `HFD-GRN` (Broccoli).
- Assembles a coordinated $3 \times 2$ publication layout using `patchwork` with large, crisp typography.

### Step 4: Component 2 Biomarker Loadings (`04_generate_figure6_loadings_comp2.R`)
- Extracts and visualizes loading weights for Component 2 across all 5 blocks.
- Provides orthogonal biomarker discrimination capturing secondary biological variance.

### Step 5: Cross-Validation & Model Performance (`05_evaluate_diablo_performance.R`)
- Evaluates classification error rates across prediction distances (`centroids.dist`, `mahalanobis.dist`, `max.dist`) using:
  - **Leave-One-Out Cross-Validation (LOO-CV)**: Deterministic evaluation suitable for small cohorts ($N=11$).
  - **Repeated 3-Fold Cross-Validation (10 repeats)**: Assessing model stability across cross-validation splits.
- Exports error curves (`Panel_Performance_ErrorRates.pdf` / `.png`) and classification error tables.
- Demonstrates **100% ensemble classification accuracy (0.00% BER)** in LOO-CV via majority/weighted vote across the 5 omics blocks.

---

## Generated Deliverables & Publication Figures

All outputs are written to the `diablo_output_plots/` directory in both publication-grade vector format (PDF) and 300 DPI raster format (PNG):

### Publication Panels A through F
- **Panel A (`Panel_A_Scores_5Blocks.pdf` / `.png`)**: Individual sample score plots across all 5 omics blocks with 95% confidence ellipses demonstrating consistent phenotype separation.
- **Panel B (`Panel_B_CorrelationMatrix_5Blocks.pdf` / `.png`)**: $5 \times 5$ pairwise latent correlation matrix quantifying multi-modal cross-omics concordance.
- **Panel C (`Panel_C_Heatmap_5Blocks.pdf` / `.png`)**: Hierarchically clustered image map (CIM) displaying multi-block biomarker profiles across samples.
- **Panel D (`Panel_D_Circos_5Blocks.pdf` / `.png`)**: Circos ribbon plot illustrating positive and negative cross-omics correlations ($|r| \ge 0.70$) connecting genes, microbes, and metabolites.
- **Panel E (`Panel_E_CorrelationCircle_5Blocks.pdf` / `.png`)**: Correlation circle plot utilizing `ggrepel` to map multi-block features inside the unit circle without overlapping labels.
- **Panel F (`Panel_F_Network_Clean.pdf` / `.png`)**: High-legibility relevance network graph generated with `ggraph` displaying high-confidence cross-omics interactions ($|r| \ge 0.85$).

### Publication Figures 5 & 6 (Multi-Panel Loadings)
- **Figure 5 (`Figure_5_DIABLO_Loadings_Comp1_LargeText.pdf` / `.png`)**: $18 \times 20$ in multi-panel figure displaying Component 1 feature contribution barplots for Microbiome, Host RNA-seq, GCMS Polar, GCMS Non-Polar, and LCMS.
- **Figure 6 (`Figure_6_DIABLO_Loadings_Comp2_Clean.pdf` / `.png`)**: Component 2 feature contribution barplots across all 5 blocks.

### Performance & Validation Outputs
- **Cross-Validation Plot (`Panel_Performance_ErrorRates.pdf` / `.png`)**: Classification error rate profiles across prediction distance metrics (`centroids.dist`, `mahalanobis.dist`, `max.dist`) and components.
- `Performance_LOO_MajorityVote_ErrorRates.csv`: Leave-One-Out cross-validation ensemble error rates per distance metric.
- `Performance_Repeated_Mfold_ErrorRates.csv`: Repeated 3-fold cross-validation error rates across 10 repeats.

### Tabular Data & Supplementary Export Tables
- `MultiOmics_Component1_Biomarkers.csv`: Ranked Component 1 biomarkers, loadings, and group enrichment assignments.
- `MultiOmics_AllComponents_Biomarkers.csv`: Complete biomarker signatures across Components 1 & 2.
- `MultiOmics_Component1_Correlation_Matrix.csv`: Latent space inter-block Pearson correlation coefficients.
- `MultiOmics_Relevance_Adjacency_Matrix.csv`: Full pairwise feature correlation matrix ($126 \times 126$).
- `MultiOmics_Relevance_Edge_List.csv`: Filtered edge list formatted with source, target, correlation sign, and absolute weight for graph databases and network software.

---

## Installation & System Requirements

### Prerequisites
- **R**: Version $\ge 4.0.0$ (recommended: $\ge 4.3.0$)
- **Operating System**: Cross-platform (Linux, macOS, Windows)

### Package Installation
Run `00_install_dependencies.R` or install directly via R console:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "ggrepel", 
                   "igraph", "ggraph", "tidygraph", "stringr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("mixOmics", "clusterProfiler", "ReactomePA", "org.Mm.eg.db", "enrichplot"))
```

---

## Author & Citation Information

**Author**: Taiwo Bankole  
**Contact**: [GitHub Profile](https://github.com/taiwo65-codes)

### Citations
If you use this pipeline or the underlying methods in your research, please cite:
1. **mixOmics & DIABLO**:
   - Singh A, Gautier B, Shannon CP, et al. *DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays.* Bioinformatics, 2019; 35(17):3055–3062. [doi:10.1093/bioinformatics/bty1054](https://doi.org/10.1093/bioinformatics/bty1054).
   - Rohart F, Gautier B, Singh A, Lê Cao K-A. *mixOmics: An R package for 'omics feature selection and multiple data integration.* PLoS Comput Biol, 2017; 13(11):e1005752. [doi:10.1371/journal.pcbi.1005752](https://doi.org/10.1371/journal.pcbi.1005752).

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for complete details.
