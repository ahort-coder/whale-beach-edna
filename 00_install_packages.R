# ============================================================
# R Package Setup Script
# Installs every package used across the scripts in "Coding Scripts"
# Generated for Annabelle — 2026-08-17
# ============================================================

# ---- 1. CRAN packages ----
cran_packages <- c(
  "tidyverse",
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "writexl",
  "stringr",
  "tibble",
  "purrr",
  "ggplot2",
  "scales",
  "RColorBrewer",
  "vegan",
  "digest",
  "BiocManager"
)

installed <- rownames(installed.packages())
to_install <- setdiff(cran_packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
} else {
  message("All CRAN packages are already installed.")
}

# ---- 2. Bioconductor packages ----
# These must be installed via BiocManager, not install.packages()
bioc_packages <- c(
  "dada2",
  "phyloseq",
  "Biostrings",
  "ShortRead"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_to_install <- setdiff(bioc_packages, installed)

if (length(bioc_to_install) > 0) {
  BiocManager::install(bioc_to_install, update = FALSE, ask = FALSE)
} else {
  message("All Bioconductor packages are already installed.")
}

# ---- 3. Base R packages ----
# "tools" ships with R automatically — nothing to install.

# ---- 4. Load everything to confirm it all works ----
all_packages <- c(cran_packages, bioc_packages)

invisible(lapply(all_packages, function(pkg) {
  library(pkg, character.only = TRUE)
}))

message("All packages installed and loaded successfully.")
