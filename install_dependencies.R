# Installation script for Shiny AD Time and Space dependencies
# This script installs all required packages for the application

cat("===========================================\n")
cat("Installing dependencies for Shiny AD app\n")
cat("===========================================\n\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", quiet = TRUE)
}

# Bioconductor packages
bioc_packages <- c(
  "ComplexHeatmap",
  "tidySummarizedExperiment"
)

cat("\n--- Installing Bioconductor packages ---\n")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

# CRAN packages
cran_packages <- c(
  "shiny",
  "DT",
  "dplyr",
  "ggpubr",
  "purrr",
  "rstatix",
  "tidyr",
  "readr",
  "knitr",
  "kableExtra",
  "stringr",
  "ggplot2",
  "openxlsx",
  "tibble"
)

cat("\n--- Installing CRAN packages ---\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, quiet = TRUE, dependencies = TRUE)
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

# Optional packages
optional_packages <- c(
  "plotly",
  "lintr",
  "renv"
)

cat("\n--- Installing optional packages ---\n")
for (pkg in optional_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, quiet = TRUE)
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

cat("\n===========================================\n")
cat("Installation complete!\n")
cat("===========================================\n")
cat("\nTo run the application:\n")
cat("  shiny::runApp('App')\n\n")
