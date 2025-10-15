# Reproducibility Guide

## Overview
This document provides instructions for reproducing the analysis and running the Shiny application for "Assessment of Spatial and Temporal Variation in the Skin Transcriptome of Atopic Dermatitis".

## System Requirements

### R Version
- R >= 4.0.0 is required
- RStudio is recommended for interactive development

### Required R Packages
All required packages are listed in the `DESCRIPTION` file. Install them using:

```r
# Install Bioconductor packages first
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "ComplexHeatmap",
    "tidySummarizedExperiment"
))

# Install CRAN packages
install.packages(c(
    "shiny", "DT", "dplyr", "ggpubr", "purrr", 
    "rstatix", "tidyr", "readr", "knitr", 
    "kableExtra", "stringr", "ggplot2", "openxlsx", 
    "tibble", "plotly"
))
```

### Using renv for Reproducibility (Recommended)
For exact package version reproducibility:

```r
# Install renv
install.packages("renv")

# Initialize renv (only needed once)
renv::init()

# This will create renv.lock with all package versions
# To restore the exact environment later:
renv::restore()
```

## Data Sources

### Raw Data
- **GEO Repository**: GSE193309
- **URL**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309

### Processed Data
All processed data files are hosted on GitHub via jsDelivr CDN:
- Time variation data
- Space variation data (biological replicate)
- Cell type deconvolution data
- LINC RNA expression data
- GSEA results
- Variance Partition results
- And more...

Files are automatically loaded from the CDN when running the application.

## Running the Application

### Local Deployment
1. Clone the repository:
```bash
git clone https://github.com/tuhulab/Shiny_AD_time_space.git
cd Shiny_AD_time_space
```

2. Open R or RStudio and set working directory:
```r
setwd("/path/to/Shiny_AD_time_space")
```

3. Run the application:
```r
# Option 1: Run from App directory
shiny::runApp("App")

# Option 2: Source the file directly
source("App/ad_in_time_space.R")
```

### Deploying to shinyapps.io
```r
# Install rsconnect
install.packages("rsconnect")

# Configure your account (first time only)
rsconnect::setAccountInfo(
    name='your-account',
    token='your-token',
    secret='your-secret'
)

# Deploy
rsconnect::deployApp(appDir = "App", appName = "AD_time_space")
```

## Prototype Analysis
The `prototype_data.R` file contains exploratory analysis code used during development. To run:

```r
# Ensure you're in the project root directory
source("prototype_data.R")
```

Note: This requires local data files in the `data/` directory.

## Session Information
To document your R session for reproducibility:

```r
# Get session info
sessionInfo()

# Or for more detailed info
devtools::session_info()
```

## Troubleshooting

### Common Issues

1. **Package Installation Errors**
   - For Bioconductor packages, ensure BiocManager is up to date
   - Check that your R version is compatible

2. **Data Loading Errors**
   - Verify internet connection for CDN access
   - Check that URLs in the code are accessible

3. **Memory Issues**
   - Some visualizations may require significant memory
   - Increase R memory limit if needed: `memory.limit(size=8000)` (Windows)

## Citation
When using this code or data, please cite:
- Hu et al., "Assessment of Spatial and Temporal Variation in the Skin Transcriptome of Atopic Dermatitis by Use of 1.5 mm Mini Punch Biopsies", Journal of Investigative Dermatology

## Data Analysis Pipeline
The complete reproducible data analysis pipeline is indexed on Zenodo:
- DOI: https://doi.org/10.5281/zenodo.5827799

## Contact
Questions regarding reproducibility should be addressed to Tu Hu (UYHDK AT leo-pharma DOT com)
