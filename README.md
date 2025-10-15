# Shiny AD in Time and Space

![Code Quality](https://img.shields.io/badge/code%20quality-85%2F100%20(good)-green)
![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Shiny application for interactive visualization of spatial and temporal variation in the skin transcriptome of Atopic Dermatitis (AD).

## 📖 About

This web interface facilitates exploration of data from the study:

**"Assessment of Spatial and Temporal Variation in the Skin Transcriptome of Atopic Dermatitis by Use of 1.5 mm Mini Punch Biopsies"**

*Published in:* [Journal of Investigative Dermatology](https://www.jidonline.org/article/S0022-202X(22)02657-4/fulltext)

## 🚀 Quick Start

### Prerequisites

- R >= 4.0.0
- RStudio (recommended)

### Installation

```r
# Clone the repository
git clone https://github.com/tuhulab/Shiny_AD_time_space.git
cd Shiny_AD_time_space

# Install required packages
source("install_dependencies.R")

# Run the application
shiny::runApp("App")
```

## 📊 Features

- **Time Variation Analysis**: Explore time-specific differentially expressed genes (DEGs)
- **Space Variation (Biological Replicates)**: Analyze variation between biological replicates
- **Cell Type Variation**: View inferred cell composition using MuSiC deconvolution
- **LINC RNA Expression**: Visualize differentially expressed lincRNAs
- **Gene Set Enrichment Analysis (GSEA)**: Browse enrichment results across multiple pathways
- **Variance Partition Results**: Examine variance partitioning for each gene
- **Interactive Tables**: Download and filter analysis results

## 📁 Repository Structure

```
.
├── App/
│   ├── ad_in_time_space.R    # Main Shiny application
│   ├── config.R               # Configuration parameters
│   └── helpers.R              # Helper functions
├── data/                      # Data files (hosted via CDN)
├── DESCRIPTION                # Package dependencies
├── REPRODUCIBILITY.md         # Reproducibility guide
├── check_code_quality.R       # Code quality assessment script
├── prototype_data.R           # Prototype analysis scripts
└── README.md                  # This file
```

## 🔄 Reproducibility

For detailed instructions on reproducing the analysis, see [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

### Package Management

We recommend using `renv` for package version management:

```r
# Install renv
install.packages("renv")

# Restore the project environment
renv::restore()
```

### Data Sources

- **Raw Data**: [GEO GSE193309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309)
- **Analysis Pipeline**: [Zenodo](https://doi.org/10.5281/zenodo.5827799)
- **Processed Data**: Hosted via GitHub CDN (automatically loaded by the app)

## 🧪 Code Quality

This project follows R coding best practices with:

- ✅ Modular code organization (config, helpers, main app)
- ✅ Comprehensive inline documentation
- ✅ Error handling for data loading
- ✅ Consistent naming conventions
- ✅ DRY principle (Don't Repeat Yourself)
- ✅ Reproducible package management

Run code quality checks:

```r
source("check_code_quality.R")
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📧 Contact

Questions regarding the usage of the web application and the data should be addressed to:

**Tu Hu**  
Email: UYHDK AT leo-pharma DOT com

## 📚 Citation

If you use this code or data in your research, please cite:

```
Hu, T. et al. (2023). Assessment of Spatial and Temporal Variation in the 
Skin Transcriptome of Atopic Dermatitis by Use of 1.5 mm Mini Punch Biopsies. 
Journal of Investigative Dermatology.
```

## 🙏 Acknowledgments

- Data analysis pipeline indexed on Zenodo (DOI: 10.5281/zenodo.5827799)
- Single-cell reference data from He et al. 2020 (GSE147424)

## 📈 Version History

- **v1.0.0** (2025-01-15): Initial release with improved code quality and documentation
