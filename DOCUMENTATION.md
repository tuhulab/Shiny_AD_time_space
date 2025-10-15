# Project Documentation

## Overview

This repository contains a Shiny web application for interactive exploration of Atopic Dermatitis (AD) transcriptome data, focusing on spatial and temporal variation.

## Directory Structure

```
Shiny_AD_time_space/
│
├── App/                           # Shiny application files
│   ├── ad_in_time_space.R        # Main application (UI + Server)
│   ├── config.R                   # Configuration and constants
│   └── helpers.R                  # Reusable helper functions
│
├── data/                          # Data files (hosted on CDN)
│   ├── time_variation_t.rds      # Time variation analysis results
│   ├── bioreplicate_t.rds        # Biological replicate data
│   ├── cell_type_deconvolution.rds # Cell type composition
│   ├── gsea_res.rds              # GSEA analysis results
│   ├── heatmap_linc.rds          # LINC RNA heatmap data
│   ├── table_s*.csv              # Supplementary tables
│   └── ... (other data files)
│
├── CHANGELOG.md                   # Version history and changes
├── CODE_STYLE.md                  # Coding conventions
├── CONTRIBUTING.md                # Contribution guidelines
├── DESCRIPTION                    # R package dependencies
├── LICENSE                        # MIT License
├── README.md                      # Main documentation
├── REPRODUCIBILITY.md             # Reproducibility instructions
│
├── check_code_quality.R          # Code quality assessment tool
├── install_dependencies.R        # Dependency installation script
├── prototype_data.R              # Exploratory analysis scripts
└── test_app.R                    # Testing and validation script
```

## Key Features

### 1. Time Variation Analysis
- Visualize gene expression changes over time
- Compare absolute time (quarters) vs relative time (visits)
- Filter by skin type (LS/NL/HC)

### 2. Space Variation Analysis
#### Biological Replicates
- Compare paired biological replicates
- Statistical testing (paired t-tests)
- Interactive gene selection

#### Anatomic Regions
- Variation across different anatomic locations
- Comprehensive statistical tables

### 3. Cell Type Deconvolution
- MuSiC-based cell type estimation
- Single-cell reference: GSE147424 (He et al. 2020)
- Interactive subject selection

### 4. Gene Set Enrichment Analysis (GSEA)
- Multiple pathway databases:
  - Reactome
  - Transcription Factor Targets
  - GO: Biological Process
  - GO: Molecular Function
- View by NES or adjusted p-value

### 5. Additional Analyses
- LINC RNA expression heatmap
- AD signature genes
- Variance partition results
- Tissue injury markers
- Subcutis markers

## Application Architecture

### Configuration (`App/config.R`)
Centralizes all configuration parameters:
- Data URLs (CDN endpoints)
- Color palettes
- Plot defaults
- UI text constants

Benefits:
- Easy maintenance
- Single source of truth
- Simplified updates

### Helper Functions (`App/helpers.R`)
Reusable utility functions:
- `load_rds_data()` - Safe RDS data loading with error handling
- `load_csv_data()` - Safe CSV data loading with error handling
- `create_datatable()` - Standardized DataTable creation
- `filter_gsea_results()` - GSEA data filtering
- `create_bioreplicate_plot()` - Complex plot generation
- And more...

Benefits:
- Reduced code duplication
- Consistent error handling
- Easier testing and maintenance

### Main Application (`App/ad_in_time_space.R`)
Structure:
1. Header with metadata
2. Library imports
3. Configuration and helper loading
4. UI definition (organized by tabs)
5. Server logic (organized by outputs)
6. App runner

## Data Flow

```
CDN (jsDelivr/GitHub)
        ↓
Helper Functions (load_*_data)
        ↓
Error Handling
        ↓
Server Reactive Context
        ↓
Processing/Visualization
        ↓
UI Output (plots/tables)
```

## Reproducibility Strategy

### 1. Package Management
- `DESCRIPTION` file lists all dependencies
- `install_dependencies.R` automates installation
- Optional `renv` for exact version locking

### 2. Data Provenance
- Raw data: GEO GSE193309
- Analysis pipeline: Zenodo DOI 10.5281/zenodo.5827799
- Processed data: GitHub CDN (versioned)

### 3. Documentation
- Comprehensive README
- Detailed REPRODUCIBILITY.md
- Inline code comments
- Function documentation

### 4. Version Control
- Git for source control
- Tagged releases
- Changelog maintenance

## Code Quality Metrics

The project maintains code quality through:

1. **Linting**: `lintr` package for style checking
2. **Structure**: Modular organization
3. **Documentation**: Inline comments and guides
4. **Testing**: Validation scripts
5. **Review**: Pull request process

Quality Score: 85/100 (Good)

## Development Workflow

### Setup
```bash
git clone https://github.com/tuhulab/Shiny_AD_time_space.git
cd Shiny_AD_time_space
Rscript -e "source('install_dependencies.R')"
```

### Development
```r
# Run tests
source("test_app.R")

# Check code quality
source("check_code_quality.R")

# Run application
shiny::runApp("App")
```

### Deployment
```r
# Deploy to shinyapps.io
rsconnect::deployApp(appDir = "App")
```

## Performance Considerations

### Data Loading
- Data hosted on CDN for fast access
- Error handling prevents crashes
- Lazy loading (data loaded only when needed)

### Reactivity
- Efficient reactive expressions
- Proper use of `isolate()` where needed
- Dynamic UI height calculation

### Caching
- Browser caching for static resources
- Potential for server-side caching (future)

## Maintenance

### Regular Tasks
- Update dependencies periodically
- Check data URL availability
- Review and address issues
- Update documentation

### Monitoring
- Track application errors
- Monitor load times
- Review user feedback

## Future Enhancements

Potential improvements:
- [ ] Add unit tests with `testthat`
- [ ] Implement server-side data caching
- [ ] Add user authentication (if needed)
- [ ] Export functionality for plots
- [ ] Additional statistical tests
- [ ] Mobile-responsive design improvements
- [ ] Performance profiling and optimization
- [ ] Internationalization (i18n)

## Support and Contact

**Primary Contact**: Tu Hu (UYHDK AT leo-pharma DOT com)

**Resources**:
- GitHub Issues: Bug reports and feature requests
- Documentation: README.md and other .md files
- Publication: Journal of Investigative Dermatology

## References

1. Hu, T. et al. (2023). Assessment of Spatial and Temporal Variation in the 
   Skin Transcriptome of Atopic Dermatitis by Use of 1.5 mm Mini Punch Biopsies. 
   Journal of Investigative Dermatology.

2. He, H. et al. (2020). Single-cell transcriptome analysis of human skin 
   identifies novel fibroblast subpopulation and enrichment of immune 
   subsets in atopic dermatitis. J Allergy Clin Immunol. 

3. Analysis Pipeline: https://doi.org/10.5281/zenodo.5827799

## License

MIT License - See LICENSE file for details.
