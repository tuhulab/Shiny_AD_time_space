# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-01-15

### Added
- Comprehensive DESCRIPTION file with all package dependencies
- REPRODUCIBILITY.md guide for ensuring reproducible analysis
- LICENSE file (MIT License)
- CODE_STYLE.md with coding conventions and best practices
- Modular code structure:
  - `App/config.R` - Centralized configuration
  - `App/helpers.R` - Reusable helper functions
- `install_dependencies.R` - Automated dependency installation script
- `test_app.R` - Basic validation and testing script
- `check_code_quality.R` - Code quality assessment tool
- Code quality badge in README
- Detailed inline documentation and comments

### Changed
- Refactored `App/ad_in_time_space.R` for better readability and maintainability
- Improved `README.md` with comprehensive documentation
- Enhanced `prototype_data.R` with better structure and comments
- Standardized naming conventions throughout the codebase
- Improved error handling for data loading operations

### Improved
- Code modularity - separated concerns into config, helpers, and main app
- Documentation - added comprehensive inline comments and function documentation
- Maintainability - reduced code duplication through helper functions
- Performance - added error handling to prevent crashes on data loading failures
- Reproducibility - documented all dependencies and setup procedures

### Fixed
- Typo: "facilitaes" → "facilitates" in UI text
- Inconsistent use of `formatPercentage` and `formatRound` in tables
- Repetitive GSEA table rendering code

## [0.1.0] - Initial Release

### Added
- Initial Shiny application for AD time and space visualization
- Data visualization for time variation, space variation, and cell types
- GSEA results browser
- Variance partition results
- Interactive plots and tables
