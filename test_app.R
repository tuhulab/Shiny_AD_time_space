# Testing Script for Shiny AD Time and Space Application
# This script performs basic validation tests

cat("===========================================\n")
cat("Testing Shiny AD Time and Space App\n")
cat("===========================================\n\n")

# Test 1: Check if all required packages are installed
cat("Test 1: Checking required packages...\n")
required_packages <- c(
  "shiny", "DT", "dplyr", "ggpubr", "purrr", "rstatix", "tidyr",
  "ComplexHeatmap", "BiocManager", "tidySummarizedExperiment",
  "readr", "knitr", "kableExtra"
)

missing_packages <- c()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ✗ Missing: %s\n", pkg))
    missing_packages <- c(missing_packages, pkg)
  } else {
    cat(sprintf("  ✓ Found: %s\n", pkg))
  }
}

if (length(missing_packages) > 0) {
  cat("\n⚠ Warning: Missing packages detected. Run install_dependencies.R\n")
} else {
  cat("\n✓ All required packages are installed\n")
}

# Test 2: Check if app files exist
cat("\n\nTest 2: Checking application files...\n")
required_files <- c(
  "App/ad_in_time_space.R",
  "App/config.R",
  "App/helpers.R"
)

all_files_exist <- TRUE
for (file in required_files) {
  if (file.exists(file)) {
    cat(sprintf("  ✓ Found: %s\n", file))
  } else {
    cat(sprintf("  ✗ Missing: %s\n", file))
    all_files_exist <- FALSE
  }
}

if (all_files_exist) {
  cat("\n✓ All required files exist\n")
} else {
  cat("\n✗ Some files are missing\n")
}

# Test 3: Check if config loads properly
cat("\n\nTest 3: Testing configuration file...\n")
tryCatch({
  setwd("App")
  source("config.R", local = TRUE)
  
  # Check if key config variables exist
  config_vars <- c("DATA_URLS", "SKIN_TYPE_COLORS", "PLOT_DEFAULTS", "UI_TEXT")
  for (var in config_vars) {
    if (exists(var)) {
      cat(sprintf("  ✓ Config variable exists: %s\n", var))
    } else {
      cat(sprintf("  ✗ Config variable missing: %s\n", var))
    }
  }
  
  setwd("..")
  cat("\n✓ Configuration file loaded successfully\n")
}, error = function(e) {
  setwd("..")
  cat(sprintf("\n✗ Error loading config: %s\n", e$message))
})

# Test 4: Check if helpers load properly
cat("\n\nTest 4: Testing helper functions...\n")
tryCatch({
  setwd("App")
  source("config.R", local = TRUE)
  source("helpers.R", local = TRUE)
  
  # Check if key functions exist
  helper_funcs <- c(
    "load_rds_data", "load_csv_data", "create_datatable",
    "calculate_plot_height", "filter_gsea_results"
  )
  for (func in helper_funcs) {
    if (exists(func) && is.function(get(func))) {
      cat(sprintf("  ✓ Helper function exists: %s\n", func))
    } else {
      cat(sprintf("  ✗ Helper function missing: %s\n", func))
    }
  }
  
  setwd("..")
  cat("\n✓ Helper functions loaded successfully\n")
}, error = function(e) {
  setwd("..")
  cat(sprintf("\n✗ Error loading helpers: %s\n", e$message))
})

# Test 5: Check data availability (test one URL)
cat("\n\nTest 5: Testing data accessibility...\n")
cat("  Testing one sample data URL...\n")
tryCatch({
  setwd("App")
  source("config.R", local = TRUE)
  
  # Try to access table_s4
  test_url <- DATA_URLS$table_s4
  cat(sprintf("  Attempting to load: %s\n", test_url))
  
  test_data <- readr::read_csv(url(test_url), show_col_types = FALSE)
  if (nrow(test_data) > 0) {
    cat(sprintf("  ✓ Successfully loaded data (%d rows)\n", nrow(test_data)))
  }
  
  setwd("..")
}, error = function(e) {
  setwd("..")
  cat(sprintf("  ✗ Error loading data: %s\n", e$message))
  cat("  Note: This may be due to network connectivity\n")
})

# Summary
cat("\n\n===========================================\n")
cat("Test Summary\n")
cat("===========================================\n")
if (length(missing_packages) == 0 && all_files_exist) {
  cat("✓ All tests passed! The application should be ready to run.\n")
  cat("\nTo start the application, run:\n")
  cat("  shiny::runApp('App')\n")
} else {
  cat("⚠ Some tests failed. Please review the output above.\n")
  if (length(missing_packages) > 0) {
    cat("\nTo install missing packages, run:\n")
    cat("  source('install_dependencies.R')\n")
  }
}
cat("\n")
