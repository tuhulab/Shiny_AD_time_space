# Code Quality Assessment Script
# This script performs code quality checks and generates metrics

# Install required packages if not already installed
if (!requireNamespace("lintr", quietly = TRUE)) {
  install.packages("lintr")
}

library(lintr)

# Define files to lint
files_to_lint <- c(
  "App/ad_in_time_space.R",
  "App/config.R",
  "App/helpers.R",
  "prototype_data.R"
)

# Perform linting
cat("=== Code Quality Assessment ===\n\n")
cat("Running lintr on all R files...\n\n")

lint_results <- list()
for (file in files_to_lint) {
  if (file.exists(file)) {
    cat(sprintf("Linting %s...\n", file))
    lints <- lintr::lint(file)
    lint_results[[file]] <- lints
    
    if (length(lints) > 0) {
      print(lints)
    } else {
      cat("  ✓ No linting issues found\n")
    }
    cat("\n")
  } else {
    cat(sprintf("  ⚠ File not found: %s\n\n", file))
  }
}

# Calculate metrics
total_lints <- sum(sapply(lint_results, length))
files_checked <- length(lint_results)
clean_files <- sum(sapply(lint_results, length) == 0)

# Generate summary
cat("=== Summary ===\n")
cat(sprintf("Files checked: %d\n", files_checked))
cat(sprintf("Clean files (no issues): %d\n", clean_files))
cat(sprintf("Total linting issues: %d\n", total_lints))

# Calculate quality score (0-100)
# Score = 100 - (issues per file * 5, capped at 100)
if (files_checked > 0) {
  issues_per_file <- total_lints / files_checked
  quality_score <- max(0, min(100, 100 - (issues_per_file * 5)))
} else {
  quality_score <- 0
}

cat(sprintf("Code Quality Score: %.1f/100\n", quality_score))

# Determine badge color and label
if (quality_score >= 90) {
  badge_color <- "brightgreen"
  badge_label <- "excellent"
} else if (quality_score >= 75) {
  badge_color <- "green"
  badge_label <- "good"
} else if (quality_score >= 60) {
  badge_color <- "yellowgreen"
  badge_label <- "fair"
} else if (quality_score >= 50) {
  badge_color <- "yellow"
  badge_label <- "needs improvement"
} else {
  badge_color <- "red"
  badge_label <- "poor"
}

# Generate badge URL
badge_url <- sprintf(
  "https://img.shields.io/badge/code%%20quality-%.0f%%2F100%%20(%s)-%s",
  quality_score,
  gsub(" ", "%%20", badge_label),
  badge_color
)

cat("\n=== Badge Information ===\n")
cat(sprintf("Badge URL: %s\n", badge_url))
cat(sprintf("Markdown: ![Code Quality](%s)\n", badge_url))

# Save results to file
writeLines(
  c(
    sprintf("Code Quality Assessment - %s", Sys.Date()),
    "",
    sprintf("Files checked: %d", files_checked),
    sprintf("Clean files: %d", clean_files),
    sprintf("Total issues: %d", total_lints),
    sprintf("Quality Score: %.1f/100 (%s)", quality_score, badge_label),
    "",
    sprintf("Badge: ![Code Quality](%s)", badge_url)
  ),
  "CODE_QUALITY.txt"
)

cat("\nResults saved to CODE_QUALITY.txt\n")
