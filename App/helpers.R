# Helper functions for Shiny AD Time and Space Application
# This file contains reusable utility functions to reduce code duplication

#' Load RDS data from URL with error handling
#'
#' @param url Character string of the data URL
#' @param description Optional description for error messages
#' @return Loaded R object or NULL if error occurs
load_rds_data <- function(url, description = "data") {
  tryCatch(
    readRDS(url(url)),
    error = function(e) {
      message(sprintf("Error loading %s from %s: %s", description, url, e$message))
      NULL
    }
  )
}

#' Load CSV data from URL with error handling
#'
#' @param url Character string of the data URL
#' @param description Optional description for error messages
#' @return Tibble or NULL if error occurs
load_csv_data <- function(url, description = "data") {
  tryCatch(
    readr::read_csv(url, show_col_types = FALSE),
    error = function(e) {
      message(sprintf("Error loading %s from %s: %s", description, url, e$message))
      NULL
    }
  )
}

#' Create a formatted DataTable with standard options
#'
#' @param data Data frame or tibble to display
#' @param page_length Number of rows per page (default: 30)
#' @param row_names Whether to show row names (default: FALSE)
#' @return DT::datatable object
create_datatable <- function(data, page_length = 30, row_names = FALSE) {
  DT::datatable(
    data,
    rownames = row_names,
    options = list(pageLength = page_length)
  )
}

#' Calculate dynamic plot height based on number of features
#'
#' @param n_features Number of features to plot
#' @param height_per_feature Height in pixels per feature (default: 200)
#' @return Character string with height in pixels (e.g., "600px")
calculate_plot_height <- function(n_features, height_per_feature = 200) {
  paste0(n_features * height_per_feature, "px")
}

#' Filter GSEA results by gene set and parameter
#'
#' @param gsea_data GSEA results data
#' @param gene_set Gene set type (e.g., "REACTOME", "TFT", "BP", "MF")
#' @param parameter Parameter type (e.g., "NES", "padj")
#' @return Filtered data frame
filter_gsea_results <- function(gsea_data, gene_set, parameter) {
  gsea_data %>%
    dplyr::filter(gs == gene_set) %>%
    dplyr::pull(data) %>%
    purrr::reduce(~ .x[[1]]) %>%
    dplyr::filter(parameter == !!parameter) %>%
    dplyr::select(-parameter)
}

#' Create GSEA DataTable with standard formatting
#'
#' @param data Filtered GSEA data
#' @return DT::datatable object with HTML formatting
create_gsea_table <- function(data) {
  DT::datatable(
    data,
    rownames = FALSE,
    escape = FALSE,  # Allow HTML in pathway names
    options = list(
      pageLength = 30,
      autoWidth = TRUE,
      columnDefs = list(list(width = '30%', targets = 0))
    )
  )
}

#' Apply standard boxplot theme and styling
#'
#' @param plot ggplot object
#' @param y_scale Y-axis scale transformation (default: "log10")
#' @return Modified ggplot object
apply_plot_theme <- function(plot, y_scale = "log10") {
  ggpubr::ggpar(plot, yscale = y_scale)
}

#' Create paired biological replicate plot
#'
#' @param data Data frame with biological replicate data
#' @param features Vector of feature names to plot
#' @param skin_colors Named vector of colors for skin types
#' @return ggplot object
create_bioreplicate_plot <- function(data, features, skin_colors) {
  # Filter data for selected features
  plot_data <- data %>%
    dplyr::filter(feature %in% features) %>%
    dplyr::arrange(biological_rep_id) %>%
    dplyr::mutate(counts_scaled = counts_scaled + 1)
  
  # Calculate statistical tests
  stat_test <- plot_data %>%
    dplyr::group_by(feature, skin_type) %>%
    rstatix::t_test(counts_scaled ~ replicate_ID, paired = TRUE) %>%
    rstatix::add_significance() %>%
    rstatix::add_xy_position(x = "replicate_ID", y.trans = log2, step.increase = 0.08) %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(y.position = max(y.position)) %>%
    dplyr::ungroup()
  
  # Create plot
  plot_data %>%
    tidyr::pivot_wider(names_from = replicate_ID, values_from = counts_scaled) %>%
    ggpubr::ggpaired(
      cond1 = "01", cond2 = "02",
      y = "counts_scaled",
      xlab = "Biological replicate",
      ylab = "",
      line.color = "grey",
      line.size = 0.4,
      fill = "skin_type",
      palette = skin_colors
    ) %>%
    ggpubr::facet(facet.by = c("feature", "skin_type"), scales = "free_y") +
    ggpubr::stat_pvalue_manual(stat_test, label = "p") +
    ggplot2::scale_y_continuous(trans = "log2", expand = ggplot2::expansion(mult = c(0.05, 0.1)))
}
