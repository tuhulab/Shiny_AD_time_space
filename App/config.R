# Configuration file for Shiny AD Time and Space Application
# This file centralizes all configuration parameters for easier maintenance

# Data URLs (CDN)
DATA_BASE_URL <- "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/"

DATA_URLS <- list(
  time_variation = paste0(DATA_BASE_URL, "time_variation_t.rds"),
  bioreplicate = paste0(DATA_BASE_URL, "bioreplicate_t.rds"),
  cell_type_deconvolution = paste0(DATA_BASE_URL, "cell_type_deconvolution.rds"),
  heatmap_linc = paste0(DATA_BASE_URL, "heatmap_linc.rds"),
  gsea_res = paste0(DATA_BASE_URL, "gsea_res.rds"),
  table_s2 = paste0(DATA_BASE_URL, "table_s2.csv"),
  table_s4 = paste0(DATA_BASE_URL, "table_s4.csv"),
  table_s5 = paste0(DATA_BASE_URL, "table_s5.csv"),
  rna_conc = paste0(DATA_BASE_URL, "table_theo_calc_RNA_conc.csv"),
  cell_composition = paste0(DATA_BASE_URL, "table_cell_composition_gse121212_genad.csv"),
  tissue_injury = paste0(DATA_BASE_URL, "tissue_injury_se_t.rds"),
  subcutis = paste0(DATA_BASE_URL, "subcutis_se_t.rds")
)

# Color palette for skin types
SKIN_TYPE_COLORS <- c(
  "LS" = "#eb2d0c",  # Lesional skin - red
  "NL" = "#eb8b9b",  # Non-lesional skin - pink
  "HC" = "#91cf60"   # Healthy control - green
)

# Default plot parameters
PLOT_DEFAULTS <- list(
  height_per_feature = 200,  # pixels per feature in plots
  page_length = 30,           # default rows per page in tables
  log_scale = "log10"         # default log scale for plots
)

# UI text constants
UI_TEXT <- list(
  app_title = "AD in Time and Space",
  intro_text = paste(
    "This web interface facilitates the readers to interact with and",
    "download the data published on Hu. et al,",
    "Assessment of Spatial and Temporal Variation in the Skin Transcriptome",
    "of Atopic Dermatitis by Use of 1.5 mm Mini Punch Biopsies,",
    "Journal of Investigative Dermatology (in press)"
  ),
  contact_email = "Tu Hu (UYHDK AT leo-pharma DOT com)"
)
