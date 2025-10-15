# ==============================================================================
# Prototype Data Analysis Script
# ==============================================================================
# Purpose: Exploratory data analysis and visualization prototypes for the
#          Shiny AD Time and Space application
# Note: This script uses local data files from the data/ directory
# ==============================================================================

# Load required libraries ----
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)

# Variance Partition Table ----
# Display variance partition results with formatted percentages
readr::read_csv("data/table_s4.csv") %>%
  dplyr::select(-1) %>%
  dplyr::arrange(-`tissue type (LS/NL/HC)`) %>%
  DT::datatable(rownames = FALSE) %>%
  DT::formatPercentage(2:9)

# GSEA Results Processing ----
# Process GSEA results from Excel file with multiple sheets
gsea_res <-
  tibble(
    sheet_name = openxlsx::getSheetNames("data/table_s3_gsea_res.xlsx")
  ) %>%
  mutate(
    gs = sheet_name %>% str_extract("REACTOME|TFT|BP|MF"),
    parameter = sheet_name %>% str_extract("NES|padj"),
    content = map(
      sheet_name, 
      ~ openxlsx::read.xlsx("data/table_s3_gsea_res.xlsx", .x) %>%
        as_tibble() %>%
        mutate(
          across(-1:-3, ~ round(., 4)),
          pathway = str_remove(pathway, "REACTOME_"),
          pathway = paste0("<a href='", url, "'>", pathway, "</a>")
        ) %>%
        select(-url, -n)
    ),
    content = map2(content, parameter, ~ mutate(.x, parameter = .y))
  ) %>%
  select(-sheet_name, -parameter) %>%
  group_by(gs) %>% 
  nest() %>%
  mutate(data = map(data, ~.x[[1]] %>% reduce(bind_rows)))

# Save processed GSEA results
saveRDS(gsea_res, "data/gsea_res.rds")

# Display GSEA results for Reactome pathways
gsea_res_display <-
  gsea_res %>% 
  filter(sheet_name == "CP_REACTOME_NES") %>%
  pull(content) %>% 
  reduce(~ .x[[1]]) %>%
  DT::datatable(
    rownames = FALSE, 
    escape = FALSE,
    options = list(
      pageLength = 15,
      autoWidth = TRUE,
      columnDefs = list(list(width = '30%', targets = c(0, 1)))
    )
  )
gsea_res_display

# Biological Replicate Visualization ----
# Load biological replicate data
bioreplicate_t <- readr::read_rds("data/bioreplicate_t.rds")

# Create visualization of biological replicate variation
replicate_g <-
  bioreplicate_t %>%
  select(
    biological_rep_id, subject, visit, skin_type, 
    feature, counts_scaled, replicate_ID
  ) %>%
  filter(feature %in% (bioreplicate_t$feature %>% unique() %>% head(5))) %>%
  ggplot(aes(
    x = replicate_ID,
    y = counts_scaled,
    group = biological_rep_id
  )) +
  geom_point(aes(color = skin_type)) +
  facet_grid(feature ~ skin_type, scales = "free") +
  geom_line(aes(color = skin_type), alpha = 0.5) +
  geom_boxplot(aes(group = replicate_ID, fill = skin_type), alpha = 0.5) +
  scale_color_manual(values = c(
    "LS" = "#eb2d0c", 
    "NL" = "#eb8b9b", 
    "HC" = "#91cf60"
  )) +
  scale_fill_manual(values = c(
    "LS" = "#eb2d0c", 
    "NL" = "#eb8b9b", 
    "HC" = "#91cf60"
  )) +
  scale_y_continuous(trans = "log2") +
  xlab("replicate") +
  theme(axis.title.y = element_blank())

# Create interactive plot
plotly::ggplotly(replicate_g)

# Time Variation Visualization ----
# Create boxplot for time variation of gene expression
time_variation_g <-
  readr::read_rds("data/time_variation_t.rds") %>%
  filter(skin_type == "LS", time_type == "quarter") %>%
  pull(d) %>% 
  reduce(~.x) %>%
  ggpubr::ggboxplot(
    x = ifelse("quarter" == "quarter", "visit_quarter", "visit"),
    y = "counts_scaled",
    add = "jitter",
    facet.by = "feature",
    ylab = FALSE,
    yscale = "log10",
    scales = "free"
  )






