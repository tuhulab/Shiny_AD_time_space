# ==============================================================================
# Shiny Application: AD in Time and Space
# ==============================================================================
# Purpose: Interactive visualization of spatial and temporal variation in the
#          skin transcriptome of Atopic Dermatitis
# Author: Tu Hu (UYHDK@leo-pharma.com)
# Publication: Journal of Investigative Dermatology
# Repository: https://github.com/tuhulab/Shiny_AD_time_space
# ==============================================================================

# Load required libraries
library(shiny)
library(DT)
library(dplyr)
library(ggpubr)
library(purrr)
library(rstatix)
library(tidyr)
library(ComplexHeatmap)
library(BiocManager)
library(tidySummarizedExperiment)

# Load configuration and helper functions
source("config.R", local = TRUE)
source("helpers.R", local = TRUE)

# ==============================================================================
# User Interface (UI) Definition
# ==============================================================================

ui <- navbarPage(
  UI_TEXT$app_title,
  
  # Introduction Tab ----
  tabPanel(
    "Introduction",
    p(UI_TEXT$intro_text),
    
    h2("Download data"),
    
    h4("Raw data:"),
    tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309", 
           "Raw data is deposited on GEO (GSE193309)"),
    br(),
    
    h4("Data analysis results (models):"),
    tags$a(href = "https://doi.org/10.5281/zenodo.5827799", 
           "The reproducible data analysis pipelines is indexed on Zenodo"),
    br(),
    tags$a(href = DATA_URLS$time_variation, "Time variation data (.rds)"),
    br(),
    tags$a(href = DATA_URLS$bioreplicate, 
           "Space variation (biological replicate) data (.rds)"),
    br(),
    tags$a(href = DATA_URLS$cell_type_deconvolution, "Cell type variation data (.rds)"),
    br(),
    tags$a(href = "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/figure_s5.jpg", 
           "LINC RNA expression"),
    br(),
    tags$a(href = DATA_URLS$table_s2, "AD signature gene (.csv)"),
    br(),
    tags$a(href = DATA_URLS$gsea_res, "GSEA result (.rds)"),
    br(),
    tags$a(href = DATA_URLS$table_s4, "variance Partition result (.csv)"),
    br(),
    tags$a(href = DATA_URLS$table_s5, "Space variation (anatomic region) data (.csv)"),
    br(),
    h3("Contact"),
    p(paste("Questions regarding the usage of the web application and the data",
            "should be addressed to", UI_TEXT$contact_email))
  ),
  
  # Time Variation Tab ----
  tabPanel(
    "Time variation",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          h5(paste("Boxplot showing time-specific DEGs. It appears that the",
                   "differential expression is driven by a few high-count samples.")),
          radioButtons(
            "time", 
            h4("Choose how to visualize time:"),
            choices = list(
              "quarter (absolute time)" = "quarter",
              "visit (relative time)" = "visit"
            ),
            selected = "visit"
          ),
          radioButtons(
            "skin_type", 
            h4("Choose skin type (LS/NL/HC) to visualize:"),
            choices = list(
              "AD lesional skin (LS)" = "LS",
              "AD non-lesional skin (NL)" = "NL",
              "Healthy control skin (HC)" = "HC"
            ),
            selected = "LS"
          )
        ),
        mainPanel(uiOutput("timevariation.ui"))
      )
    )
  ),
  
  # Space Variation (Biological Replicate) Tab ----
  tabPanel(
    "Space variation (biological replicate)",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          h5(paste("Boxplot showing gene expression for selected AD biomarkers from",
                   "biological replicates. The biological replicate pairs, which were",
                   "taken from the same anatomical region from the same subject at the",
                   "same time point, are connected by lines. Although the group means",
                   "are similar, the gene expression from two biological replicates",
                   "shows considerable variation.")),
          checkboxGroupInput(
            "gene_to_plot_duplicate",
            "Choose the genes to plot:",
            choices = load_rds_data(DATA_URLS$bioreplicate, "bioreplicate data")$feature %>% unique(),
            selected = load_rds_data(DATA_URLS$bioreplicate, "bioreplicate data")$feature %>% head(5)
          )
        ),
        mainPanel(uiOutput("bioreplicate.ui"))
      )
    )
  ),
  
  # Cell Type Variation Tab ----
  tabPanel(
    "Cell type variation",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          h5(paste("The inferred cell composition for each sample. Each bar represents",
                   "the proportion (%) of cell composition, as estimated by MuSiC",
                   "(Multi-Subject Single Cell deconvolution) using GSE147424 (He et al. 2020)",
                   "single-cell data as the reference.")),
          checkboxGroupInput(
            "subject_to_plot",
            "Choose the subject IDs to plot:",
            choices = load_rds_data(DATA_URLS$cell_type_deconvolution, "cell type data") %>% 
              pull(subject) %>% 
              unique(),
            selected = c("AD_02", "AD_06")
          )
        ),
        mainPanel(uiOutput("celltypevariation.ui"))
      )
    )
  ),
  
  # LINC RNA Tab ----
  tabPanel(
    "LINC RNA",
    basicPage(
      h2("Heatmap for LINC RNA"),
      br(),
      h5("Heatmap showing 177 differentially expressed lincRNAs between skin conditions."),
      imageOutput("heatmap_linc")
    )
  ),
  
  # AD Signature Gene Tab ----
  tabPanel(
    "",
    basicPage(
      h2("AD signature gene"),
      p("AD signature gene"),
      br(),
      DT::dataTableOutput("ad_signature_table")
    )
  ),
  
  # GSEA Result Tab ----
  tabPanel(
    "GSEA result",
    fluidPage(
      title = "Gene set enrichment analysis results",
      sidebarPanel(
        checkboxGroupInput(
          "study_to_show", 
          h4("Choose studies to show: (TO BE IMPLEMENTED, all by default)"),
          c("Hu et al (GSE193309)", "Tsoi et al (GSE121212)", "Acute_chronic", 
            "GSE107361_adult", "GSE107361_children", "GSE12721"),
          selected = c("Hu et al (GSE193309)", "Tsoi et al (GSE121212)", "Acute_chronic", 
                      "GSE107361_adult", "GSE107361_children", "GSE12721")
        ),
        radioButtons(
          "gsea_res_parameter", 
          h4("Choose parameter to show:"),
          choices = list(
            "Normalized enrichment score" = "NES",
            "p-value adjusted" = "padj"
          ),
          selected = "NES"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = "dataset",
          tabPanel("Reactome", DT::dataTableOutput("gsea_table_reactome")),
          tabPanel("Transcription Factor Target", DT::dataTableOutput("gsea_table_tft")),
          tabPanel("GO: Biological Process", DT::dataTableOutput("gsea_table_bp")),
          tabPanel("GO: Molecular Function", DT::dataTableOutput("gsea_table_mf"))
        )
      )
    )
  ),
  
  # Variance Partition Tab ----
  tabPanel(
    "Variance Partition",
    basicPage(
      h2("VP results"),
      p("variance partition results for each gene"),
      br(),
      DT::dataTableOutput("vp_table")
    )
  ),
  
  # Space Variation (Anatomic Region) Tab ----
  tabPanel(
    "Space variation (anatomic region)",
    basicPage(
      h2("Anatomic region variation"),
      br(),
      DT::dataTableOutput("space_anatomic_region")
    )
  ),
  
  # Theoretical RNA Concentration Tab ----
  tabPanel(
    "Theoretical calculation of RNA conc.",
    basicPage(
      h4("Theoretical calculation of RNA concentration from different biopsy diameters"),
      tableOutput("RNA_conc_kable")
    )
  ),
  
  # Cell Composition Tab ----
  tabPanel(
    "Cell composition (GSE121212 and GENAD)",
    basicPage(
      p("Cell composition of different skin types (HC/NL/LS) for",
        "the", strong("GSE121212"), "and", strong("GENAD"),
        "study inferred by", em("in silico"), "cell sorting analysis"),
      DT::dataTableOutput("table_cell_composition_gse121212_genad")
    )
  ),
  
  # Tissue Injury Gene Expression Tab ----
  tabPanel(
    "Tissue injury gene expression",
    basicPage(
      p("Longitudinal (visit) variation of gene expression for MGST1, MUC1, PTGS2, and SAA2"),
      imageOutput("tissue_injury_g")
    )
  ),
  
  # Subcutis Gene Expression Tab ----
  tabPanel(
    "Subcutis gene expression",
    basicPage(
      p("Longitudinal (visit) variation of gene expression for CIDEC, FABP4, and PLIN4"),
      imageOutput("subcuits_g")
    )
  )
)


# ==============================================================================
# Server Logic
# ==============================================================================

server <- function(input, output) {
  
  # Variance Partition Table ----
  output$vp_table <- DT::renderDataTable({
    load_csv_data(DATA_URLS$table_s4, "variance partition data") %>%
      dplyr::select(-1) %>%
      dplyr::arrange(-`tissue type (LS/NL/HC)`) %>%
      create_datatable(page_length = PLOT_DEFAULTS$page_length) %>%
      DT::formatPercentage(2:9)
  })
  
  # AD Signature Gene Table ----
  output$ad_signature_table <- DT::renderDataTable({
    load_csv_data(DATA_URLS$table_s2, "AD signature data") %>%
      create_datatable(page_length = PLOT_DEFAULTS$page_length) %>%
      DT::formatRound(2:4, digits = 3) %>%
      DT::formatSignif(5:6)
  })
  
  # GSEA Tables ----
  # Helper function to render GSEA table for a specific gene set
  render_gsea_table <- function(gene_set_type) {
    DT::renderDataTable({
      gsea_data <- load_rds_data(DATA_URLS$gsea_res, "GSEA results")
      if (is.null(gsea_data)) return(NULL)
      
      filter_gsea_results(gsea_data, gene_set_type, input$gsea_res_parameter) %>%
        create_gsea_table()
    })
  }
  
  output$gsea_table_reactome <- render_gsea_table("REACTOME")
  output$gsea_table_tft <- render_gsea_table("TFT")
  output$gsea_table_bp <- render_gsea_table("BP")
  output$gsea_table_mf <- render_gsea_table("MF")
  
  # Space Variation (Anatomic Region) Table ----
  output$space_anatomic_region <- DT::renderDataTable({
    load_csv_data(DATA_URLS$table_s5, "anatomic region data") %>%
      create_datatable(page_length = PLOT_DEFAULTS$page_length) %>%
      DT::formatRound(4:6, digits = 3) %>%
      DT::formatSignif(7:8)
  })
  
  # Biological Replicate Plot ----
  output$bioreplicate_g <- renderPlot({
    data <- load_rds_data(DATA_URLS$bioreplicate, "bioreplicate data")
    if (is.null(data)) return(NULL)
    
    data_subset <- data %>%
      dplyr::select(biological_rep_id, subject, visit, skin_type, 
                   feature, counts_scaled, replicate_ID)
    
    create_bioreplicate_plot(
      data_subset, 
      input$gene_to_plot_duplicate, 
      SKIN_TYPE_COLORS
    )
  })
  
  output$bioreplicate.ui <- renderUI({
    plotOutput(
      "bioreplicate_g", 
      width = "100%",
      height = calculate_plot_height(length(input$gene_to_plot_duplicate))
    )
  })
  
  # Cell Type Variation Plot ----
  output$cell_type_variation_g <- renderPlot({
    data <- load_rds_data(DATA_URLS$cell_type_deconvolution, "cell type data")
    if (is.null(data)) return(NULL)
    
    data %>%
      dplyr::filter(subject %in% input$subject_to_plot) %>%
      ggpubr::ggbarplot(
        x = "plot_name", 
        y = "prop", 
        fill = "cell_type",
        xlab = "Sample", 
        ylab = "Proportion(%)", 
        palette = "npg"
      ) %>%
      ggpubr::facet(facet.by = "subject", scales = "free", ncol = 2) %>%
      ggpubr::ggpar(xtickslab.rt = 90, font.tickslab = 6)
  })
  
  output$celltypevariation.ui <- renderUI({
    n_subjects <- length(input$subject_to_plot)
    height <- 300 * round((n_subjects + 0.01) / 2)
    plotOutput("cell_type_variation_g", width = "100%", height = paste0(height, "px"))
  })
  
  # LINC RNA Heatmap ----
  output$heatmap_linc <- renderPlot({
    heatmap <- load_rds_data(DATA_URLS$heatmap_linc, "LINC RNA heatmap")
    if (is.null(heatmap)) return(NULL)
    
    ComplexHeatmap::draw(
      heatmap, 
      heatmap_legend_side = "right",
      annotation_legend_side = "top", 
      merge_legend = FALSE
    )
  }, height = 1500)
  
  # Time Variation Plot ----
  output$time_variation_g <- renderPlot({
    data <- load_rds_data(DATA_URLS$time_variation, "time variation data")
    if (is.null(data)) return(NULL)
    
    data %>%
      dplyr::filter(skin_type == input$skin_type, time_type == input$time) %>%
      dplyr::pull(d) %>% 
      purrr::reduce(~.x) %>%
      ggpubr::ggboxplot(
        x = ifelse(input$time == "quarter", "visit_quarter", "visit"),
        y = "counts_scaled",
        add = "jitter", 
        add.params = list(alpha = 0.5),
        facet.by = "feature",
        ylab = FALSE,
        yscale = PLOT_DEFAULTS$log_scale,
        scales = "free",
        ncol = 3
      )
  })
  
  output$timevariation.ui <- renderUI({
    # Calculate height based on conditions
    height_multiplier <- dplyr::case_when(
      input$time == "quarter" & input$skin_type == "LS" ~ 33,
      input$time == "quarter" & input$skin_type == "NL" ~ 3,
      input$time == "quarter" & input$skin_type == "HC" ~ 1.3,
      input$time == "visit" & input$skin_type == "LS" ~ 12,
      input$time == "visit" & input$skin_type == "NL" ~ 24,
      input$time == "visit" & input$skin_type == "HC" ~ 8
    )
    
    width <- ifelse(input$time == "quarter" & input$skin_type == "HC", "50%", "100%")
    height <- paste0(height_multiplier * PLOT_DEFAULTS$height_per_feature, "px")
    
    plotOutput("time_variation_g", width = width, height = height)
  })
  
  # RNA Concentration Table ----
  output$RNA_conc_kable <- function() {
    load_csv_data(DATA_URLS$rna_conc, "RNA concentration data") %>%
      knitr::kable("html") %>%
      kableExtra::kable_styling("striped", full_width = FALSE) %>%
      kableExtra::add_footnote(
        "Tsoi et al. 2019 has reported that from a 5 mm punch biopsy, 50 ng/uL RNA can be obtained for high quality sequencing."
      )
  }
  
  # Cell Composition Table ----
  output$table_cell_composition_gse121212_genad <- DT::renderDataTable({
    load_csv_data(DATA_URLS$cell_composition, "cell composition data") %>%
      DT::datatable(rownames = FALSE)
  })
  
  # Tissue Injury Gene Expression Plot ----
  output$tissue_injury_g <- renderPlot({
    data <- load_rds_data(DATA_URLS$tissue_injury, "tissue injury data")
    if (is.null(data)) return(NULL)
    
    data %>%
      ggpubr::ggline(
        x = "visit",
        y = "counts_scaled",
        facet.by = "feature",
        add = c("mean", "dotplot"),
        color = "skin_type",
        palette = SKIN_TYPE_COLORS,
        xlab = "Visit",
        ylab = "Counts"
      ) %>%
      apply_plot_theme(PLOT_DEFAULTS$log_scale)
  })
  
  # Subcutis Gene Expression Plot ----
  output$subcuits_g <- renderPlot({
    data <- load_rds_data(DATA_URLS$subcutis, "subcutis data")
    if (is.null(data)) return(NULL)
    
    data %>%
      ggpubr::ggline(
        x = "visit",
        y = "counts_scaled",
        facet.by = "feature",
        add = c("mean", "dotplot"),
        color = "skin_type",
        palette = SKIN_TYPE_COLORS,
        xlab = "Visit",
        ylab = "Counts"
      ) %>%
      apply_plot_theme(PLOT_DEFAULTS$log_scale)
  })
}


# ==============================================================================
# Run the Shiny Application
# ==============================================================================
shinyApp(ui = ui, server = server)
