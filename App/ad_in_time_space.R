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
# options(repos = BiocManager::repositories())
# Sys.setenv(CMAKE_BIN="/usr/bin/cmake3")

# Define UI ----
ui <- navbarPage("AD in Time and Space",
                 tabPanel("Introduction",
                          p("This web interface facilitaes the readers to interact with and download the data published on Hu. et al,
                          Assessment of Spatial and Temporal Variation in the Skin Transcriptome of Atopic Dermatitis by Use of 1.5 mm Mini Punch Biopsies,
                            Journal of Investigative Dermatology (in press)"),

                          h2("Download data"),

                          h4("Raw data:"),
                          tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309", "Raw data is deposited on GEO (GSE193309)"),
                          br(),

                          h4("Data analysis results (models):"),
                          tags$a(href = "https://doi.org/10.5281/zenodo.5827799", "The reproducible data analysis pipelines is indexed on Zenodo"),
                          br(),
                          tags$a(href = "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/time_variation_t.rds", "Time variation data (.rds)"),
                          br(),
                          tags$a(href = "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/bioreplicate_t.rds", "Space variation (biological replicate) data (.rds)"),
                          br(),
                          tags$a(href = "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/cell_type_deconvolution.rds", "Cell type variation data (.rds)"),
                          br(),
                          tags$a(href = "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/figure_s5.jpg", "LINC RNA expression"),
                          br(),
                          tags$a(href = "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/table_s2.csv", "AD signature gene (.csv)"),
                          br(),
                          tags$a(href = "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/gsea_res.rds", "GSEA result (.rds)"),
                          br(),
                          tags$a(href = "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/table_s4.csv", "variance Parition result (.csv)"),
                          br(),
                          tags$a(href = "https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/table_s5.csv", "Space variation (anatomic region) data (.csv)"),
                          br(),
                          h3("Contact"),
                          p("Questions regarding the usage of the web application and the data should be addressed to Tu Hu (UYHDK AT leo-pharma DOT com)")
                          ),
                 tabPanel("Time variation", # Figure S2
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel(h5("Boxplot showing time-specific DEGs. It appears that the differential expression is driven by a few high-count samples."),
                                           radioButtons(
                                             "time", h4("Choose how to visualize time:"),
                                             choices = list("quarter (absolute time)" = "quarter",
                                                            "visit (relative time)" = "visit"),
                                             selected = "visit"),
                                           radioButtons(
                                             "skin_type", h4("Choose skin type (LS/NL/HC) to visualize:"),
                                             choices = list("AD lesional skin (LS)" = "LS",
                                                            "AD non-lesional skin (NL)" = "NL",
                                                            "Healthy control skin (HC)" = "HC"),
                                             selected = "LS"
                                           )
                                           ),
                              mainPanel(
                                # plotOutput("time_variation_g", width = "100%", height = "800px")
                                uiOutput("timevariation.ui")
                              )
                            )
                          )),
                 tabPanel("Space variation (biological replicate)",
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel(
                                  h5("Boxplot showing gene expression for selected AD biomarkers from biological replicates. The biological replicate pairs, which were taken from the same anatomical region from the same subject at the same time point, are connected by lines. Although the group means are similar, the gene expression from two biological replicates shows considerable variation."),
                                  checkboxGroupInput("gene_to_plot_duplicate",
                                                     "Choose the genes to plot:",
                                                     choices = readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/bioreplicate_t.rds"))$feature %>% unique,
                                                     selected = readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/bioreplicate_t.rds"))$feature %>% head(5))
                                           ),
                              mainPanel(
                                uiOutput("bioreplicate.ui")
                              )
                            )
                          )), # Figure S3
                 tabPanel("Cell type variation",
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel(h5("The inferred cell composition for each sample. Each bar represents the proportion (%) of cell composition, as estimated by MuSiC (Multi-Subject Single Cell deconvolution) using GSE147424 (He et al. 2020) single-cell data as the reference."),
                                           checkboxGroupInput(
                                             "subject_to_plot",
                                             "Choose the subject IDs to plot:",
                                             choices = readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/cell_type_deconvolution.rds")) %>% pull(subject) %>% unique(),
                                             selected = c("AD_02", "AD_06")
                                           )),
                              mainPanel(
                                # plotOutput("cell_type_variation_g", height = "800px")
                                uiOutput("celltypevariation.ui")

                              )
                            )
                          )), # Figure S4
                 tabPanel("LINC RNA",
                          basicPage(
                            h2("Heatmap for LINC RNA"),
                            br(),
                            h5("Heatmap showing 177 differentially expressed lincRNAs between skin conditions."),
                            imageOutput("heatmap_linc")
                          )), # Figure S5
                 tabPanel("",
                          basicPage(
                            h2("AD signature gene"),
                            p("AD signature gene"),
                            br(),
                            DT::dataTableOutput("ad_signature_table")
                          )), # Table S2
                 tabPanel("GSEA result",
                          fluidPage(
                            title = "Gene set enrichment analysis results",
                            sidebarPanel(
                              checkboxGroupInput("study_to_show", h4("Choose studies to show: (TO BE IMPLEMENTED, all by default)"),
                                                 c("Hu et al (GSE193309)", "Tsoi et al (GSE121212)", "Acute_chronic", "GSE107361_adult", "GSE107361_children", "GSE12721"),
                                                 selected =c("Hu et al (GSE193309)", "Tsoi et al (GSE121212)", "Acute_chronic", "GSE107361_adult", "GSE107361_children", "GSE12721")),
                              radioButtons(
                                "gsea_res_parameter", h4("Choose parameter to show:"),
                                choices = list("Normalized enrichment score" = "NES",
                                               "p-value adjusted" = "padj"),
                                selected = "NES")),
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
                          ), # Table S3
                 tabPanel("Variance Partition",
                          basicPage(
                            h2("VP results"),
                            p("variance parition results for each gene"),
                            br(),
                            DT::dataTableOutput("vp_table")
                          )
                          ), # Table S4
                 tabPanel("Space variation (anatomic region)",
                          basicPage(
                            h2("Anatomic region variation"),
                            br(),
                            DT::dataTableOutput("space_anatomic_region")
                          )), # Table S5)
                 tabPanel("Theoritical calculation of RNA conc.",
                          basicPage(
                            h4("Theoretical calculation of RNA concentration from different biopsy diameters"),
                            tableOutput("RNA_conc_kable")
                          )),
                 tabPanel("Cell composition (GSE121212 and GENAD)",
                          basicPage(
                            p("Cell composition of different skin types (HC/NL/LS) for",
                              "the",
                              strong("GSE121212"),
                              "and",
                              strong("GENAD"),
                              "study inferred by",
                              em("in silico"),
                              "cell sorting analysis"),
                          DT::dataTableOutput("table_cell_composition_gse121212_genad"))
                          ),
                 tabPanel("Tissue injury gene expression",
                          basicPage(
                            p("Longitudinal (visit) variation of gene expression for MGST1, MUC1, PTGS2, and SAA2"),
                            imageOutput("tissue_injury_g")
                          )),
                 tabPanel("Subcutis gene expression",
                          basicPage(
                            p("Longitudinal (visit) variation of gene expression for CIDEC, FABP4, and PLIN4"),
                            imageOutput("subcuits_g")
                          )),
)

# Define server logic ----
server <- function(input, output) {
  output$vp_table <- DT::renderDataTable({
    readr::read_csv(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/table_s4.csv")) %>%
      dplyr::select(-1) %>%
      arrange(-`tissue type (LS/NL/HC)`) %>%
      DT::datatable(rownames = FALSE, options = list(pageLength = 30)) %>%
      formatPercentage(2:9)
  })

  output$ad_signature_table <-
    DT::renderDataTable({
      readr::read_csv(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/table_s2.csv")) %>%
        DT::datatable(rownames = FALSE, options = list(pageLength = 30)) %>%
        formatRound(2:4, digits = 3) %>%
        formatSignif(5:6)
    })

  output$gsea_table_reactome <- DT::renderDataTable({
    readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/gsea_res.rds")) %>%
      filter(gs == "REACTOME") %>%
      pull(data) %>% purrr::reduce(~ .x[[1]]) %>%
      filter(parameter == input$gsea_res_parameter) %>%
      select(-parameter) %>%
      DT::datatable(rownames = FALSE, escape = F,
                    options = list(
                      pageLength = 30,
                      autoWidth = T,
                      columnDefs = list(list(width = '30%', targets = 0))
                      ))
    })

  output$gsea_table_tft <- DT::renderDataTable({
    readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/gsea_res.rds")) %>%
      filter(gs == "TFT") %>%
      pull(data) %>% purrr::reduce(~ .x[[1]]) %>%
      filter(parameter == input$gsea_res_parameter) %>%
      select(-parameter) %>%
      DT::datatable(rownames = FALSE, escape = F,
                    options = list(
                      pageLength = 30,
                      autoWidth = T,
                      columnDefs = list(list(width = '30%', targets = 0))
                    ))
  })

  output$gsea_table_bp <- DT::renderDataTable({
    readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/gsea_res.rds")) %>%
      filter(gs == "BP") %>%
      pull(data) %>% purrr::reduce(~ .x[[1]]) %>%
      filter(parameter == input$gsea_res_parameter) %>%
      select(-parameter) %>%
      DT::datatable(rownames = FALSE, escape = F,
                    options = list(
                      pageLength = 30,
                      autoWidth = T,
                      columnDefs = list(list(width = '30%', targets = 0))
                    ))
  })

  output$gsea_table_mf <- DT::renderDataTable({
    readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/gsea_res.rds")) %>%
      filter(gs == "MF") %>%
      pull(data) %>% purrr::reduce(~ .x[[1]]) %>%
      filter(parameter == input$gsea_res_parameter) %>%
      select(-parameter) %>%
      DT::datatable(rownames = FALSE, escape = F,
                    options = list(
                      pageLength = 30,
                      autoWidth = T,
                      columnDefs = list(list(width = '30%', targets = 0))
                    ))
  })

  output$space_anatomic_region <- DT::renderDataTable({
    readr::read_csv(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/table_s5.csv")) %>%
      DT::datatable(rownames = FALSE, options = list(pageLength = 30)) %>%
      formatRound(4:6, digits = 3) %>%
      formatSignif(7:8)
  })

  output$bioreplicate_g <-
    renderPlot(
      {
        t <-
          readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/bioreplicate_t.rds")) %>%
          select(biological_rep_id,
                 subject, visit, skin_type, feature,
                 counts_scaled, replicate_ID) %>%
          arrange(biological_rep_id) %>%
          mutate(counts_scaled = counts_scaled + 1) %>%
          filter(feature %in% input$gene_to_plot_duplicate)  #, "CCL18", "IL10", "IL22"

        stat_test <-
          t %>% group_by(feature, skin_type) %>%
          t_test(counts_scaled ~ replicate_ID, paired = T) %>%
          add_significance() %>%
          add_xy_position(x = "replicate_ID",
                          y.trans = log2, step.increase = .08) %>%
          group_by(feature) %>% mutate(y.position = max(y.position)) %>%
          ungroup()

          t %>%
          pivot_wider(names_from = replicate_ID, values_from = counts_scaled) %>%
          ggpaired(cond1 = "01", cond2 = "02", y = "counts_scaled",
                   xlab = "Biological replicate",
                   ylab = "",
                   line.color = "grey", line.size = .4,
                   fill = "skin_type",
                   palette = c("LS" = "#eb2d0c", "NL" = "#eb8b9b", "HC" = "#91cf60")) %>%
          facet(facet.by = c("feature", "skin_type"), scales = "free_y") +
          stat_pvalue_manual(stat_test, label = "p") +
          scale_y_continuous(trans = "log2", expand = expansion(mult = c(.05, .1)))
        }
      # res = 96
    )

  output$bioreplicate.ui <- renderUI(
    {plotOutput("bioreplicate_g", width = "100%",
                height = paste0(200*length(input$gene_to_plot_duplicate), "px"))}
  )

  output$cell_type_variation_g <-
    renderPlot(
      readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/cell_type_deconvolution.rds")) %>%
        filter(subject %in% input$subject_to_plot) %>%
        ggbarplot(x = "plot_name", y = "prop", fill = "cell_type",
                  xlab = "Sample", ylab = "Proportion(%)", palette = "npg") %>%
        facet(facet.by = "subject", scales = "free", ncol = 2) %>%
        ggpar(xtickslab.rt = 90, font.tickslab = 6)
      # res = 96
    )

  output$celltypevariation.ui <- renderUI(
    {plotOutput("cell_type_variation_g", width = "100%",
                height = paste0(300*round( (length(input$subject_to_plot) + .01)/2), "px"))}
  )

  output$heatmap_linc <-
    renderPlot(
      ComplexHeatmap::draw(readRDS(
        url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/heatmap_linc.rds")), heatmap_legend_side = "right",
                           annotation_legend_side = "top", merge_legend = FALSE), height = 1500
    )

    # renderImage({
    #   list(
    #     src = file.path(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/figure_s5.jpg")),
    #     contentType = "image/jpg",
    #     width = 1000
    #   )
    # }, deleteFile = F)

  output$time_variation_g <-
    renderPlot(
      readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/time_variation_t.rds")) %>%
        filter(skin_type == input$skin_type, time_type == input$time) %>%
        pull(d) %>% purrr::reduce(~.x) %>%
        ggboxplot(x = ifelse(input$time == "quarter", "visit_quarter", "visit"),
                  y = "counts_scaled",
                  add = "jitter", add.params = list(alpha = .5),
                  facet.by = "feature",
                  ylab = FALSE,
                  yscale = "log10",
                  scales = "free",
                  ncol = 3)
    )
  output$timevariation.ui <- renderUI(
    {plotOutput("time_variation_g",
                width = ifelse(input$time == "quarter" & input$skin_type == "HC", "50%", "100%"),
                paste0(case_when(input$time == "quarter" & input$skin_type == "LS" ~ 33,
                                 input$time == "quarter" & input$skin_type == "NL" ~ 3,
                                 input$time == "quarter" & input$skin_type == "HC" ~ 1.3,
                                 input$time == "visit" & input$skin_type == "LS" ~ 12,
                                 input$time == "visit" & input$skin_type == "NL" ~ 24,
                                 input$time == "visit" & input$skin_type == "HC" ~ 8)*200,
                       "px")
                )}
                # height = paste0(200*5, "px"))}
  )
  output$RNA_conc_kable <- function(){
    readr::read_csv("../data/table_theo_calc_RNA_conc.csv") %>%
      knitr::kable("html") %>%
      kableExtra::kable_styling("striped", full_width = F) %>%
      kableExtra::add_footnote("Tsoi et al. 2019 has reported that from a 5 mm punch biopsy, 50 ng/uL RNA can be obtained for high quality sequencing.")
  }
  output$table_cell_composition_gse121212_genad <-
    DT::renderDataTable({
      readr::read_csv("../data/table_cell_composition_gse121212_genad.csv") %>%
      DT::datatable(rownames = FALSE)
      })

  output$tissue_injury_g <-
    renderPlot(
      readRDS("../data/tissue_injury_se_t.rds") %>%
        ggline(x = "visit",
               y = "counts_scaled",
               facet.by = "feature",
               add = c("mean", "dotplot"),
               color = "skin_type",
               palette = c("LS" = "#eb2d0c", "NL" = "#eb8b9b", "HC" = "#91cf60"),
               xlab = "Visit",
               ylab = "Counts") %>%
        ggpar(yscale = "log10")
    )

  output$subcuits_g <-
    renderPlot(
      readRDS("../data/subcutis_se_t.rds") %>%
        ggline(x = "visit",
               y = "counts_scaled",
               facet.by = "feature",
               add = c("mean", "dotplot"),
               color = "skin_type",
               palette = c("LS" = "#eb2d0c", "NL" = "#eb8b9b", "HC" = "#91cf60"),
               xlab = "Visit",
               ylab = "Counts") %>%
        ggpar(yscale = "log10")
    )
}

# Run the app ----
shinyApp(ui = ui, server = server)
