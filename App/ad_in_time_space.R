library(shiny)
library(DT)
library(dplyr)
library(ggpubr)
library(purrr)
library(rstatix)
library(tidyr)


# Define UI ----
ui <- navbarPage("AD in Time and Space",
                 tabPanel("Introduction",
                          p("This web interface facilitaes the readers to interact with and download the data published on Hu. et al, Journal of Investigative Dermatology"),
                          h3("Download"),
                          downloadLink("TimeVariation_download", label = "Time variation data"),
                          br(),
                          downloadLink("SpaceVariation_bioreplicate_download", label = "Space variation (biological replicate) data"), # if two widgets get same name. it will not work without a warning message
                          br(),
                          downloadLink("CellTypeVariation_download", label = "Cell type variation data"),
                          br(),
                          downloadLink("LINCRNA_download", label = "LINC RNA expression data"),
                          br(),
                          downloadLink("ADSignatureGene_download", label = "AD signature gene data"),
                          br(),
                          downloadLink("GSEAResult_download", label = "GSEA result data"),
                          br(),
                          downloadLink("variancePartition_download", label = "Variance Parition data"),
                          br(),
                          downloadLink("SpaceVariation_download", label = "Space variation (anatomic region) data"),
                          br()),
                 tabPanel("Time variation", # Figure S2
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel("",
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
                                plotOutput("time_variation_g", width = "100%", height = "800px")
                              )
                            )
                          )),
                 tabPanel("Space variation (duplicate)",
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel("",
                                  checkboxGroupInput("gene_to_plot_duplicate",
                                                     "Choose the genes to plot:",
                                                     choices = readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/bioreplicate_t.rds"))$feature %>% unique,
                                                     selected = readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/bioreplicate_t.rds"))$feature %>% head(5))
                                           ),
                              mainPanel(
                                # plotOutput("bioreplicate_g", height = paste0(800, "px"))
                                uiOutput("bioreplicate.ui")
                              )
                            )
                          )), # Figure S3
                 tabPanel("Cell type variation",
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel("",
                                           checkboxGroupInput(
                                             "subject_to_plot",
                                             "Choose the subject IDs to plot:",
                                             choices = readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/cell_type_deconvolution.rds")) %>% pull(subject) %>% unique(),
                                             selected = c("AD_02", "AD_06")
                                           )),
                              mainPanel(
                                plotOutput("cell_type_variation_g", height = "800px")
                              )
                            )
                          )), # Figure S4
                 tabPanel("LINC RNA",
                          basicPage(
                            h2("Heatmap for LINC RNA"),
                            br(),
                            imageOutput("heatmap_linc")
                          )), # Figure S5
                 tabPanel("AD signature gene",
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
                 tabPanel("varianceParition",
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
                          )) # Table S5)
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
      pull(data) %>% reduce(~ .x[[1]]) %>%
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
      pull(data) %>% reduce(~ .x[[1]]) %>%
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
      pull(data) %>% reduce(~ .x[[1]]) %>%
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
      pull(data) %>% reduce(~ .x[[1]]) %>%
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
        facet(facet.by = "subject", scales = "free") %>%
        ggpar(xtickslab.rt = 90, font.tickslab = 6)
      # res = 96
    )

  output$heatmap_linc <-
    renderImage({
      list(
        src = file.path(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/figure_s5.jpg")),
        contentType = "image/jpg",
        width = 1000
      )
    }, deleteFile = F)

  output$time_variation_g <-
    renderPlot(
      readRDS(url("https://cdn.jsdelivr.net/gh/tuhulab/Shiny_AD_time_space@master/data/time_variation_t.rds")) %>%
        filter(skin_type == input$skin_type,
               time_type == input$time) %>%
        pull(d) %>% reduce(~.x) %>%
        ggboxplot(x = ifelse(input$time == "quarter", "visit_quarter", "visit"),
                  y = "counts_scaled",
                  add = "jitter",
                  facet.by = "feature",
                  ylab = FALSE,
                  yscale = "log10",
                  scales = "free")

    )
}

# Run the app ----
shinyApp(ui = ui, server = server)
