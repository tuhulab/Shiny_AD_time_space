library(shiny)
library(DT)
library(dplyr)
library(ggpubr)
library(purrr)


# Define UI ----
ui <- navbarPage("AD in Time and Space",
                 tabPanel("Introduction",
                          "This applications facilitaes the readers to interact with the data published on Hu. et al, Journal of Investigative Dermatology"),
                 tabPanel("Time variation",
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
                                plotOutput("time_variation_g")
                              )
                            )
                          )), # Figure S2
                 tabPanel("Space variation (duplicate)",
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel("",
                                  checkboxGroupInput("gene_to_plot_duplicate",
                                                     "Choose the genes to plot:",
                                                     choices = readRDS(url("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/bioreplicate_t.rds"))$feature %>% unique,
                                                     selected = readRDS(url("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/bioreplicate_t.rds"))$feature %>% head(5))
                                           ),
                              mainPanel(
                                plotOutput("bioreplicate_g")
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
                                             choices = readr::read_rds("../data/cell_type_deconvolution.rds") %>% pull(subject) %>% unique(),
                                             selected = c("AD_02", "AD_06")
                                           )),
                              mainPanel(
                                plotOutput("cell_type_variation_g")
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
    readr::read_csv("../data/table_s4.csv") %>%
      dplyr::select(-1) %>%
      arrange(-`tissue type (LS/NL/HC)`) %>%
      DT::datatable(rownames = FALSE, options = list(pageLength = 30)) %>%
      formatPercentage(2:9)
  })

  output$ad_signature_table <-
    DT::renderDataTable({
      readr::read_csv("../data/table_s2.csv") %>%
        DT::datatable(rownames = FALSE, options = list(pageLength = 30)) %>%
        formatRound(2:4, digits = 3) %>%
        formatSignif(5:6)
    })

  output$gsea_table_reactome <- DT::renderDataTable({
    readRDS(url("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/gsea_res.rds")) %>%
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
    readRDS("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/gsea_res.rds") %>%
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
    readRDS("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/gsea_res.rds") %>%
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
    readRDS("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/gsea_res.rds") %>%
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
    readr::read_csv("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/table_s5.csv") %>%
      DT::datatable(rownames = FALSE, options = list(pageLength = 30)) %>%
      formatRound(4:6, digits = 3) %>%
      formatSignif(7:8)
  })

  output$bioreplicate_g <-
    renderPlot(
      readr::read_rds("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/bioreplicate_t.rds") %>%
        select(biological_rep_id,
               subject, visit, skin_type, feature,
               counts_scaled, replicate_ID) %>%
        filter(feature %in% input$gene_to_plot_duplicate) %>%
        ggplot(aes(x = replicate_ID,
                   y = counts_scaled,
                   group = biological_rep_id)) +
        geom_point(aes(color = skin_type)) +
        facet_grid(feature ~ skin_type, scales = "free") +
        geom_line(aes(color = skin_type), alpha = .5) +
        geom_boxplot(aes(group = replicate_ID, fill = skin_type), alpha = .5) +
        scale_color_manual(values = c("LS" = "#eb2d0c", "NL" = "#eb8b9b", "HC" = "#91cf60")) +
        scale_fill_manual(values = c("LS" = "#eb2d0c", "NL" = "#eb8b9b", "HC" = "#91cf60")) +
        scale_y_continuous(trans = "log2") +
        xlab("replicate") +
        theme(axis.title.y = element_blank())
      # res = 96
    )

  output$cell_type_variation_g <-
    renderPlot(
      readr::read_rds("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/cell_type_deconvolution.rds") %>%
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
        src = file.path("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/figure_s5.jpg"),
        contentType = "image/jpg",
        width = 1000
      )
    }, deleteFile = F)

  output$time_variation_g <-
    renderPlot(
      readr::read_rds("https://rawcdn.githack.com/tuhulab/Shiny_AD_time_space/4fec4d74e2b31b40d936c7abb49441a3f3e7ca10/data/time_variation_t.rds") %>%
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