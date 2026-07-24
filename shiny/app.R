library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(tidyr)

# ---- Load data ----
# Expects med_all_res.RData in the same folder, containing a data frame `res`
load("med_all_res.RData")

# CV and deff columns follow a "<metric>_<target>" naming pattern, e.g. cv_y_1, deff_b_1
cv_cols   <- grep("^cv_",   names(res), value = TRUE)
deff_cols <- grep("^deff_", names(res), value = TRUE)
targets   <- sub("^cv_", "", cv_cols)  # e.g. "y_1", "y_2", "b_1", "g_1"

ui <- fluidPage(
  titlePanel("Mixed variable stratification results explorer"),

  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Method", choices = sort(unique(res$method)), multiple = TRUE,
                  selected = unique(res$method)),
      selectInput("dataset", "Dataset", choices = sort(unique(res$dataset)), multiple = TRUE,
                  selected = unique(res$dataset)),
      selectInput("type", "Type", choices = sort(unique(res$type)), multiple = TRUE,
                  selected = unique(res$type)),
      sliderInput("strata", "Number of strata",
                  min = min(res$strata), max = max(res$strata),
                  value = c(min(res$strata), max(res$strata)), step = 1),
      selectInput("target", "Target variable (for CV / deff plots)",
                  choices = targets, selected = targets[1]),
      hr(),
      downloadButton("download_filtered", "Download filtered data (.csv)")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("CV by strata",
                 plotOutput("cv_plot", height = "450px")),
        tabPanel("Design effect (deff) by strata",
                 plotOutput("deff_plot", height = "450px")),
        tabPanel("Sample size (n) by strata",
                 plotOutput("n_plot", height = "450px")),
        tabPanel("Summary table",
                 DTOutput("summary_table")),
        tabPanel("Raw filtered data",
                 DTOutput("raw_table"))
      )
    )
  )
)

server <- function(input, output, session) {

  filtered <- reactive({
    res %>%
      filter(
        method %in% input$method,
        dataset %in% input$dataset,
        type %in% input$type,
        strata >= input$strata[1],
        strata <= input$strata[2]
      )
  })

  cv_col_selected   <- reactive(paste0("cv_",   input$target))
  deff_col_selected <- reactive(paste0("deff_", input$target))

  cv_plot_data <- reactive({
    df <- filtered()
    col <- cv_col_selected()
    validate(need(col %in% names(df), paste("No CV column for target", input$target)))
    df %>%
      mutate(cv_val = .data[[col]]) %>%
      filter(!is.na(cv_val))
  })

  deff_plot_data <- reactive({
    df <- filtered()
    col <- deff_col_selected()
    validate(need(col %in% names(df), paste("No deff column for target", input$target)))
    df %>%
      mutate(deff_val = .data[[col]]) %>%
      filter(!is.na(deff_val))
  })

  output$cv_plot <- renderPlot({
    df <- cv_plot_data()
    validate(need(nrow(df) > 0, "No data for this combination of filters."))
    ggplot(df, aes(x = factor(strata), y = cv_val, colour = method)) +
      geom_boxplot(position = position_dodge(width = 0.7)) +
      facet_wrap(~ type + dataset) +
      labs(x = "Number of strata", y = paste("CV -", input$target),
           title = paste("Coefficient of variation across seeds -", input$target)) +
      theme_minimal(base_size = 13)
  })

  output$deff_plot <- renderPlot({
    df <- deff_plot_data()
    validate(need(nrow(df) > 0, "No data for this combination of filters."))
    ggplot(df, aes(x = factor(strata), y = deff_val, colour = method)) +
      geom_boxplot(position = position_dodge(width = 0.7)) +
      facet_wrap(~ type + dataset) +
      labs(x = "Number of strata", y = paste("Design effect -", input$target),
           title = paste("Design effect across seeds -", input$target)) +
      theme_minimal(base_size = 13)
  })

  output$n_plot <- renderPlot({
    df <- filtered()
    validate(need(nrow(df) > 0, "No data for this combination of filters."))
    ggplot(df, aes(x = factor(strata), y = n, colour = method)) +
      geom_boxplot(position = position_dodge(width = 0.7)) +
      facet_wrap(~ type + dataset) +
      labs(x = "Number of strata", y = "Required sample size (n)",
           title = "Sample size across seeds") +
      theme_minimal(base_size = 13)
  })

  output$summary_table <- renderDT({
    df <- filtered()
    validate(need(nrow(df) > 0, "No data for this combination of filters."))
    summary_df <- df %>%
      group_by(method, dataset, type, strata) %>%
      summarise(
        n_runs    = n(),
        mean_n    = mean(n, na.rm = TRUE),
        mean_deff = mean(deff, na.rm = TRUE),
        mean_time = mean(time, na.rm = TRUE),
        across(all_of(cv_cols), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
        .groups = "drop"
      )
    numeric_cols <- names(summary_df)[sapply(summary_df, is.numeric)]
    datatable(summary_df, options = list(pageLength = 15), rownames = FALSE) %>%
      formatRound(columns = numeric_cols, digits = 3)
  })

  output$raw_table <- renderDT({
    datatable(filtered(), options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })

  output$download_filtered <- downloadHandler(
    filename = function() paste0("filtered_results_", Sys.Date(), ".csv"),
    content = function(file) write.csv(filtered(), file, row.names = FALSE)
  )
}

shinyApp(ui, server)
