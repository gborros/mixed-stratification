
library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(tidyr)

# ---- Load data ----
load("app_res.RData")
app_res <- res


# ---- Add display method labels ----
# Keep the original `method` variable unchanged, but create a labelled
# version for plots/tables.

app_res <- app_res %>%
  mutate(
    method_label = case_when(
      method == "GGA-BC" & type == "continuous" ~ "GGA-BC (cont)",
      method == "GGA-BC" & type == "atomic" ~ "GGA-BC (atomic)",
      TRUE ~ method
    )
  )


# ---- Detect variable columns dynamically ----
cv_cols_all   <- grep("^cv_", names(app_res), value = TRUE)
cv_fpc_cols   <- grep("_fpc$", cv_cols_all, value = TRUE)
cv_base_cols  <- setdiff(cv_cols_all, cv_fpc_cols)
deff_var_cols <- grep("^deff_", names(app_res), value = TRUE)

targets <- sub("^cv_", "", cv_base_cols)
has_fpc_variant <- length(cv_fpc_cols) > 0


# ---- Choices for filters ----
ref_method_choices <- sort(unique(app_res$ref_method))
method_choices     <- sort(unique(app_res$method_label))
fpc_choices        <- sort(unique(app_res$fpc))
strata_choices     <- sort(unique(app_res$strata))


# ---- UI ----
ui <- fluidPage(
  
  # ---- CSS ----
  tags$head(
    tags$style(HTML("
      html, body {
        height: 100%;
        margin: 0;
      }

      #page-container {
        min-height: 100vh;
        display: flex;
        flex-direction: column;
      }

      #content-wrap {
        flex: 1;
      }

      #footer {
        text-align: center;
        font-size: 12px;
        color: grey;
        padding: 10px;
        border-top: 1px solid #ddd;
        background-color: #fafafa;
      }

      .plot-title {
        margin-top: 35px;
        margin-bottom: 5px;
      }

      .plot-description {
        margin-bottom: 15px;
      }
    "))
  ),
  
  div(
    id = "page-container",
    
    # ---- MAIN CONTENT ----
    div(
      id = "content-wrap",
      
      # ---- Title ----
      div(
        style = "text-align: center; margin-bottom: 10px;",
        tags$h2(
          "Mixed variable stratified sampling method comparison on the 2024 GHS Data"
        ),
        tags$p(
          "Georgi Borros, Şebnem Er, Sulaiman Salau",
          style = "font-size: 15px; color: #555; margin-top: -10px;"
        )
      ),
      
      sidebarLayout(
        
        # ---- SIDEBAR ----
        sidebarPanel(
          
          selectInput(
            "method",
            "Method",
            choices = method_choices,
            multiple = TRUE,
            selected = method_choices
          ),
          
          selectInput(
            "ref_method",
            "Reference method",
            choices = ref_method_choices,
            multiple = TRUE,
            selected = ref_method_choices
          ),
          
          selectInput(
            "fpc",
            "FPC status",
            choices = fpc_choices,
            multiple = TRUE,
            selected = fpc_choices
          ),
          
          selectInput(
            "strata",
            "Number of strata",
            choices = strata_choices,
            multiple = TRUE,
            selected = strata_choices
          ),
          
          # ---- Target variable only shown on variable-specific page ----
          conditionalPanel(
            condition = "input.main_tabs == 'variable_specific'",
            
            hr(),
            
            selectInput(
              "target",
              "Target variable",
              choices = targets,
              selected = targets[1]
            )
          ),
          
          hr(),
          
          if (has_fpc_variant) {
            radioButtons(
              "cv_variant",
              "CV variant for totals/plots",
              choices = c(
                "Standard" = "base",
                "FPC-adjusted" = "fpc"
              ),
              selected = "base"
            )
          },
          
          hr(),
          
          downloadButton(
            "download_filtered",
            "Download filtered data (.csv)"
          )
        ),
        
        
        # ---- MAIN PANEL ----
        mainPanel(
          width = 8,
          
          tabsetPanel(
            id = "main_tabs",
            
            # ==========================================================
            # PAGE 1: OVERALL RESULTS
            # ==========================================================
            tabPanel(
              title = "Overall results",
              value = "overall",
              icon = icon("chart-line"),
              
              # ---- Overall Deff ----
              tags$div(
                class = "plot-title",
                tags$h4("Deff by method")
              ),
              
              tags$p(
                class = "plot-description",
                "Distribution of the overall design effect across seeds, by method and number of strata."
              ),
              
              plotOutput(
                "deff_plot",
                height = "500px"
              ),
              
              hr(),
              
              # ---- Total CV ----
              tags$div(
                class = "plot-title",
                tags$h4("Total CV by method")
              ),
              
              tags$p(
                class = "plot-description",
                "Distribution of the total coefficient of variation across seeds, by method and number of strata."
              ),
              
              plotOutput(
                "total_cv_plot",
                height = "500px"
              ),
              
              hr(),
              
              # ---- Total n ----
              tags$div(
                class = "plot-title",
                tags$h4("Total sample size by method")
              ),
              
              tags$p(
                class = "plot-description",
                "Distribution of the total sample size across seeds, by method and number of strata."
              ),
              
              plotOutput(
                "total_n_plot",
                height = "500px"
              )
            ),
            
            
            # ==========================================================
            # PAGE 2: VARIABLE-SPECIFIC
            # ==========================================================
            tabPanel(
              title = "Variable-specific",
              value = "variable_specific",
              icon = icon("chart-bar"),
              
              tags$div(
                class = "plot-title",
                tags$h4("Variable-specific CV")
              ),
              
              tags$p(
                class = "plot-description",
                "Distribution of the coefficient of variation for the selected target variable, by method and number of strata."
              ),
              
              plotOutput(
                "var_cv_plot",
                height = "500px"
              ),
              
              hr(),
              
              tags$div(
                class = "plot-title",
                tags$h4("Variable-specific Deff")
              ),
              
              tags$p(
                class = "plot-description",
                "Distribution of the design effect for the selected target variable, by method and number of strata."
              ),
              
              plotOutput(
                "var_deff_plot",
                height = "500px"
              )
            ),
            
            
            # ==========================================================
            # SUMMARY TABLE
            # ==========================================================
            tabPanel(
              title = "Summary table",
              value = "summary",
              icon = icon("table"),
              
              tags$p(
                "Aggregated across seeds for the current filter selection.",
                style = "margin-top: 15px;"
              ),
              
              DTOutput("summary_table")
            ),
            
            
            # ==========================================================
            # RAW DATA
            # ==========================================================
            tabPanel(
              title = "Raw filtered data",
              value = "raw",
              icon = icon("list"),
              
              DTOutput("raw_table")
            )
          )
        )
      )
    ),
    
    
    # ---- FOOTER ----
    div(
      id = "footer",
      paste(
        "Last updated:",
        format(Sys.time(), "%d %B %Y, %H:%M")
      )
    )
  )
)


# =====================================================================
# SERVER
# =====================================================================

server <- function(input, output, session) {
  
  # ---- Which CV columns to use ----
  active_cv_cols <- reactive({
    
    if (
      has_fpc_variant &&
      isTRUE(input$cv_variant == "fpc")
    ) {
      cv_fpc_cols
    } else {
      cv_base_cols
    }
    
  })
  
  
  # ---- Apply filters ----
  filtered <- reactive({
    
    df <- app_res %>%
      filter(
        method_label %in% input$method,
        ref_method %in% input$ref_method,
        fpc %in% input$fpc,
        strata %in% input$strata
      )
    
    # Calculate total CV
    cols <- active_cv_cols()
    
    df$total_cv <- rowSums(
      df[, cols, drop = FALSE],
      na.rm = TRUE
    )
    
    df
  })
  
  
  # ---- Selected CV column ----
  cv_col_selected <- reactive({
    
    cols <- active_cv_cols()
    
    hit <- cols[
      sub(
        "_fpc$",
        "",
        sub("^cv_", "", cols)
      ) == input$target
    ]
    
    if (length(hit) == 0) {
      NA_character_
    } else {
      hit[1]
    }
    
  })
  
  
  # ---- Selected Deff column ----
  deff_col_selected <- reactive({
    
    paste0(
      "deff_",
      input$target
    )
    
  })
  
  
  # ===================================================================
  # OVERALL DEFF
  # ===================================================================
  
  output$deff_plot <- renderPlot({
    
    df <- filtered()
    
    validate(
      need(
        nrow(df) > 0,
        "No data for this combination of filters."
      )
    )
    
    ggplot(
      df,
      aes(
        x = factor(seed),
        y = deff,
        colour = method_label
      )
    ) +
      
      geom_boxplot(
        position = position_dodge(width = 0.7)
      ) +
      
      scale_colour_brewer(
        palette = "Set1"
      ) +
      
      labs(
        x = "Number of strata",
        y = "Design effect (overall)",
        colour = "Method",
        title = "Overall design effect by method"
      ) +
      
      theme_bw(
        base_size = 15
      )
    
  })
  
  
  # ===================================================================
  # TOTAL CV
  # ===================================================================
  
  output$total_cv_plot <- renderPlot({
    
    df <- filtered()
    
    validate(
      need(
        nrow(df) > 0,
        "No data for this combination of filters."
      )
    )
    
    variant_label <-
      if (
        has_fpc_variant &&
        isTRUE(input$cv_variant == "fpc")
      ) {
        "FPC-adjusted"
      } else {
        "Standard"
      }
    
    ggplot(
      df,
      aes(
        x = factor(seed),
        y = total_cv,
        colour = method_label
      )
    ) +
      
      geom_boxplot(
        position = position_dodge(width = 0.7)
      ) +
      
      scale_colour_brewer(
        palette = "Set1"
      ) +
      
      labs(
        x = "Number of strata",
        y = paste(
          "Total CV (summed,",
          variant_label,
          ")"
        ),
        colour = "Method",
        title = paste(
          "Total CV across target variables by method -",
          variant_label
        )
      ) +
      
      theme_bw(
        base_size = 15
      )
    
  })
  
  
  # ===================================================================
  # TOTAL SAMPLE SIZE
  # ===================================================================
  
  output$total_n_plot <- renderPlot({
    
    df <- filtered()
    
    validate(
      need(
        nrow(df) > 0,
        "No data for this combination of filters."
      )
    )
    
    ggplot(
      df,
      aes(
        x = factor(seed),
        y = n,
        colour = method_label
      )
    ) +
      
      geom_boxplot(
        position = position_dodge(width = 0.7)
      ) +
      
      scale_colour_brewer(
        palette = "Set1"
      ) +
      
      labs(
        x = "Number of strata",
        y = "Total sample size",
        colour = "Method",
        title = "Total sample size by method"
      ) +
      
      theme_bw(
        base_size = 15
      )
    
  })
  
  
  # ===================================================================
  # VARIABLE-SPECIFIC CV
  # ===================================================================
  
  output$var_cv_plot <- renderPlot({
    
    df <- filtered()
    
    col <- cv_col_selected()
    
    validate(
      need(
        !is.na(col),
        paste(
          "No CV column found for target",
          input$target,
          "in the selected variant."
        )
      )
    )
    
    validate(
      need(
        nrow(df) > 0,
        "No data for this combination of filters."
      )
    )
    
    df <- df %>%
      mutate(
        cv_val = .data[[col]]
      ) %>%
      filter(
        !is.na(cv_val)
      )
    
    ggplot(
      df,
      aes(
        x = factor(strata),
        y = cv_val,
        colour = method_label
      )
    ) +
      
      geom_boxplot(
        position = position_dodge(width = 0.7)
      ) +
      
      scale_colour_brewer(
        palette = "Set1"
      ) +
      
      labs(
        x = "Number of strata",
        y = paste(
          "CV -",
          input$target
        ),
        colour = "Method",
        title = paste(
          "CV by method -",
          input$target
        )
      ) +
      
      theme_bw(
        base_size = 15
      )
    
  })
  
  
  # ===================================================================
  # VARIABLE-SPECIFIC DEFF
  # ===================================================================
  
  output$var_deff_plot <- renderPlot({
    
    df <- filtered()
    
    col <- deff_col_selected()
    
    validate(
      need(
        col %in% names(df),
        paste(
          "No deff column for target",
          input$target
        )
      )
    )
    
    validate(
      need(
        nrow(df) > 0,
        "No data for this combination of filters."
      )
    )
    
    df <- df %>%
      mutate(
        deff_val = .data[[col]]
      ) %>%
      filter(
        !is.na(deff_val)
      )
    
    ggplot(
      df,
      aes(
        x = factor(strata),
        y = deff_val,
        colour = method_label
      )
    ) +
      
      geom_boxplot(
        position = position_dodge(width = 0.7)
      ) +
      
      scale_colour_brewer(
        palette = "Set1"
      ) +
      
      labs(
        x = "Number of strata",
        y = paste(
          "Deff -",
          input$target
        ),
        colour = "Method",
        title = paste(
          "Design effect by method -",
          input$target
        )
      ) +
      
      theme_bw(
        base_size = 15
      )
    
  })
  
  
  # ===================================================================
  # SUMMARY TABLE
  # ===================================================================
  
  output$summary_table <- renderDT({
    
    df <- filtered()
    
    validate(
      need(
        nrow(df) > 0,
        "No data for this combination of filters."
      )
    )
    
    cols <- active_cv_cols()
    
    summary_df <- df %>%
      
      group_by(
        method_label,
        ref_method,
        fpc,
        strata
      ) %>%
      
      summarise(
        
        n_runs = n(),
        
        mean_n =
          mean(
            n,
            na.rm = TRUE
          ),
        
        mean_deff =
          mean(
            deff,
            na.rm = TRUE
          ),
        
        mean_total_cv =
          mean(
            total_cv,
            na.rm = TRUE
          ),
        
        mean_time =
          mean(
            time,
            na.rm = TRUE
          ),
        
        across(
          all_of(cols),
          ~ mean(
            .x,
            na.rm = TRUE
          ),
          .names = "mean_{.col}"
        ),
        
        across(
          all_of(deff_var_cols),
          ~ mean(
            .x,
            na.rm = TRUE
          ),
          .names = "mean_{.col}"
        ),
        
        .groups = "drop"
      )
    
    numeric_cols <-
      names(summary_df)[
        sapply(
          summary_df,
          is.numeric
        )
      ]
    
    datatable(
      summary_df,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      ),
      rownames = FALSE
    ) %>%
      
      formatRound(
        columns = numeric_cols,
        digits = 3
      )
    
  })
  
  
  # ===================================================================
  # RAW TABLE
  # ===================================================================
  
  output$raw_table <- renderDT({
    
    datatable(
      filtered(),
      options = list(
        pageLength = 15,
        scrollX = TRUE
      ),
      rownames = FALSE
    )
    
  })
  
  
  # ===================================================================
  # DOWNLOAD
  # ===================================================================
  
  output$download_filtered <- downloadHandler(
    
    filename = function() {
      paste0(
        "filtered_results_",
        Sys.Date(),
        ".csv"
      )
    },
    
    content = function(file) {
      write.csv(
        filtered(),
        file,
        row.names = FALSE
      )
    }
    
  )
  
}


# ---- Run app ----
shinyApp(
  ui,
  server
)

