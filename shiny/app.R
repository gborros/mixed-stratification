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
# version for plots/tables that folds the reference method into the
# label for GGA-BC and med-ga rows.

app_res <- app_res %>%
  mutate(
    method_label = case_when(
      method == "GGA-BC" & ref_method == "manual"  ~ "GGA-BC (manual)",
      method == "GGA-BC" & ref_method == "med-ga"   ~ "GGA-BC (MED-GA)",
      method == "MED-GA" & ref_method == "gga-bc"   ~ "MED-GA (GGA-BC)",
      method == "MED-GA" & ref_method == "bc"       ~ "MED-GA (BC)",
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
method_choices  <- sort(unique(app_res$method_label))
fpc_choices     <- sort(unique(app_res$fpc))
type_choices    <- sort(unique(app_res$type))

# Default selection excludes all GGA-BC variants
method_choices_default <- method_choices[!grepl("^GGA-BC", method_choices)]

# Default FPC selection = "no fpc"
# Falls back to the first available option if "no fpc" is not present.
fpc_default <- if ("no fpc" %in% fpc_choices) {
  "no fpc"
} else {
  fpc_choices[1]
}


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
      
      .method-explainer {
  background-color: #f7f7f7;
  border: 1px solid #ddd;
  border-radius: 5px;
  padding: 12px 15px;
  margin-bottom: 20px;
  font-size: 13px;
  line-height: 1.5;
}

.method-explainer h4 {
  margin-top: 0;
  margin-bottom: 8px;
  font-size: 15px;
}

.method-explainer p {
  margin-bottom: 6px;
}

.method-explainer strong {
  font-weight: 600;
}

.deff-warning {
  background-color: #fff3cd;
  border: 1px solid #ffeeba;
  color: #856404;
  border-radius: 5px;
  padding: 10px 15px;
  margin-bottom: 10px;
  font-size: 13px;
  line-height: 1.5;
}

.deff-warning strong {
  font-weight: 600;
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
          # ---- Method explainer ----
          tags$div(
            class = "method-explainer",
            
            tags$h4("Understanding the method labels"),
            
            tags$p(
              "The method in parentheses indicates the method used to generate ",
              "the inputs for the method being evaluated."
            ),
            
            tags$p(
              tags$strong("GGA-BC (manual): "),
              "GGA-BC results using manually specified CVs (0.05)."
            ),
            
            tags$p(
              tags$strong("MED-GA (GGA-BC): "),
              "MED-GA results using GGA-BC sample size and strata as inputs."
            ),
            
            tags$p(
              tags$strong("GGA-BC (MED-GA): "),
              "GGA-BC results using MED-GA CVs as inputs."
            ),
            
            tags$p(
              tags$strong("BC: "),
              "Bethel-Chromy results using MED-GA CVs and stratification."
            ),
            
            tags$p(
              tags$strong("MED-GA (BC): "),
              "MED-GA results using BC sample sizes as inputs."
            )
          ),
          
          selectInput(
            "method",
            "Method",
            choices = method_choices,
            multiple = TRUE,
            selected = method_choices_default
          ),
          
          selectInput(
            "fpc",
            "FPC status",
            choices = fpc_choices,
            multiple = FALSE,
            selected = fpc_default
          ),
          
          selectInput(
            "type",
            "Type",
            choices = type_choices,
            multiple = FALSE,
            selected = type_choices
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
              
              # ---- Total n ----
              tags$div(
                class = "plot-title",
                tags$h4("Total sample size by method")
              ),
              
              tags$p(
                class = "plot-description",
                "Sample size for each seed, by method."
              ),
              
              plotOutput(
                "total_n_plot",
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
                "Total coefficient of variation for each seed, by method."
              ),
              
              plotOutput(
                "total_cv_plot",
                height = "500px"
              ),
              
              hr(),
              
              # ---- Overall Deff ----
              tags$div(
                class = "plot-title",
                tags$h4("Deff by method")
              ),
              
              tags$p(
                class = "plot-description",
                "Overall design effect for each seed, by method."
              ),
              
              uiOutput("deff_warning_overall"),
              
              plotOutput(
                "deff_plot",
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
                "Coefficient of variation for each seed, by method."
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
                "Design effect for each seed, by method."
              ),
              
              uiOutput("deff_warning_var"),
              
              plotOutput(
                "var_deff_plot",
                height = "500px"
              )
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
        fpc %in% input$fpc,
        type %in% input$type
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
  
  
  # ---- Whether GGA-BC (MED-GA) + atomic is present in filtered data ----
  # This method approximates a full census, so the SRS design variance used
  # in the DEFF denominator becomes very close to zero (e.g. dividing something
  # like 0.0003 by 0.000006), which makes DEFF blow up and become unstable
  # rather than reflecting a genuine efficiency loss.
  deff_warning_needed <- reactive({
    
    df <- filtered()
    
    if (nrow(df) == 0) {
      return(FALSE)
    }
    
    any(
      df$method_label == "GGA-BC (MED-GA)" &
        df$type == "atomic"
    )
    
  })
  
  
  deff_warning_ui <- function() {
    tags$div(
      class = "deff-warning",
      tags$strong("Note: "),
      "Since GGA-BC (MED-GA) approximates a full census, the SRS design ",
      "variance becomes very close to 0, making DEFF unstable."
    )
  }
  
  
  output$deff_warning_overall <- renderUI({
    
    if (isTRUE(deff_warning_needed())) {
      deff_warning_ui()
    }
    
  })
  
  
  output$deff_warning_var <- renderUI({
    
    if (isTRUE(deff_warning_needed())) {
      deff_warning_ui()
    }
    
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
      
      geom_point(
        position = position_dodge(width = 0.1),
        size = 4
      ) +
      
      scale_colour_brewer(
        palette = "Set2"
      ) +
      
      labs(
        x = "Seed",
        y = "Design effect (overall)",
        colour = "Method",
        title = "Overall design effect by method"
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
      
      geom_point(
        position = position_dodge(width = 0.1),
        size = 4
      ) +
      
      scale_colour_brewer(
        palette = "Set2"
      ) +
      
      labs(
        x = "Seed",
        y = "Total sample size",
        colour = "Method",
        title = "Total sample size by method"
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
      
      geom_point(
        position = position_dodge(width = 0.1),
        size = 4
      ) +
      
      scale_colour_brewer(
        palette = "Set2"
      ) +
      
      labs(
        x = "Seed",
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
        x = factor(seed),
        y = cv_val,
        colour = method_label
      )
    ) +
      
      geom_point(
        position = position_dodge(width = 0.1),
        size = 4
      ) +
      
      scale_colour_brewer(
        palette = "Set2"
      ) +
      
      labs(
        x = "Seed",
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
        x = factor(seed),
        y = deff_val,
        colour = method_label
      )
    ) +
      
      geom_point(
        position = position_dodge(width = 0.1),
        size = 4
      ) +
      
      scale_colour_brewer(
        palette = "Set2"
      ) +
      
      labs(
        x = "Seed",
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