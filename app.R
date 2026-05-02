# Import Libraries
# Install and load the shinydashboard package, which simultaneously loads shiny.
#install.packages("shinydashboard")
library(shiny)
library(bslib)
library(shinydashboard)
#library(gt)
library(DT)
library(dplyr)
library(tidyr)
library(plotly)
#library(ggplot2)
#library(tseries)
library(readr)
#library(MVN)
#########################

# ---- Chargement des fonctions du projet ----
# Toutes les fonctions statistiques et utilitaires sont dans le dossier R/
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(f, encoding = "UTF-8")
}


##############################################################

# Set Up UI Components
header <- dashboardHeader(title = "Rbcm", 
                          titleWidth = 200) 

sidebar <- dashboardSidebar(
  width = 200,
  sidebarMenu(
    menuItem("Main Dashboard", tabName = "dashboard", icon = icon("gauge-high")),
    menuItem("Methodological Framework", tabName = "methods", icon = icon("book-open-reader")),
    
    menuItem("Exploratory Analysis", icon = icon("searchengin"),
             menuSubItem("General Overview", tabName = "overview", icon = icon("eye")),
             menuSubItem("Statistical Summaries", tabName = "stat_summary", icon = icon("list-check")),
             menuSubItem("Visual Diagnostics", tabName = "graphs", icon = icon("chart-pie"))
    ),
    
    menuItem("Causal mediation", tabName = "mediate", icon = icon("diagram-project")),
    menuItem("External Link", href = "https://github.com/komiayi/dna_mediation/blob/main/docs/document_final.pdf", icon = icon("external-link")),
    
    # Search Form: For enhanced user interaction
    sidebarSearchForm(textId = "search", buttonId = "searchButton", label = "Search...")
  )
)

body <- dashboardBody(
  withMathJax(),
  tags$script(HTML("
    Shiny.addCustomMessageHandler('mathjax-typeset', function(message) {
      if (window.MathJax) {
        MathJax.Hub.Queue(['Typeset', MathJax.Hub]);
      }
    });
  ")),
  tabItems(
    tabItem(tabName = "methods",
            withMathJax(),
            uiOutput("markdown_ui")
    ),
    tabItem(tabName = "overview",
            h2("Descriptive analysis tab content"),
      
            fluidRow(
              column(6, fileInput("data_file", 
                                  "Data source acquisition (CSV or RData)", 
                                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".RData")
                        )
              ),
              column(3, 
                     actionButton("load_data", 
                                  "Initialize & Process Data", 
                                  class = "btn-info", style = "width: 100%; margin-bottom: 20px;"
                     )
              )
            ),
            fluidRow(
              column(3, selectInput("treat", "Exposure", choices = NULL)),
              column(3, selectInput("outcome", "Outcome", choices = NULL)),
              column(3, 
                     selectInput("mediators", "Mediators", choices = NULL, multiple = TRUE),
                     tags$small(style = "color: grey;", 
                                "Please select mediators in order: 1st = Primary, 2nd = Secondary.")
              )
            ),
            br(),
            fluidRow(
              column(3, selectInput("out_cov",
                                    "Outcome covariates",
                                    choices = NULL, 
                                    multiple = TRUE)
              ),
              column(3, selectInput("intmed_cov", 
                                    "Primary mediator covariates",
                                    choices = NULL, 
                                    multiple = TRUE
                        )
              ),
              column(3, selectInput("sed_cov", 
                                    "Secondary mediator covariates",
                                    choices = NULL, 
                                    multiple = TRUE
                        )
              )
            ),
            br(),
            fluidRow(
              column(3, numericInput("B", "Number of Bootstrap Resamples", value = 20)),
              column(6, radioButtons("interaction", "Interaction terms in the outcome model ?", 
                                     choices = list("Yes" = TRUE, "No" = FALSE), selected = FALSE)),
              
              conditionalPanel(
                condition = "input.interaction == 'TRUE'",
                column(12,
                       checkboxGroupInput("interaction_types", "Select interaction terms for the Outcome model:",
                                          choices = c(
                                            "Exposure \u00D7 Primary Mediator" = "treat_m1",
                                            "Exposure \u00D7 Secondary Mediator" = "treat_m2",
                                            "Primary Mediator \u00D7 Secondary Mediator" = "m1_m2",
                                            "Exposure \u00D7 Primary \u00D7 Secondary Mediator" = "treat_m1_m2"
                                          )),
                       tags$small(style = "color: grey;", 
                                  "Note: These interactions will be included in the regression model where the Outcome is the dependent variable.")
                )
              )
            ),
            
            fluidRow(
              br(),
              column(12, actionButton("run", "Run", class = "btn-primary")),
              verbatimTextOutput("debug_output")  
            )
    ),
    tabItem(tabName = "stat_summary",
            h2("Statistical summary"),
            fluidRow(
              box(title = "Data preview", width = 12, dataTableOutput("data_table")),
              box(title = "Residuals preview", width = 12, dataTableOutput("residual_table")),
              box(title = "Descriptive statistics for quantitative variables", width = 12, 
                  dataTableOutput("summary_stats_quant")),
              box(title = "Descriptive statistics for qualitative variables", width = 12,
                  dataTableOutput("summary_stats_qual"))
            )
    ),
    tabItem(tabName = "graphs",
            h2("Data Visualization and statistical hypothesis testing  of residuals "),
            
            fluidRow(
              tabBox(
                title = "Boxplot",
                side = "right", height = "450px",
                selected = textOutput("Outcomebox"),
                tabPanel(textOutput("Outcomebox"), plotlyOutput("boxplotoutcome")),
                tabPanel(textOutput("Primarymediatorbox"), plotlyOutput("boxplotmedint")),
                tabPanel(textOutput("Secondarymediatorbox"), plotlyOutput("boxplotmedsed"))
              ),
              tabBox(
                title = "Histogram",
                side = "right", height = "450px",
                selected = textOutput("Outcomehist"),
                tabPanel(textOutput("Outcomehist"), plotlyOutput("histogramoutcome")),
                tabPanel(textOutput("Primarymediatorhist"), plotlyOutput("histogrammedint")),
                tabPanel(textOutput("Secondarymediatorhist"), plotlyOutput("histogrammedsed"))
              )
            ),
            
            fluidRow(
              tabBox(
                title = "Normal Q-Q Plot",
                side = "right", height = "450px",
                selected = textOutput("Outcomeqq"),
                tabPanel(textOutput("Outcomeqq"), plotlyOutput("qqplotoutcome")),
                tabPanel(textOutput("Primarymediatorqq"), plotlyOutput("qqplotmedint")),
                tabPanel(textOutput("Secondarymediatorqq"), plotlyOutput("qqplotmedsed"))
              ),
              tabBox(
                title = "Normality test",
                side = "right", height = "350px",
                selected = textOutput("Outcomenor"),
                tabPanel(textOutput("Outcomenorm"), tableOutput("testoutcome")),
                tabPanel(textOutput("Primarymediatornorm"), tableOutput("testmedint")),
                tabPanel(textOutput("Secondarymediatornorm"), tableOutput("testmedsed")),
                tabPanel(textOutput("Multivariatenorm"), tableOutput("testmulti"))
              )
            )
    ),
    tabItem(tabName = "mediate",
            div(class = "d-flex justify-content-between align-items-center mb-4",
                h2(tags$i(class = "fas fa-chart-line me-2"), "Mediations results")
            ),
            fluidRow(
              column(width = 12,
                     card(
                       full_screen = TRUE,
                       card_header(
                         # Navigation par onglets (Pills pour un look moderne)
                         navset_card_pill(
                           id = "mediation_tabs",
                           nav_panel(
                             title = "Constant correlation (CC)",
                             div(class = "p-3",
                                 hr(),
                                 layout_column_wrap(
                                   width = 1/3, 
                                   heights_equal = "row",
                                   
                                   # direct
                                   card(
                                     card_header(
                                       div(class = "d-flex justify-content-center align-items-center",
                                           tags$b(withMathJax("\\(\\hat{\\rho}(0,1)\\)"))
                                       )
                                     ),
                                     uiOutput("cc_rho")
                                   ),
                                   # direct
                                   card(
                                     card_header(
                                       div(class = "d-flex justify-content-center align-items-center",
                                           tags$b(withMathJax("\\(\\hat{\\zeta}\\) (direct)"))
                                       )
                                     ),
                                     uiOutput("cc_direct_grid")
                                   ),
                                   # indirect
                                   card(
                                     card_header(
                                       div(class = "d-flex justify-content-center align-items-center",
                                           tags$b(withMathJax("\\(\\hat{\\delta}\\) (indirect)"))#,
                                           #bsicons::bs_icon("info-circle", title = "Each card represents a different correlation (rho) assumption")
                                       )
                                     ),
                                     uiOutput("cc_indirect_grid"),
                                   )
                                 )
                             )
                           ),
                           nav_panel(
                             title = "Non constant correlation (CNC)",
                             div(class = "p-3",
                                 hr(),
                                 layout_column_wrap(
                                   width = 1/3, 
                                   heights_equal = "row",
                                   
                                   # direct
                                   card(
                                     card_header(
                                       div(class = "d-flex justify-content-center align-items-center",
                                           tags$b(withMathJax("\\(\\hat{\\rho}(0,1)\\)"))
                                       )
                                     ),
                                     uiOutput("rho")
                                   ),
                                   # direct
                                   card(
                                     card_header(
                                       div(class = "d-flex justify-content-center align-items-center",
                                           tags$b(withMathJax("\\(\\hat{\\zeta}\\) (direct)"))
                                       )
                                     ),
                                     uiOutput("cnc_direct_grid")
                                   ),
                                   # indirect
                                   card(
                                     card_header(
                                       div(class = "d-flex justify-content-center align-items-center",
                                           tags$b(withMathJax("\\(\\hat{\\delta}\\) (indirect)"))#,
                                           #bsicons::bs_icon("info-circle", title = "Each card represents a different correlation (rho) assumption")
                                       )
                                     ),
                                      uiOutput("cnc_indirect_grid"),
                                   )
                                 )
                             )
                           )
                         )
                       ),
                       br(),br(),
                       # Pied de page informatif (Optionnel)
                       card_footer(
                         tags$small(class = "text-muted",
                               "Note : Results are based on bootstrap resampling (B = ", textOutput("current_B", inline = TRUE)," iterations)")
                       )
                     )
              )
            )
    )
  )
)  # Main body for content

# Assemble UI
ui <- dashboardPage(header, sidebar, body)

# Serveur
server <- function(input, output, session) {
  
  output$markdown_ui <- renderUI({
    tagList(
      includeHTML("docs/Methods.Rhtml"),#chapter_translation.Rhtml #Methodes.html
      tags$script(HTML("
      setTimeout(function() {
        if (window.MathJax && window.MathJax.Hub) {
          MathJax.Hub.Queue(['Typeset', MathJax.Hub]);
        }
      }, 100);
    "))
    )
  })
  
  # Variable pour stocker les données (persiste pendant la session)
  data <- reactiveVal(NULL)
  
  observeEvent(data(), {
    cols <- colnames(data())
    
    updateSelectInput(session, "treat", choices = cols)
    updateSelectInput(session, "outcome", choices = cols)
    updateSelectInput(session, "mediators", choices = cols)
    updateSelectInput(session, "intmed", choices = c("NULL" = " ", cols))
    updateSelectInput(session, "out_cov", choices = c("Aucun" = "", cols))
  })
  
  observeEvent(input$load_data, {
    req(input$data_file)
    ext <- tools::file_ext(input$data_file$name)
    
    output$debug_output <- renderPrint({
      paste("Fichier téléchargé :", input$data_file$name, "Extension :", ext)
    })
    
    loaded_data <- NULL  
    
    if (ext == "csv") {
      loaded_data <- tryCatch({
        read.csv(input$data_file$datapath)
      }, error = function(e) {
        output$debug_output <- renderPrint({ paste("Erreur lors du chargement du fichier CSV :", e$message) })
        NULL
      })
      
    } else if (ext == "RData") {
      tryCatch({
        load(input$data_file$datapath)  
        obj_names <- ls()[1]  
        if (length(obj_names) == 0) {
          stop("Aucun objet trouvé dans le fichier RData.")
        } else if (length(obj_names) == 1) {
          loaded_data <- get(obj_names[1])
        } else {
          output$debug_output <- renderPrint({
            paste("Le fichier RData contient plusieurs objets :", paste(obj_names, collapse = ", "), 
                  ". Veuillez spécifier lequel utiliser.")
          })
          return()
        }
      }, error = function(e) {
        output$debug_output <- renderPrint({ paste("Erreur lors du chargement du fichier RData :", e$message) })
      })
    }
    
    if (!is.null(loaded_data)) {
      CL <- colnames(loaded_data)
      if (is.array(loaded_data)) {
        loaded_data <- as.data.frame(loaded_data)  
        colnames(loaded_data) <- CL
      } else if (is.list(loaded_data)) {
        loaded_data <- as.data.frame(loaded_data)  
        colnames(loaded_data) <- CL
      }else {
        loaded_data <- as.data.frame(loaded_data)  
        colnames(loaded_data) <- CL
      }
      
      data(loaded_data)  
    }
    
    if (!is.null(data())) {
      output$debug_output <- renderPrint({
        paste("Dataset successfully instantiated. Total observations :", nrow(data()))
      })
    } else {
      output$debug_output <- renderPrint({cat("Status: No dataset currently residing in memory.")})
    }
  })
  
  selected_vars_list <- eventReactive(input$run,{
    req(input$data_file, input$outcome, input$treat)
    
    validate(
      validate(
        need(length(input$mediators) == 2, 
             "Please select exactly 2 mediators (Primary and Secondary) to proceed.")
      )
    )
    
    clean_input <- function(x) {
      if (is.null(x) || x == "") return(NULL)
      if (length(x) > 1) return(x)
      strsplit(x, ",")[[1]] %>% trimws()
    }
    
    list(
      outcome = input$outcome,
      exposure = input$treat,
      mediators = input$mediators,
      out_cov = clean_input(input$out_cov),
      intmed_cov = clean_input(input$intmed_cov),
      sed_cov = clean_input(input$sed_cov),
      interactions = input$interaction,
      all = unique(c(input$outcome, input$treat, input$mediators, 
                     clean_input(input$out_cov), clean_input(input$intmed_cov), 
                     clean_input(input$sed_cov)))
    )
  })
  get_residuals_df <- eventReactive(input$run,{
    req(data())
    
    outcome <- selected_vars_list()$outcome
    exposure <- selected_vars_list()$exposure
    mediators <- selected_vars_list()$mediators
    
    out_cov_vars <- selected_vars_list()$out_cov
    intmed_cov_vars <- selected_vars_list()$intmed_cov
    sed_cov_vars <- selected_vars_list()$sed_cov
    
    covariate <- unique(c(intmed_cov_vars, sed_cov_vars))
    
    inter <- as.logical(selected_vars_list()$interactions)
    selected_interactions <- input$interaction_types
    interaction_terms <- c()
    
    if (isTRUE(inter) && !is.null(selected_interactions)) {
      if ("treat_m1" %in% selected_interactions) {
        interaction_terms <- c(interaction_terms, paste(exposure, mediators[1], sep = ":"))
      }
      if ("treat_m2" %in% selected_interactions && length(mediators) > 1) {
        interaction_terms <- c(interaction_terms, paste(exposure, mediators[2], sep = ":"))
      }
      if ("m1_m2" %in% selected_interactions && length(mediators) > 1) {
        interaction_terms <- c(interaction_terms, paste(mediators[1], mediators[2], sep = ":"))
      }
      if ("treat_m1_m2" %in% selected_interactions && length(mediators) > 1) {
        triple <- paste(exposure, mediators[1], mediators[2], sep = ":")
        interaction_terms <- c(interaction_terms, triple)
      }
      lmfor <- paste(c(exposure, mediators, covariate, interaction_terms), collapse = " + ")
    }else{
      lmfor <- paste(c(exposure, mediators, covariate), collapse = " + ")
    }
    
    req(nzchar(outcome), nzchar(exposure), length(mediators) > 0)
    all_vars <- c(outcome, mediators)
    
    resid_mat <- sapply(all_vars, function(var) {
      if (var == outcome) {
        formula_str <- paste(var, "~", lmfor)
      } else {
        predictors <- c(exposure, covariate)
        formula_str <- paste(var, "~", paste(predictors, collapse = " + "))
      }
      model <- tryCatch({
        lm(as.formula(formula_str), data = data())
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(model)){
        return(rep(NA, nrow(data())))
      }
      residuals(model)
    })
    
    resid_df <- as.data.frame(resid_mat)
    colnames(resid_df) <- all_vars
    resid_df
  })

  output$data_table <- renderDataTable({
    req(data())  
    datatable(data(), options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px"))  
  })
  
  output$residual_table <- renderDataTable({
    df <- get_residuals_df()
    req(nrow(df) > 0)
    datatable(df, options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px"))
  })
  
  output$summary_stats_quant <- renderDataTable({
    vars <- selected_vars_list()$all
    req(length(vars) > 0)
    
    stats <- data() %>%
      select(any_of(vars)) %>%
      summarise(across(where(is.numeric), list(
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        min = ~min(., na.rm = TRUE),
        Q1 = ~quantile(., 0.25, na.rm = TRUE),
        median = ~median(., na.rm = TRUE),
        Q3 = ~quantile(., 0.75, na.rm = TRUE),
        max = ~max(., na.rm = TRUE),
        missing = ~sum(is.na(.))
      ))) 
    
    if (ncol(stats) == 0) {
      return(NULL) 
    }
    
   
    stats_tidy <- stats %>%
      pivot_longer(cols = everything(), names_to = "variable_statistic", values_to = "value") %>%
      separate(variable_statistic, into = c("variable", "statistic"), sep = "_(?=[^_]+$)", extra = "merge")%>%
      pivot_wider(names_from = statistic, values_from = value)%>%
      mutate(across(where(is.numeric), ~ round(.x, 5)))  
    
    
    stats_filtered <- stats_tidy
    if(length(vars)!= 0){
      stats_filtered <- stats_tidy %>%
        mutate(variable = factor(variable, levels = c(vars, unique(variable[!variable %in% vars]))))%>%
        arrange(variable)
    }
    datatable(stats_tidy, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$summary_stats_qual <- renderDataTable({
    
    qualitative_summary <- data() %>%
      mutate(across(where(is.numeric), ~ if (length(unique(.)) < 4){as.factor(.)}else{.})) %>%
      select(where(is.character), where(is.factor), where(is.ordered))
    
    if (ncol(qualitative_summary) == 0) {
      return(NULL) 
    }
    
    
    stats_filtered_quali <- qualitative_summary %>%
      drop_na() %>%  
      pivot_longer(cols = everything(), names_to = "Variable", values_to = "Valeur") %>%
      group_by(Variable, Valeur) %>%
      summarise(
        Frequence = n(),
        .groups = 'drop'
      )%>%
      group_by(Variable) %>% 
      mutate(Pourcentage = round((Frequence / sum(Frequence)) * 100,2)) 
    
    datatable(stats_filtered_quali, options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px"))
  })
  
  ###########################################################
  render_residual_plotly <- function(type, var_index) {
    renderPlotly({
      # 1. Préparation des données
      res_df <- get_residuals_df()
      vars <- colnames(res_df)
      req(length(vars) >= var_index)
      
      var_name <- vars[var_index]
      data_plot <- res_df[[var_name]]
      df <- data.frame(val = data_plot)
      
      # 2. Logique de génération des graphiques
      if (type == "box") {
        # --- BOXPLOT INTERACTIF ---
        p <- plot_ly(df, y = ~val, type = "box", 
                     name = var_name,
                     marker = list(color = '#337ab7'),
                     fillcolor = 'lightblue') %>%
          layout(yaxis = list(title = "Residual Value"))
        
      } else if (type == "hist") {
        # --- HISTOGRAMME + DENSITÉ ---
        # Calcul de la densité en amont pour la superposition
        dens <- density(data_plot)
        
        p <- plot_ly(df) %>%
          add_histogram(x = ~val, name = "Histogram", 
                        nbinsx = 30, histnorm = "probability density",
                        marker = list(color = 'lightgray', 
                                      line = list(color = 'white', width = 1))) %>%
          add_lines(x = dens$x, y = dens$y, name = "Density", 
                    line = list(color = 'blue', width = 2)) %>%
          layout(xaxis = list(title = "Residuals"),
                 yaxis = list(title = "Density"),
                 showlegend = FALSE)
        
      } else if (type == "qq") {
        # --- Q-Q PLOT INTERACTIF ---
        # Calcul des quantiles théoriques
        probs <- ppoints(length(data_plot))
        theo_quantiles <- qnorm(probs)
        sample_quantiles <- sort(data_plot)
        
        p <- plot_ly(x = theo_quantiles, y = sample_quantiles, 
                     type = 'scatter', mode = 'markers',
                     marker = list(color = '#337ab7', opacity = 0.6),
                     name = "Residuals") %>%
          # Ajout de la ligne de référence (Q-Q Line)
          add_lines(x = theo_quantiles, 
                    y = theo_quantiles * sd(sample_quantiles) + mean(sample_quantiles),
                    line = list(color = 'red', width = 2),
                    name = "Reference") %>%
          layout(xaxis = list(title = "Theoretical Quantiles"),
                 yaxis = list(title = "Sample Quantiles"),
                 showlegend = FALSE)
      }
      
      # 3. Personnalisation finale (Configuration et Layout)
      p %>% layout(
        margin = list(t = 40),
        hovermode = "closest",
        plot_bgcolor = "#f8f9fa"
      ) %>% 
        config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d"))
    })
  }
  
  # Boxplot
  output$boxplotoutcome <- render_residual_plotly ("box", 1)
  output$boxplotmedint <- render_residual_plotly ("box", 2)
  output$boxplotmedsed <- render_residual_plotly ("box", 3)
  
  # Histogram
  output$histogramoutcome <- render_residual_plotly ("hist", 1)
  output$histogrammedint <- render_residual_plotly ("hist", 2)
  output$histogrammedsed <- render_residual_plotly ("hist", 3)
  
  # Quantile
  output$qqplotoutcome <- render_residual_plotly ("qq", 1)
  output$qqplotmedint <- render_residual_plotly ("qq", 2)
  output$qqplotmedsed <- render_residual_plotly ("qq", 3)
  
  ##### Test de normalité
  render_norm_test <- function(var_index) {
    renderTable({
      res_df <- get_residuals_df()
      req(ncol(res_df) >= var_index)
      test_normalite_clean(res_df[[var_index]])
    }, striped = TRUE, hover = TRUE, bordered = TRUE)
  }
  #
  output$testoutcome <- render_norm_test(1)
  output$testmedint <- render_norm_test(2)
  output$testmedsed <- render_norm_test(3)
  output$testmulti <- renderTable({
    # Extraire les résidus du médiateur secondaire
    residuals_df <- get_residuals_df()
    
    mediator_vars <- input$mediators
    mediator_resid <- residuals_df[,mediator_vars]
    
    test_result <- MVN::mvn(mediator_resid, mvn_test = "hz", univariate_test = "SW",
                       show_new_data = FALSE) 
    print(test_result$multivariate_normality)
  }, striped = TRUE, hover = TRUE)
  
  
  output$Outcomebox <- renderText({selected_vars_list()$outcome})
  output$Outcomehist <- renderText({selected_vars_list()$outcome})
  output$Outcomeqq <- renderText({selected_vars_list()$outcome})
  output$Outcomenorm <- renderText({selected_vars_list()$outcome})
  
  output$Primarymediatorbox <- renderText({selected_vars_list()$mediators[1]})
  output$Primarymediatorhist <- renderText({selected_vars_list()$mediators[1]})
  output$Primarymediatorqq <- renderText({selected_vars_list()$mediators[1]})
  output$Primarymediatornorm <- renderText({selected_vars_list()$mediators[1]})
  
  output$Secondarymediatorbox <- renderText({selected_vars_list()$mediators[2]})
  output$Secondarymediatorhist <- renderText({selected_vars_list()$mediators[2]})
  output$Secondarymediatorqq <- renderText({selected_vars_list()$mediators[2]})
  output$Secondarymediatornorm <- renderText({selected_vars_list()$mediators[2]})
  
  
  output$Multivariatenorm <- renderText({selected_vars_list()$mediators})
  
  ####### Corrélation constante
  

  resultsCC <- eventReactive(input$run, {
    req(data(), selected_vars_list())
    
    vars        <- selected_vars_list()
    treat_var   <- vars$exposure
    outcome_var <- vars$outcome
    mediator_vars <- vars$mediators 
    out_cov_vars  <- vars$out_cov
    intmed_cov_vars <- vars$intmed_cov
    sed_cov_vars    <- vars$sed_cov
    
    Br <- input$B
    inter <- as.logical(vars$interactions)
    selected_interactions <- input$interaction_types
    cor_cste <- 1
    rh <- 0.5
    
    covariate <- unique(c(intmed_cov_vars, sed_cov_vars))
    interaction_terms <- c()
    names_vec <- NULL
    
    if (isTRUE(inter) && !is.null(selected_interactions)) {
      if ("treat_m1" %in% selected_interactions) 
        interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[1], sep = ":"))
      if ("treat_m2" %in% selected_interactions) 
        interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[2], sep = ":"))
      if ("m1_m2" %in% selected_interactions) 
        interaction_terms <- c(interaction_terms, paste(mediator_vars[1], mediator_vars[2], sep = ":"))
      if ("treat_m1_m2" %in% selected_interactions) 
        interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[1], mediator_vars[2], sep = ":"))
      
      lmfor <- paste(c(treat_var, mediator_vars, covariate, interaction_terms), collapse = " + ")
      names_vec <- c(paste(treat_var, mediator_vars, sep=":"), 
                     paste(mediator_vars, collapse = ":"), 
                     paste(treat_var, paste(mediator_vars, collapse = ":"), sep=":"))
    } else {
      lmfor <- paste(c(treat_var, mediator_vars, covariate), collapse = " + ")
    }
    
    run_analysis_core <- function(current_data) {
      tableCC <- initialParams(
        treat = treat_var, mediators = mediator_vars, intmed = mediator_vars[1],
        outcome = outcome_var, intmed_cov = intmed_cov_vars, sed_cov = sed_cov_vars, 
        out_cov = out_cov_vars, inter_treat_cov = FALSE, cor_cste = cor_cste, 
        data = current_data, formula_one3 = lmfor, names_vec = names_vec
      )
      
      q1 <- 2 + length(intmed_cov_vars)
      q2 <- 2 + length(sed_cov_vars)
      
      coefs_med <- data.frame(
        matrix(c(tableCC$med[1:q1], tableCC$med[(q1+1):(q1+q2)]), 2, q1, byrow = TRUE), 
        mediator_vars
      )
      colnames(coefs_med) <- c("inter", treat_var, intmed_cov_vars, "name")
      
      bet_names <- if(isTRUE(inter)) {
        c("inter", treat_var, mediator_vars, out_cov_vars, names_vec, "sd")
      } else {
        c("inter", treat_var, mediator_vars, out_cov_vars, "sd")
      }
      Bet <- data.frame(Beta = tableCC$outc, name = bet_names)
      
      cor_coefs <- as.vector(tableCC$med[-c(1:(q1+q2))])
      var_covar <- varcovarEstimes(cor_coefs, cor_cste = cor_cste)
      r011 <- as.vector(var_covar$r01) 
      
      sol <- effectdirectindirct(
        alpha = coefs_med, beta = Bet, treat = treat_var, mediators = mediator_vars,
        intmed_cov = intmed_cov_vars, sed_cov = sed_cov_vars, inter = inter, ro = r011,
        corC = cor_coefs, names_vec = names_vec, cor_cste = cor_cste, data = current_data
      )
      return(c(DE = sol$DE1sk, IE = sol$IE1sk))
    }
    
    # 5. Calcul sur les données réelles
    SolVrai <- run_analysis_core(data())
    
    #6. Bootstrap (Vectorisé)
    set.seed(145)
    boot_list <- lapply(1:Br, function(i) {
      sample_data <- data()[sample(nrow(data()), replace = TRUE), ]
      res_boot <- run_analysis_core(sample_data)$res
      return(as.numeric(res_boot)) 
    })
    
    SolB <- do.call(rbind, boot_list)
    
    # 7. Résultats finaux
    SolBi <- apply(SolB, 2, var)
    SolBcant <- apply(SolB, 2, quantile, probs = c(0.025, 0.975))
    
    boot_means <- colMeans(SolB, na.rm = TRUE)
    biais_bootstrap <- boot_means - SolVrai
    
    tibble::tibble(
      Effect = c("Direct Effect", "Indirect Effect"),
      Estimate = round(SolVrai, 2),
      StdError = round(sqrt(SolBi), 2),
      BiasBoost = round(biais_bootstrap,2),
      `95% CI` = paste0("[", round(SolBcant[1, ], 2), ", ", round(SolBcant[2, ], 2), "]")
    )
  })
  
  output$table_results_cc <- renderDataTable({
    df <- resultsCC()

    datatable(df, 
              caption = "Direct and Indirect Effects via Causal Categorical approach",
              options = list(dom = 't',
                             columnDefs = list(list(className = 'dt-center', targets = "_all"))),
              selection = 'none') %>%
      formatStyle('Effect', fontWeight = 'bold') 
  })
  
  # # 1. Valeur Effet Direct
  # output$cc_direct_val <- renderText({
  #   df <- resultsCC()
  #   # On récupère l'estimate de la première ligne (ajustez l'index selon votre df)
  #   round(df$Estimate[1], 3)
  # })
  # 
  # # 2. CI Effet Direct
  # output$cc_direct_ci <- renderText({
  #   df <- resultsCC()
  #   paste("CI 95%:", df$CI[1])
  # })
  # 
  # # 3. Valeur Effet Indirect
  # output$cc_indirect_val <- renderText({
  #   df <- resultsCC()
  #   # On récupère l'estimate de l'effet indirect (ex: ligne 7)
  #   round(df$Estimate[7], 3)
  # })
  # 
  # # 4. CI Effet Indirect
  # output$cc_indirect_ci <- renderText({
  #   df <- resultsCC()
  #   paste("CI 95%:", df$CI[7])
  # })
  
  
  resultsCNC <- eventReactive(input$run, {
    req(data(), selected_vars_list())
    
    vars        <- selected_vars_list()
    treat_var   <- vars$exposure
    outcome_var <- vars$outcome
    mediator_vars <- vars$mediators
    out_cov_vars  <- vars$out_cov
    intmed_cov_vars <- vars$intmed_cov
    sed_cov_vars    <- vars$sed_cov
    
    Br <- input$B
    inter <- as.logical(vars$interactions)
    selected_interactions <- input$interaction_types
    cor_cste <- 3
    rh <- 0.5
    
    covariate <- unique(c(intmed_cov_vars, sed_cov_vars))
    interaction_terms <- c()
    names_vec <- NULL
    
    if (isTRUE(inter)) {
      if (!is.null(selected_interactions)) {
        if ("treat_m1" %in% selected_interactions) 
          interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[1], sep = ":"))
        if ("treat_m2" %in% selected_interactions && length(mediator_vars) > 1) 
          interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[2], sep = ":"))
        if ("m1_m2" %in% selected_interactions && length(mediator_vars) > 1) 
          interaction_terms <- c(interaction_terms, paste(mediator_vars[1], mediator_vars[2], sep = ":"))
        if ("treat_m1_m2" %in% selected_interactions && length(mediator_vars) > 1) 
          interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[1], mediator_vars[2], sep = ":"))
      }
      
      names_vec <- c(
        paste(treat_var, mediator_vars, sep = ":"), 
        paste(mediator_vars, collapse = ":"), 
        paste(treat_var, paste(mediator_vars, collapse = ":"), sep = ":")
      )
      lmfor <- paste(c(treat_var, mediator_vars, covariate, interaction_terms), collapse = " + ")
    } else {
      lmfor <- paste(c(treat_var, mediator_vars, covariate), collapse = " + ")
    }
    
    run_analysis_core <- function(current_data) {
      tableCNC <- initialParams(
        treat = treat_var, mediators = mediator_vars, intmed = mediator_vars[1],
        outcome = outcome_var, intmed_cov = intmed_cov_vars, sed_cov = sed_cov_vars, 
        out_cov = out_cov_vars, inter_treat_cov = FALSE, cor_cste = cor_cste, 
        data = current_data, formula_one3 = lmfor, names_vec = names_vec
      )
      
      q1 <- 2 + length(intmed_cov_vars)
      q2 <- 2 + length(sed_cov_vars)
      
      coefs_med <- data.frame(
        matrix(c(tableCNC$med[1:q1], tableCNC$med[(q1+1):(q1+q2)]), 2, q1, byrow = TRUE), 
        mediator_vars
      )
      colnames(coefs_med) <- c("inter", treat_var, intmed_cov_vars, "name")
      
      bet_names <- if(isTRUE(inter)) {
        c("inter", treat_var, mediator_vars, out_cov_vars, names_vec, "sd")
      } else {
        c("inter", treat_var, mediator_vars, out_cov_vars, "sd")
      }
      Bet <- data.frame(Beta = tableCNC$outc, name = bet_names)
      
      cor_coefs <- as.vector(tableCNC$med[-c(1:(q1+q2))])
      var_covar <- varcovarEstimes(cor_coefs, cor_cste = cor_cste)
      
      r011 <- c(as.vector(var_covar$r01), (cor_coefs[5] + cor_coefs[6]) / 2, rh)
      
      sol <- effectdirectindirct(
        alpha = coefs_med, beta = Bet, treat = treat_var, mediators = mediator_vars,
        intmed_cov = intmed_cov_vars, sed_cov = sed_cov_vars, inter = inter, ro = r011,
        corC = cor_coefs, names_vec = names_vec, cor_cste = cor_cste, data = current_data
      )

      return(list(res = c(sol$DE1sk, sol$IE1sk), rho = r011))
    }
    
    vrai_obj <- run_analysis_core(data())
    SolVrai  <- vrai_obj$res
    RhosVrai <- vrai_obj$rho
    
    set.seed(145)
    boot_list <- lapply(1:Br, function(i) {
      sample_data <- data()[sample(nrow(data()), replace = TRUE), ]
      res_boot <- run_analysis_core(sample_data)$res
      return(as.numeric(res_boot)) 
    })

    SolB <- do.call(rbind, boot_list) 
    
    SolBi <- apply(SolB, 2, var, na.rm = TRUE)
    SolBcant <- apply(SolB, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    
    boot_means <- colMeans(SolB, na.rm = TRUE)
    biais_bootstrap <- boot_means - SolVrai
    
    dtf <- tibble::tibble(
      Effect = rep(c("Direct Effect", "Indirect Effect"), each = 6),
      Method = rep(c("CNCr1", "CNCr2", "CNCr3", "CNCr4", "CNCm", "CNC"), 2),
      rho    = rep(round(RhosVrai, 2), 2),
      Estimate = round(SolVrai, 2),
    StdError = round(sqrt(SolBi), 2),
    BiasBoost = round(biais_bootstrap,2),
    `95% CI` = paste0("[", round(SolBcant[1, ], 2), ", ", round(SolBcant[2, ], 2), "]")
    )
    
    # Affiche les 6 premières lignes dans la console R pour inspection
    # Réorganisation des lignes selon votre logique slice(c(5,1,2,3,4,11,7,8,9,10))
    # Note : Vérifiez bien l'ordre des lignes si vos fonctions retournent plus ou moins d'effets
    dtf %>% slice(c(5, 1, 2, 3, 4, 11, 7, 8, 9, 10))
    #dtf
  })
  
  output$current_B <- renderText({
    req(input$B) 
    input$B
  })
  
  
  render_cnc_boxes <- function(df, effect_name, color_theme) {
    df_sub <- df[grepl(effect_name, df$Effect), ]
    
    layout_column_wrap(
      width = "250px",
      fixed_width = TRUE,
      gap = "1.5rem", 
      lapply(1:nrow(df_sub), function(i) {
        is_mean <- df_sub$Method[i] == "CNCm"
        value_box(
          title = df_sub$Method[i],
          value = span(format(round(as.numeric(df_sub$Estimate[i]), 3), nsmall = 2), 
                       style = "font-size: 4.01rem;"),
          theme = if(is_mean) color_theme else "light",
          
          tags$div(
            #style = "font-size: 0.85rem; line-height: 1.4;",
            tags$div(
              class = "d-flex justify-content-between",
              style = "font-size: 1.55rem;",
              tags$span(tags$b("Bias: "), format(df_sub$BiasBoost[i], nsmall = 3),", "),
              tags$span(tags$b("SE: "), format(df_sub$StdError[i], nsmall = 2)),
            ),
            tags$div(
              style = "font-size: 1.31rem;",
              tags$b("95% CI: "), df_sub$`95% CI`[i]
            )
          )
        )
      })
    )
  }
  
  output$cnc_direct_grid <- renderUI({
    req(resultsCNC())
    render_cnc_boxes(resultsCNC(), "Direct", "primary")
  })
  
  output$cnc_indirect_grid <- renderUI({
    req(resultsCNC())
    render_cnc_boxes(resultsCNC(), "Indirect", "info") 
  })
  
  render_rho_params <- function(df) {
    df_sub <- df[df$Effect == "Indirect Effect", ]
    
    layout_column_wrap(
      width = "250px", 
      fixed_width = TRUE,
      gap = "1.5rem",
      
      lapply(1:nrow(df_sub), function(i) {
        is_mean <- df_sub$Method[i] == "CNCm"
        
        value_box(
          title = df_sub$Method[i],
          value = span(format(round(as.numeric(df_sub$rho[i]), 2), nsmall = 2), 
                       style = "font-size: 3.95rem;"),
          theme = if(is_mean) "primary" else "secondary",
          br(),
          br(),
          br()
        )
      })
    )
  }
 
  output$rho <- renderUI({
    req(resultsCNC())
    render_rho_params(resultsCNC()) 
  })
  
  
  output$cc_rho <- renderUI({
    req(resultsCC())
    render_rho_params(resultsCC()) 
  })
  
  output$cc_direct_grid <- renderUI({
    req(resultsCC())
    render_cc_boxes(resultsCC(), "Direct", "primary")
  })
  
  output$cc_indirect_grid <- renderUI({
    req(resultsCC())
    render_cc_boxes(resultsCC(), "Indirect", "info") 
  })
  
  
}


# Launch Dashboard
shinyApp(ui, server)  # Initialize the shinydashboard
