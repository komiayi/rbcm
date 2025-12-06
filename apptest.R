library(shiny)
library(DT)

ui <- fluidPage(
  titlePanel("Test résidus"),
  sidebarLayout(
    sidebarPanel(
      selectInput("outcome", "Outcome", choices = names(mtcars), selected = "mpg"),
      selectInput("mediators", "Médiateurs", choices = names(mtcars), selected = "hp", multiple = TRUE),
      selectInput("exposure", "Exposition", choices = names(mtcars), selected = "wt")
    ),
    mainPanel(
      dataTableOutput("residual_table")
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({ mtcars })
  
  get_residuals_df <- reactive({
    req(data())
    outcome <- input$outcome
    mediators <- input$mediators
    exposure <- input$exposure
    
    req(nzchar(outcome), nzchar(exposure), length(mediators) > 0)
    
    all_vars <- c(outcome, mediators)
    req(all(all_vars %in% names(data())), exposure %in% names(data()))
    
    resid_mat <- sapply(all_vars, function(var) {
      predictors <- if (var == outcome) c(exposure, mediators) else c(exposure)
      formula_str <- paste(var, "~", paste(predictors, collapse = " + "))
      model <- tryCatch(lm(as.formula(formula_str), data = data()), error = function(e) return(NULL))
      if (is.null(model)) return(rep(NA, nrow(data())))
      residuals(model)
    })
    
    resid_df <- as.data.frame(resid_mat)
    colnames(resid_df) <- all_vars
    resid_df
  })
  
  output$residual_table <- renderDataTable({
    df <- get_residuals_df()
    req(nrow(df) > 0)
    datatable(df, options = list(pageLength = 10))
  })
}

shinyApp(ui, server)
