library(shiny)
library(shinydashboard)
library(plotly)
library(leaflet)
library(reshape2)
library(leaflet.extras)
library(shinyWidgets)
library(Polychrome)
library(dplyr)
library(tseries)
library(DT)
source("Estimations.R")

loop.detect <- as.vector(sapply(routes, function(x) strsplit(x, "_")[[1]][1]))
app_methods <- c("MTE", "LSA")
criteria <- c("AIC", "BIC")
lm.reg  <- lm(paste0("L101_volume", f), data = data)
X <- as.matrix(cbind("(Intercept)"=1,lm.reg$model[,-1]))

extract_coef_helper <- function(extract = c("Mean", "2.5q", "97.5q","Median"),coef_ = c("MEB", "Actual", "ADL1LASSO")){
  func <- function(x, extr) {
    switch(
      extr,
      "Median" = quantile(x,0.5),
      "2.5q" = quantile(x,0.025),
      "97.5q" = quantile(x,0.975),
      mean(x)
    )
  }
  RES <- lapply(routes,function(loc) {
    res <- data.frame(
      apply(
        cbind(
          "L1MTE_AIC" = c(T,T,"AIC"),
          "L1MTE_BIC" = c(T,T,"BIC"),
          "L1LSA_AIC" = c(T,F,"AIC"),
          "L1LSA_BIC" = c(T,F,"BIC"),
          "L2MTE_AIC" = c(F,T,"AIC"),
          "L2MTE_BIC" = c(F,T,"BIC"),
          "L2LSA_AIC" = c(F,F,"AIC"),
          "L2LSA_BIC" = c(F,F,"BIC")), 2,
        function(x) {
          if (coef_=="MEB") {apply(plotData_meb(loc,c(as.logical(x[1]),as.logical(x[2])),x[3],show_ ="Coefficients")[,-21],1,function(z) func(z,extract))}
          else if (coef_=="Actual") {plotData(loc,c(as.logical(x[1]),as.logical(x[2])),x[3],F,show_ ="Coefficients")[,gsub("_volume","",loc)]}
          else {plotData(loc,c(as.logical(x[1]),as.logical(x[2])),x[3],T,show_ ="Coefficients")[,1]}
          }),
      "Predictor" = coeff_names)
    res <- res %>%
      arrange(across(1:(length(res)-1), ~ desc(.)))
    res[rowSums(abs(res[,-length(res)]))>0,]  
  })
  names(RES) <- routes
  RES
}
extract_coef <-function(loc){
  Mmean     <- extract_coef_helper("Mean","MEB")[[loc]]
  Mmed      <- extract_coef_helper("Median","MEB")[[loc]]
  M97.5     <- extract_coef_helper("2.5q","MEB")[[loc]]
  M2.5      <- extract_coef_helper("97.5q","MEB")[[loc]]
  act       <- extract_coef_helper("","Actual")[[loc]]
  adl1lasso <- extract_coef_helper("","ADL1LASSO")[[loc]]
  adl1lasso <- adl1lasso[adl1lasso$Predictor %in% act$Predictor,]

  columnNames <- colnames(Mmean)
  numCols <-  length(columnNames)
  dfs <- list(Mmean, Mmed, M97.5, M2.5, act, adl1lasso)
  joined_dt <- reduce(dfs, dplyr::full_join, by = "Predictor")
  joined_dt[is.na(joined_dt)] <- 0
  preds <- joined_dt$Predictor
  joined_dt <- joined_dt[,!(names(joined_dt) %in% "Predictor")]
  numCols <-  length(columnNames)-1
  
  Mmean       <- cbind(joined_dt[,1:numCols], preds)
  Mmed        <- cbind(joined_dt[,(numCols+1):(2*numCols)], preds)
  M97.5       <- cbind(joined_dt[,(2*numCols+1):(3*numCols)], preds)
  M2.5        <- cbind(joined_dt[,(3*numCols+1):(4*numCols)], preds)
  act         <- cbind(joined_dt[,(4*numCols+1):(5*numCols)], preds)
  adl1lasso   <- cbind(joined_dt[,(5*numCols+1):(6*numCols)], preds)
  
  colnames(Mmean) =  colnames(Mmed) =  colnames(M97.5) =  colnames(M2.5) =  colnames(act) =  colnames(adl1lasso) = columnNames
  list("Mmean" = Mmean, "Mmed" = Mmed, "M97.5" = M97.5, "M2.5" = M2.5, "act" = act, "adl1lasso" = adl1lasso)
}

df_locations <- data.frame(
  name = c("L101", "L102", "L103", "L104", "L106", "L107", "L108"),
  lat = c(37.9916, 37.98942, 37.988885, 37.98721, 37.9894, 37.9890, 37.98719),
  lng = c(23.73255, 23.74586, 23.74892, 23.75793, 23.74548, 23.74767, 23.75743)
)


ui <- fluidPage(
  title = "A Comparative Analysis of Robust Penalized Estimators for Periodic Time Series",
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "coding.png")
  ),
  div(id = "dark-mode-switch-container",
      style="padding: 10px; background-color: #f8f9fa; border-bottom: 1px solid #dee2e6;", 
      shinyWidgets::materialSwitch(
        inputId = "dark_mode",
        label = "Dark Theme",
        status = "primary",
        right = TRUE,
        inline = TRUE,
        value = FALSE
      )
  ),
  dashboardPage(
    skin = "blue",
    dashboardHeader(title = textOutput("dashboard_title")), 
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        hr(style="border-top: 1px solid #ccc;"), 
        menuItem("Settings", tabName = "settings_tab", icon = icon("sliders"), selected = TRUE),
        selectInput("method", "Estimation Method:", choices = app_methods, selected = app_methods[1]),
        checkboxInput("robust", "Use Robust Method", value = FALSE),
        selectInput("criterion", "Information Criterion:", choices = criteria, selected = criteria[1]),
        selectInput("route", "Loop Detector:", choices = loop.detect, selected = loop.detect[1]),
        selectInput("plot_type", "Select Plot Type:",
                    choices = c("Coefficients", "Fitted Values", "Residuals", "Predict Values", "Prediction Errors"),
                    selected = "Fitted Values"),
        menuItem("Analysis Plots", tabName = "plots", icon = icon("chart-line")),
        menuItem("Coeficients: MEB against Initial data ", tabName = "MEB_vs_ACTUAL", icon = icon("balance-scale")),
        menuItem("Coeficients: Location-Wise", tabName = "locationWise", icon = icon("compass")),
        menuItem("Variable Importance", tabName = "vips", icon = icon("ranking-star")),
        menuItem("Normality and Stationarity Tests", tabName = "validation_tab", icon = icon("vial")),
        menuItem("Detector Map", tabName = "map", icon = icon("map-location-dot"))#,
        # menuItem("Detailed Report (PDF)", tabName = "pdf", icon = icon("file-alt", class = "fas"))
      )
    ),
    dashboardBody(
      tags$head(
        tags$style(HTML("
              body { margin: 0; font-family: sans-serif; }
              .plotly { font-family: sans-serif !important; } /* Ensure plotly uses consistent font */
              
              /* Center h2 and h3 in the main content */
              .content-wrapper h2,
              .content-wrapper h3 {
                text-align: center;
                margin-left: 0; /* default margin */
                transition: margin-left 0.3s ease;
              }
              
              /* Shift h2 and h3 1cm to the right when sidebar is visible (not collapsed) */
              body:not(.sidebar-collapse) .content-wrapper h2,
              body:not(.sidebar-collapse) .content-wrapper h3 {
                margin-left: 1cm;
              }  
        ")),
        tags$style(HTML("
              :root { --sidebar-width: 300px; }
              .main-sidebar {
                 width: var(--sidebar-width) !important;
              }
             
              .main-header .logo {
                 width: var(--sidebar-width) !important;
              }
             
             .main-header .navbar {
                 margin-left: var(--sidebar-width) !important;
             }
             .content-wrapper > .content { /* Target the inner content */
                 max-width: none !important;
                 width: 95% !important; /* Use 95% of available space */
                 margin-left: auto !important;
                 margin-right: auto !important;
                 padding-left: 30px; /* Add some padding */
                 padding-right: 30px;
             }
             .box {
                 width: 100% !important;
             }
             
             body.sidebar-collapse .main-sidebar {
                 transform: translateX(-100%) !important;
                 overflow: hidden !important;
             }
             body.sidebar-collapse .content-wrapper,
             body.sidebar-collapse .right-side,
             body.sidebar-collapse .main-footer {
                 margin-left: 0 !important;
             }
             body.sidebar-collapse .main-header .logo {
                 width: 0 !important;
                 display: none !important; /* Ensure it's fully hidden */
             }
             body.sidebar-collapse .main-header .navbar {
                 margin-left: 0 !important;
             }
             
             
             body { margin: 0; font-family: sans-serif; }
             .plotly { font-family: sans-serif !important; }
              #dark-mode-switch-container { transition: background-color 0.3s ease; } /* Smooth transition */

              body.dark-mode #dark-mode-switch-container {
                   background-color: #1a2226 !important; /* Match dark header */
                   border-bottom: 1px solid #374850 !important;
              }
              body.dark-mode #dark-mode-switch-container label { /* Ensure label text is visible */
                   color: #ecf0f1 !important;
              }
              body.dark-mode {
                background-color: #222d32 !important; /* Dark background */
                color: #bdc3c7 !important; /* Light grey text */
              }
              body.dark-mode .content-wrapper,
              body.dark-mode .main-sidebar, /* Sidebar already styled below, but good fallback */
              body.dark-mode .box,
              body.dark-mode .info-box,
              body.dark-mode .nav-tabs-custom > .nav-tabs > li.active {
                background-color: #2c3b41 !important; /* Slightly lighter dark for components */
                color: #ecf0f1 !important; /* Lighter text */
              }
              body.dark-mode .box-header { border-bottom: 1px solid #374850 !important; }
              body.dark-mode h1, body.dark-mode h2, body.dark-mode h3,
              body.dark-mode h4, body.dark-mode h5, body.dark-mode h6,
              body.dark-mode .box-title, body.dark-mode label {
                color: #ecf0f1 !important; /* Ensure headers and labels are light */
              }
              body.dark-mode .form-control, body.dark-mode .selectize-input {
                  background-color: #374850 !important;
                  color: #ecf0f1 !important;
                  border: 1px solid #5a6a72 !important;
              }
               body.dark-mode .selectize-dropdown-content .active {
                  background-color: #5F9EA0 !important; /* Match hover color */
              }
              body.dark-mode .main-header .logo,
              body.dark-mode .main-header .navbar {
                background-color: #1a2226 !important; /* Very dark header */
              }
               body.dark-mode .leaflet-tile, body.dark-mode .leaflet-control-layers-base label span{
                  filter: brightness(0.7) invert(1) contrast(1.1) hue-rotate(180deg); /* Dark mode tiles */
               }
               body.dark-mode .leaflet-control-attribution { background-color: rgba(0,0,0,0.7) !important; color: #ccc !important;}
               body.dark-mode .leaflet-popup-content-wrapper { background-color: #2c3b41 !important; color: #ecf0f1 !important;}
               body.dark-mode .leaflet-popup-tip { background-color: #2c3b41 !important;}
               body.dark-mode .modebar { background-color: rgba(44, 59, 65, 0.9) !important; }
               body.dark-mode .modebar-btn path { fill: #ecf0f1 !important; }
               body.dark-mode .js-plotly-plot .plotly .cursor-crosshair { stroke: white !important; }
               body.dark-mode .js-plotly-plot .plotly .xaxislayer-above .xtick text,
               body.dark-mode .js-plotly-plot .plotly .yaxislayer-above .ytick text { fill: #ecf0f1 !important; }
               body.dark-mode .js-plotly-plot .plotly .legendtext { fill: #ecf0f1 !important; }
               body.dark-mode .js-plotly-plot .plotly .annotation-text { fill: #ecf0f1 !important; }
               body.dark-mode .radio label, body.dark-mode .checkbox label { color: #ecf0f1 !important; } /* Fix radio button labels */

               /* Explicit light mode box override if needed */
               body.dark-mode .keep-light {
                  background-color: white !important;
                  color: black !important;
                }
                body.dark-mode .keep-light .box-header, body.dark-mode .keep-light h1, body.dark-mode .keep-light h2, body.dark-mode .keep-light h3, body.dark-mode .keep-light h4, body.dark-mode .keep-light label {
                  color: black !important;
                }
          ")),
        tags$style(HTML("
              .skin-blue .main-sidebar {
                background-color: steelblue !important;
              }
              .skin-blue .main-sidebar .sidebar a,
              .skin-blue .main-sidebar .sidebar .header,
              .skin-blue .main-sidebar .sidebar .radio label, /* Style radio label */
              .skin-blue .main-sidebar .sidebar .checkbox label,
              .skin-blue .main-sidebar .sidebar .form-group label {
                color: white !important;
              }
              .skin-blue .main-sidebar .sidebar .sidebar-menu .active a {
                background-color: #3c8dbc !important;
              }
              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover {
                background-color: #5F9EA0 !important;
              }
              .skin-blue .main-sidebar .sidebar .form-control,
              .skin-blue .main-sidebar .sidebar .selectize-input {
                  background-color: rgba(255,255,255, 0.2); /* Light background for inputs on sidebar */
                  color: white;
                  border: 1px solid #a9d1ec;
              }
              .skin-blue .main-sidebar .sidebar .selectize-dropdown-content .active {
                  background-color: #5F9EA0 !important; /* Match hover color */
              }
          ")),
        tags$script(HTML("
              Shiny.addCustomMessageHandler('toggle-dark-theme', function(dark) {
                if (dark) {
                  document.body.classList.add('dark-mode');
                } else {
                  document.body.classList.remove('dark-mode');
                }
              });
          "))
      ),
      
      box(
        title = "A Comparative Analysis of Robust Penalized Estimators for Periodic Time Series",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        collapsed = FALSE,
        div(style = "font-size: 14px; line-height: 1.6;",
            tags$p("This interactive dashboard presents a comprehensive comparison between two advanced estimation techniques,",
                   tags$b("Maximum Tangent Likelihood Estimation (MTE)"), " and ",
                   tags$b("Least Squares Approximation (LSA)"), ", within the framework of the Unified Lasso methodology."),
            tags$p("The application is designed to facilitate robust model assessment, variable selection, and predictive performance evaluation across multiple traffic sensor locations."),
            tags$hr(),
            tags$p("Key features of this tool include:"),
            tags$ul(
              tags$li("Selection of estimation method: ", tags$b("MTE"), " or ", tags$b("LSA")),
              tags$li("Optional activation of ", tags$b("robustification procedures")),
              tags$li("Model evaluation using information criteria: ", tags$b("AIC"), " or ", tags$b("BIC")),
              tags$li("Focused analysis ", tags$b("across multiple loop detector locations")," (L101-L108)"),
              tags$li("Visualization of estimated ", tags$b("coefficients"),", ",tags$b("residual diagnostics"),", ", tags$b("fitted values"),", and ", tags$b("forecast errors")),
              tags$li("Computation and comparison of ", tags$b("Variable Importance Proportions (VIP)")),
              tags$li(tags$b("Inspection of residual")," normality and stationarity"),
              tags$li(tags$b("Interactive mapping")," of detector locations for spatial context (on Alexandras\' Avenue)")
            ),
            tags$hr(),
            tags$p("A ", tags$b("dark mode toggle"), " is also available to improve readability in low-light settings."),
            tags$p("To begin, use the sidebar to configure your analytical preferences. The dashboard content will adapt accordingly.")
        )
      )
      
      ,
      tabItems(
        tabItem("settings_tab"),
        tabItem("plots",
                fluidRow(
                  box(title = textOutput("plot_title"), status = "primary", solidHeader = TRUE, width = 12,
                      plotlyOutput("main_plot", height = "600px")
                  )
                )
        ),
        tabItem("MEB_vs_ACTUAL",
                fluidRow(
                  box(title = textOutput("MEB_vs_ACTUAL_title"), status = "primary", solidHeader = TRUE, width = 12,
                      plotlyOutput("comparison_plot", height = "600px") 
                  )
                )
        ),
        tabItem("locationWise",
                fluidRow(
                  box(
                    title = textOutput("locationWise_title"), 
                    status = "info", 
                    solidHeader = TRUE, width = 12,
                    radioButtons("data_source", "Select which coefficients to see:",
                                 choices = c("Average from MEB data", "From Initial Data", "Post-hoc Adaptive L1 Lasso on Initial Data"),
                                 selected = "MEB",
                                 inline = TRUE),
                    hr(), 
                    plotlyOutput("locationWise_plot", height = "600px")
                  )
                )
        ),
        tabItem("vips",
                fluidRow(
                  box(title = textOutput("vip_title"), status = "primary", solidHeader = TRUE, width = 12,
                      plotlyOutput("vip_plot", height = "600px"),
                      plotlyOutput("heatmap_plot", height = "600px")
                  )
                )
        ),
        tabItem(tabName = "validation_tab",
                fluidRow(
                  box(title = "Violin Plots for all Metrics", status = "primary", solidHeader = TRUE, width = 12,
                      plotlyOutput("metrics_violin_plot"),
                  ),
                  box(title = "Violin data table for all Metrics", status = "primary", solidHeader = TRUE, width = 12,
                      DTOutput("artifacts_table")
                  )
                ),
                fluidRow(
                  box(title = "ACF/PACF for Residuals and Prediction Errors", status = "primary", solidHeader = TRUE, width = 12,
                      plotlyOutput("acf_pacf_plot", height = "600px")
                  )
                ),
                fluidRow(
                  box(title = "Normality and Stationarity Tests", status = "primary", solidHeader = TRUE, width = 12,
                      verbatimTextOutput("tests_summary")
                  )
                )
        ),
        tabItem("map",
                fluidRow(
                  box(title = "Loop Detector Locations", status = "primary", solidHeader = TRUE, width = 12,
                      leafletOutput("map_output", height = "600px")
                  )
                )
        ),
        tabItem(tabName = "pdf",
                fluidRow(
                  box(title = "Thesis PDF: Detailed report of the afore-shown analysis", status = "primary", solidHeader = TRUE, width = 12,
                      tags$iframe(
                        style = "height:800px; width:100%; border:none;",
                        src = "A_Comparative_Analysis_of_Robust_Penalized_Estimators_for_Periodic_Time_Series.pdf"
                      )
                  )
                )
        )
      ) 
    ) 
  ) 
) 

server <- function(input, output, session) {
  observe({
    session$sendCustomMessage("toggle-dark-theme", input$dark_mode)
  })
  
  reactive_plot_data <- reactive({
    req(input$method, input$criterion, input$route, input$plot_type)
    
    location_full_name <- paste0(input$route,"_volume") 
    if(length(location_full_name) == 0) {
      showNotification("Selected route mapping failed.", type="error")
      return(NULL)
    }
    
    is_robust_param <- input$robust
    is_mte_param <- (input$method == "MTE")
    method_param <- c(is_robust_param, as.integer(is_mte_param))
    
    current_plot_type <- input$plot_type
    
    tryCatch({
      result <- plotData_meb(
        location = location_full_name,
        method = method_param,
        crit = input$criterion,
        show_ = current_plot_type
      )
      validate(
        need(result, "Data loading function returned NULL."),
        need(is.data.frame(result) && nrow(result) > 0, "No data returned for the selected parameters.")
      )
      result
    }, error = function(e) {
      showNotification(paste("Data Load Error:", e$message), type = "error", duration = 10)
      NULL 
    })
  })
  
  reactive_plotData_loc  <- reactive({
    req(input$data_source, input$method, input$criterion) 
    is_robust_param <- input$robust
    is_mte_param <- (input$method == "MTE")
    method_param <- c(is_robust_param, as.integer(is_mte_param))
    selected_criterion <- input$criterion
    num_locations <- length(routes)
    predictor_names <- colnames(X)
    num_predictors <- length(predictor_names)
    
    result_matrix <- matrix(NA, nrow = num_predictors, ncol = num_locations + 1)
    colnames(result_matrix) <- c(loop.detect, "Predictor") 
    
    tryCatch({
      data_source_type <- input$data_source
      if (data_source_type == "Average from MEB data") {
        all_coef_matrices <- list()
        for (iter in 1:num_locations) {
          coefs_df <- plotData_meb(
            location = routes[iter], 
            method = method_param,
            crit = selected_criterion,
            show_ = "Coefficients"
          )
          all_coef_matrices[[iter]] <- rowMeans(coefs_df[,-21])
        }
        combined_values <- do.call(cbind, all_coef_matrices)
        result_matrix[, 1:num_locations] <- combined_values
      } else if (data_source_type == "From Initial Data") {
        coefs_df <- plotData(
          location = paste0(input$route,"_volume"),
          method = method_param,
          crit = selected_criterion,
          meb_ = FALSE,
          show_ = "Coefficients"
        )
        result_df <- coefs_df
        return(result_df) 
      } else { 
        for (iter in 1:num_locations) {
          coefs_df <- plotData(
            location = routes[iter],
            method = method_param,
            crit = selected_criterion,
            meb_ = TRUE, 
            show_ = "Coefficients"
          )
          
          result_matrix[, iter] <- coefs_df[,1]
        }
      }
      result_matrix[, "Predictor"] <- predictor_names
      result_df <- as.data.frame(result_matrix)
      for(col in loop.detect) {
        result_df[[col]] <- as.numeric(result_df[[col]])
      }
      result_df[["Predictor"]] <- predictor_names
      validate(need(result_df, "Final result processing failed."))
      return(result_df)
      
    }, error = function(e) {
      showNotification(paste("Location-wise Coef Data Error:", e$message), type = "error", duration = 10)
      NULL 
    })
  })
  
  reactive_vip_data <- reactive({
    req(input$method, input$criterion, input$route)
    
    location_full_name <- paste0(input$route,"_volume")
    is_robust_param <- input$robust
    is_mte_param <- (input$method == "MTE")
    method_param <- c(is_robust_param, as.integer(is_mte_param))
    
    tryCatch({
      dt <- plotData_VIPs(
        location = location_full_name,
        method = method_param,
        crit = input$criterion
      )
      validate(
        need(dt, "VIP data loading function returned NULL."),
        need(is.data.frame(dt) && nrow(dt) > 0, "No VIP data returned.")
      )
      dt
    }, error = function(e) {
      showNotification(paste("VIP Load Error:", e$message), type = "error", duration=10)
      NULL
    })
  })
  
  reactive_res <- reactive({
    location <- paste0(input$route,"_volume")
    method <- c(input$robust, input$method=="MTE")
    crit <- input$criterion
    pDm <-cbind(plotData_meb(location, method, crit, show_="Residuals")[,1:20],plotData(location,method,crit, F,show_="Residuals")[,input$route],plotData(location,method,crit, T,show_="Residuals")[,1])
    names(pDm) <- c(paste0("MEB_",c(1:20)),"Actual", "ADL1LASSO")
    pDm
  })
  
  reactive_error <- reactive({
    location <- paste0(input$route,"_volume")
    method <- c(input$robust, input$method=="MTE")
    crit <- input$criterion
    pDm <- cbind(plotData_meb(location, method, crit, show_="Prediction Errors")[,1:20],plotData(location,method,crit, F,show_="Prediction Errors")[,input$route],plotData(location,method,crit, T,show_="Prediction Errors")[,1])
    names(pDm) <- c(paste0("MEB_",c(1:20)),"Actual", "ADL1LASSO")
    pDm
  })
  
  reactive_metrics <- reactive({
    location <- paste0(input$route,"_volume")
    split_loc <- input$route
    env <- new.env() 
    load(paste0(split_loc,"/synth_dataset.RData"),envir = env); meb_data <- env$MEB
    rm(env)
    train_rows <- 1:(dim(data)[1]/2)
    test_rows <- (dim(data)[1]/2 + 1):(dim(data)[1]*3/4)
    
    y_pred_tr <- cbind(plotData_meb(location, method, crit, show_="Fitted Values")[,-21], plotData(location,method,crit, F,show_="Fitted Values")[,split_loc], plotData(location,method,crit, T,show_="Fitted Values")[,1])
    y_pred_te <- cbind(plotData_meb(location, method, crit, show_="Predict Values")[,-21], plotData(location,method,crit, F,show_="Predict Values")[,split_loc], plotData(location,method,crit, T,show_="Predict Values")[,1])
    y_true_tr <- cbind(meb_data[train_rows,], data[train_rows,location], data[train_rows,location])
    y_true_te <- cbind(meb_data[test_rows,], data[test_rows,location], data[test_rows,location])
    
    mae_tr <- colMeans(abs(y_pred_tr-y_true_tr))
    rmse_tr <- sqrt(colMeans((y_pred_tr-y_true_tr)^2))
    huber_loss_tr <- mapply(
      function(col_true, col_pred) huber_loss(col_true, col_pred),
      as.data.frame(y_true_tr),
      as.data.frame(y_pred_tr)
    )
    mape_tr <- mapply(
      function(col_true, col_pred) mape(col_true, col_pred),
      as.data.frame(y_true_tr),
      as.data.frame(y_pred_tr)
    )
    mdape_tr <- mapply(
      function(col_true, col_pred) mdape(col_true, col_pred),
      as.data.frame(y_true_tr),
      as.data.frame(y_pred_tr)
    )
    smape_tr <- mapply(
      function(col_true, col_pred) smape(col_true, col_pred),
      as.data.frame(y_true_tr),
      as.data.frame(y_pred_tr)
    )
    smdape_tr <- mapply(
      function(col_true, col_pred) smdape(col_true, col_pred),
      as.data.frame(y_true_tr),
      as.data.frame(
        y_pred_tr)
    )
    
    
    mae_te <- colMeans(abs(y_pred_te-y_true_te))
    rmse_te <- sqrt(colMeans((y_pred_tr-y_true_tr)^2))
    huber_loss_te <- mapply(
      function(col_true, col_pred) huber_loss(col_true, col_pred),
      as.data.frame(y_true_te),
      as.data.frame(y_pred_te)
    )
    mape_te <- mapply(
      function(col_true, col_pred) mape(col_true, col_pred),
      as.data.frame(y_true_te),
      as.data.frame(y_pred_te)
    )
    mdape_te <- mapply(
      function(col_true, col_pred) mdape(col_true, col_pred),
      as.data.frame(y_true_te),
      as.data.frame(y_pred_te)
    )
    smape_te <- mapply(
      function(col_true, col_pred) smape(col_true, col_pred),
      as.data.frame(y_true_te),
      as.data.frame(y_pred_te)
    )
    smdape_te <- mapply(
      function(col_true, col_pred) smdape(col_true, col_pred),
      as.data.frame(y_true_te),
      as.data.frame(y_pred_te)
    )
    list(
      data = data.frame(
        Variable = c(rep("MAE_TRAIN",20),rep("MAE_TEST",20), rep("RMSE_TRAIN",20),rep("RMSE_TEST",20),
                     rep("HUBER_LOSS_TRAIN",20),rep("HUBER_LOSS_TEST",20), rep("MAPE_TRAIN",20),rep("MAPE_TEST",20),
                     rep("MdAPE_TRAIN",20),rep("MdAPE_TEST",20), rep("SMAPE_TRAIN",20),rep("SMAPE_TEST",20),
                     rep("SMdAPE_TRAIN",20),rep("SMdAPE_TEST",20)),
        Value = c(mae_tr[-(21:22)], mae_te[-(21:22)], rmse_tr[-(21:22)], rmse_te[-(21:22)],
                  huber_loss_tr[-(21:22)], huber_loss_te[-(21:22)], mape_tr[-(21:22)],mape_tr[-(21:22)],
                  mdape_tr[-(21:22)], mdape_te[-(21:22)],smape_tr[-(21:22)], smape_tr[-(21:22)],
                  smdape_tr[-(21:22)], smdape_te[-(21:22)])),
      artifacts = data.frame(
        Actual = c(mae_tr[21], mae_te[21], rmse_tr[21], rmse_te[21],
                   huber_loss_tr[21], huber_loss_te[21], mape_tr[21],mape_te[21],
                   mdape_tr[21], mdape_te[21],smape_tr[21], smape_te[21],
                   smdape_tr[21], smdape_te[21]),
        ADL1LASSO = c(mae_tr[22], mae_te[22], rmse_tr[22], rmse_te[22],
                      huber_loss_tr[22], huber_loss_te[22], mape_tr[22],mape_te[22],
                      mdape_tr[22], mdape_te[22],smape_tr[22], smape_te[22],
                      smdape_tr[22], smdape_te[22])
      ) 
    )
    
  })
  
  reactive_acf_pacf <- reactive({
    meb <- rowMeans(plotData_meb(location, method, crit, show_="Coefficients")[,1:20])
    actual <- plotData(location,method,crit, F,show_="Coefficients")[,input$route]
    adl1lasso <- plotData(location,method,crit, T,show_="Coefficients")[,1]
    dt <- data[,paste0(input$route,"_volume")]
    test_rows <- (length(dt)/2 + 1):(length(dt)*3/4)
    
    residuals_list <- data.frame(MEB = c(dt[train_rows]-meb%*%t(X[train_rows,])), Actual = c(dt[train_rows]-actual%*%t(X[train_rows,])), adl1lasso = c(dt[train_rows]-adl1lasso%*%t(X[train_rows,])))
    pred_errors_list <- data.frame(MEB = c(dt[test_rows]-meb%*%t(X[test_rows,])), Actual = c(dt[test_rows]-actual%*%t(X[test_rows,])), adl1lasso = c(dt[test_rows]-adl1lasso%*%t(X[test_rows,])))
    
    list(
      residuals = residuals_list,
      prediction_errors = pred_errors_list
    )
  })
  
  reactive_tests_results <- reactive({
    results <- list()
    residuals <- reactive_res()
    errors    <- reactive_error()
    for (col in colnames(residuals)) {
      res <- residuals[, col]
      pred <- errors[, col]
      
      normality_res <- list(
        shapiro_resid = shapiro.test(res), jberra_resid = jarque.bera.test(res), ksmirnov_resid = ks.test(res,"pnorm"),
        shapiro_pred = shapiro.test(pred), jberra_pred = jarque.bera.test(pred), ksmirnov_pred = ks.test(pred,"pnorm")
      )
      
      stationarity_res <- list(
        adf_resid = adf.test(res), kpss_res =kpss.test(res), pp_res =pp.test(res),
        adf_pred = adf.test(pred), kpss_pred =kpss.test(pred), pp_pred =pp.test(pred)
      )
      
      results[[col]] <- list(
        normality = normality_res,
        stationarity = stationarity_res
      )
    }
    
    return(results)
  })
  
  output$metrics_violin_plot <- renderPlotly({
    dist_data <- reactive_metrics()
    req(dist_data)
    plot_data <- dist_data$data
    CI95_mae_tr         <- round(quantile(plot_data[plot_data$Variable=="MAE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_mae_te         <- round(quantile(plot_data[plot_data$Variable=="MAE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_rmse_tr        <- round(quantile(plot_data[plot_data$Variable=="RMSE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_rmse_te        <- round(quantile(plot_data[plot_data$Variable=="RMSE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_huber_loss_tr  <- round(quantile(plot_data[plot_data$Variable=="HUBER_LOSS_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_huber_loss_te  <- round(quantile(plot_data[plot_data$Variable=="HUBER_LOSS_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_mape_tr        <- round(quantile(plot_data[plot_data$Variable=="MAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_mape_te        <- round(quantile(plot_data[plot_data$Variable=="MAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_mdape_tr       <- round(quantile(plot_data[plot_data$Variable=="MdAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_mdape_te       <- round(quantile(plot_data[plot_data$Variable=="MdAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_smape_tr       <- round(quantile(plot_data[plot_data$Variable=="SMAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_smape_te       <- round(quantile(plot_data[plot_data$Variable=="SMAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_smdape_tr      <- round(quantile(plot_data[plot_data$Variable=="SMdAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_smdape_te      <- round(quantile(plot_data[plot_data$Variable=="SMdAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
    
    artifacts <- data.frame(
      "MAE_TRAIN"         = c(CI95_mae_tr, mean(plot_data[plot_data$Variable=="MAE_TRAIN","Value"]), round(dist_data$artifacts$Actual[1],3), round(dist_data$artifacts$ADL1LASSO[1],3)),
      "MAE_TEST"          = c(CI95_mae_te, mean(plot_data[plot_data$Variable=="MAE_TEST","Value"]), round(dist_data$artifacts$Actual[2],3), round(dist_data$artifacts$ADL1LASSO[2],3)),
      "RMSE_TRAIN"        = c(CI95_rmse_tr, mean(plot_data[plot_data$Variable=="RMSE_TRAIN","Value"]), round(dist_data$artifacts$Actual[3],3), round(dist_data$artifacts$ADL1LASSO[3],3)),
      "RMSE_TEST"         = c(CI95_rmse_te, mean(plot_data[plot_data$Variable=="RMSE_TEST","Value"]), round(dist_data$artifacts$Actual[4],3), round(dist_data$artifacts$ADL1LASSO[4],3)),
      "HUBER_LOSS_TRAIN"  = c(CI95_huber_loss_tr, mean(plot_data[plot_data$Variable=="HUBER_LOSS_TRAIN","Value"]), round(dist_data$artifacts$Actual[5],3), round(dist_data$artifacts$ADL1LASSO[5],3)),
      "HUBER_LOSS_TEST"   = c(CI95_huber_loss_te, mean(plot_data[plot_data$Variable=="HUBER_LOSS_TEST","Value"]), round(dist_data$artifacts$Actual[6],3), round(dist_data$artifacts$ADL1LASSO[6],3)),
      "MAPE_TRAIN"        = c(CI95_mape_tr, mean(plot_data[plot_data$Variable=="MAPE_TRAIN","Value"]), round(dist_data$artifacts$Actual[7],3), round(dist_data$artifacts$ADL1LASSO[7],3)),
      "MAPE_TEST"         = c(CI95_mape_te, mean(plot_data[plot_data$Variable=="MAPE_TEST","Value"]), round(dist_data$artifacts$Actual[8],3), round(dist_data$artifacts$ADL1LASSO[8],3)),
      "MdAPE_TRAIN"       = c(CI95_mdape_tr, mean(plot_data[plot_data$Variable=="MdAPE_TRAIN","Value"]), round(dist_data$artifacts$Actual[9],3), round(dist_data$artifacts$ADL1LASSO[9],3)),
      "MdAPE_TEST"        = c(CI95_mdape_te, mean(plot_data[plot_data$Variable=="MdAPE_TEST","Value"]), round(dist_data$artifacts$Actual[10],3), round(dist_data$artifacts$ADL1LASSO[10],3)),
      "SMAPE_TRAIN"       = c(CI95_smape_tr, mean(plot_data[plot_data$Variable=="SMAPE_TRAIN","Value"]), round(dist_data$artifacts$Actual[11],3), round(dist_data$artifacts$ADL1LASSO[11],3)),
      "SMAPE_TEST"        = c(CI95_smape_te, mean(plot_data[plot_data$Variable=="SMAPE_TEST","Value"]), round(dist_data$artifacts$Actual[12],3), round(dist_data$artifacts$ADL1LASSO[12],3)),
      "SMdAPE_TRAIN"      = c(CI95_smdape_tr, mean(plot_data[plot_data$Variable=="SMdAPE_TRAIN","Value"]), round(dist_data$artifacts$Actual[13],3), round(dist_data$artifacts$ADL1LASSO[13],3)),
      "SMdAPE_TEST"       = c(CI95_smdape_te, mean(plot_data[plot_data$Variable=="SMdAPE_TEST","Value"]), round(dist_data$artifacts$Actual[14],3), round(dist_data$artifacts$ADL1LASSO[14],3))
    )
    
    artifacts_long <- artifacts %>%
      t() %>%
      as.data.frame() %>%
      mutate(Variable = rownames(.)) %>%
      rename(
        `5%` = V1,
        `50%` = V2,
        `95%` = V3,
        mean = V4,
        Actual = V5,
        ADL1LASSO = V6
      )
    var_levels <- unique(plot_data$Variable)
    var_map <- setNames(seq_along(var_levels), var_levels)
    plot_data <- plot_data %>%
      mutate(x = var_map[Variable])
    artifacts_long <- artifacts_long %>%
      mutate(x = var_map[Variable])
    tryCatch({
      p <- plot_ly()
      for (var in var_levels) {
        p <- p %>%
          add_trace(
            type = "violin",
            y = plot_data$Value[plot_data$Variable == var],
            x = rep(var_map[var], sum(plot_data$Variable == var)),
            name = var,
            box = list(visible = FALSE),
            meanline = list(visible = FALSE),
            points = 'all',
            opacity = 0.6,
            spanmode = 'hard'
          )
      }
      p <- p %>%
        add_segments(
          x = artifacts_long$x,
          xend = artifacts_long$x,
          y = artifacts_long$`5%`,
          yend = artifacts_long$`95%`,
          line = list(color = 'black', width = 2),
          name = "95% CI"
        )
      p <- p %>%
        add_segments(
          x = artifacts_long$x - 0.15,
          xend = artifacts_long$x + 0.15,
          y = artifacts_long$`50%`,
          yend = artifacts_long$`50%`,
          line = list(color = '#D7263D', width = 3),
          name = "Median"
        )
      p <- p %>%
        add_segments(
          x = artifacts_long$x - 0.15,
          xend = artifacts_long$x + 0.15,
          y = artifacts_long$mean,
          yend = artifacts_long$mean,
          line = list(color = '#F46036', width = 3),
          name = "Mean"
        )
      p <- p %>%
        add_markers(
          x = artifacts_long$x,
          y = artifacts_long$Actual,
          marker = list(color = '#1B2A41', size = 10, symbol = 'diamond'),
          name = "Actual"
        )
      p <- p %>%
        add_markers(
          x = artifacts_long$x,
          y = artifacts_long$ADL1LASSO,
          marker = list(color = '#247BA0', size = 10, symbol = 'diamond'),
          name = "ADL1LASSO"
        )
      p <- p %>%
        layout(
          title = "Distribution of all the metrics (Train/Test) with Artifacts",
          yaxis = list(title = "Error Value"),
          xaxis = list(
            title = "",
            tickmode = 'array',
            tickvals = seq_along(var_levels),
            ticktext = var_levels
          ),
          showlegend = TRUE
        )
      return(p)
      
    }, error = function(e) {
      showNotification(paste("Error generating violin plot:", e$message), type = "error")
      plot_ly() %>% layout(title = "Error Generating Plot")
    })
  })
  reactive_VIPs_tot <- reactive({
    req(input$method, input$criterion, input$route)
    location <- paste0(input$route,"_volume")
    coeff_names <- c("(Intercept)", regressors)
    res <- data.frame(
      apply(
        cbind(
          "L1MTE_AIC" = c(T,T,"AIC"),
          "L1MTE_BIC" = c(T,T,"BIC"),
          "L1LSA_AIC" = c(T,F,"AIC"),
          "L1LSA_BIC" = c(T,F,"BIC"),
          "L2MTE_AIC" = c(F,T,"AIC"),
          "L2MTE_BIC" = c(F,T,"BIC"),
          "L2LSA_AIC" = c(F,F,"AIC"),
          "L2LSA_BIC" = c(F,F,"BIC")), 2,
        function(x) plotData_VIPs(location,
                                  c(as.logical(x[1]),as.logical(x[2])),
                                  x[3])[,1]),
      "Predictor" = coeff_names)
    res <- res %>%
      arrange(across(1:8, ~ desc(.)))
    res[rowSums(res[,-9])>0,]
  })
  reactive_vip_data <- reactive({
    req(input$method, input$criterion, input$route)
    location_full_name <- paste0(input$route,"_volume")
    is_robust_param <- input$robust
    is_mte_param <- (input$method == "MTE")
    method_param <- c(is_robust_param, as.integer(is_mte_param))
    tryCatch({
      dt <- plotData_VIPs(
        location = location_full_name,
        method = method_param,
        crit = input$criterion
      )
      validate(
        need(dt, "VIP data loading function returned NULL."),
        need(is.data.frame(dt) && nrow(dt) > 0, "No VIP data returned.")
      )
      dt
    }, error = function(e) {
      showNotification(paste("VIP Load Error:", e$message), type = "error", duration=10)
      NULL
    })
  })
  output$artifacts_table <- renderDT({
    dist_data <- reactive_metrics()
    req(dist_data)
    plot_data <- dist_data$data
    CI95_mae_tr         <- round(quantile(plot_data[plot_data$Variable=="MAE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_mae_te         <- round(quantile(plot_data[plot_data$Variable=="MAE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_rmse_tr        <- round(quantile(plot_data[plot_data$Variable=="RMSE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_rmse_te        <- round(quantile(plot_data[plot_data$Variable=="RMSE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_huber_loss_tr  <- round(quantile(plot_data[plot_data$Variable=="HUBER_LOSS_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_huber_loss_te  <- round(quantile(plot_data[plot_data$Variable=="HUBER_LOSS_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_mape_tr        <- round(quantile(plot_data[plot_data$Variable=="MAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_mape_te        <- round(quantile(plot_data[plot_data$Variable=="MAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_mdape_tr       <- round(quantile(plot_data[plot_data$Variable=="MdAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_mdape_te       <- round(quantile(plot_data[plot_data$Variable=="MdAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_smape_tr       <- round(quantile(plot_data[plot_data$Variable=="SMAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_smape_te       <- round(quantile(plot_data[plot_data$Variable=="SMAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
    CI95_smdape_tr      <- round(quantile(plot_data[plot_data$Variable=="SMdAPE_TRAIN","Value"],c(0.05,0.5,0.95)),3)
    CI95_smdape_te      <- round(quantile(plot_data[plot_data$Variable=="SMdAPE_TEST","Value"],c(0.05,0.5,0.95)),3)
    
    artifacts <- data.frame(
      "MAE_TRAIN"         = c(CI95_mae_tr, mean(plot_data[plot_data$Variable=="MAE_TRAIN","Value"]), round(dist_data$artifacts$Actual[1],3), round(dist_data$artifacts$ADL1LASSO[1],3)),
      "MAE_TEST"          = c(CI95_mae_te, mean(plot_data[plot_data$Variable=="MAE_TEST","Value"]), round(dist_data$artifacts$Actual[2],3), round(dist_data$artifacts$ADL1LASSO[2],3)),
      "RMSE_TRAIN"        = c(CI95_rmse_tr, mean(plot_data[plot_data$Variable=="RMSE_TRAIN","Value"]), round(dist_data$artifacts$Actual[3],3), round(dist_data$artifacts$ADL1LASSO[3],3)),
      "RMSE_TEST"         = c(CI95_rmse_te, mean(plot_data[plot_data$Variable=="RMSE_TEST","Value"]), round(dist_data$artifacts$Actual[4],3), round(dist_data$artifacts$ADL1LASSO[4],3)),
      "HUBER_LOSS_TRAIN"  = c(CI95_huber_loss_tr, mean(plot_data[plot_data$Variable=="HUBER_LOSS_TRAIN","Value"]), round(dist_data$artifacts$Actual[5],3), round(dist_data$artifacts$ADL1LASSO[5],3)),
      "HUBER_LOSS_TEST"   = c(CI95_huber_loss_te, mean(plot_data[plot_data$Variable=="HUBER_LOSS_TEST","Value"]), round(dist_data$artifacts$Actual[6],3), round(dist_data$artifacts$ADL1LASSO[6],3)),
      "MAPE_TRAIN"        = c(CI95_mape_tr, mean(plot_data[plot_data$Variable=="MAPE_TRAIN","Value"]), round(dist_data$artifacts$Actual[7],3), round(dist_data$artifacts$ADL1LASSO[7],3)),
      "MAPE_TEST"         = c(CI95_mape_te, mean(plot_data[plot_data$Variable=="MAPE_TEST","Value"]), round(dist_data$artifacts$Actual[8],3), round(dist_data$artifacts$ADL1LASSO[8],3)),
      "MdAPE_TRAIN"       = c(CI95_mdape_tr, mean(plot_data[plot_data$Variable=="MdAPE_TRAIN","Value"]), round(dist_data$artifacts$Actual[9],3), round(dist_data$artifacts$ADL1LASSO[9],3)),
      "MdAPE_TEST"        = c(CI95_mdape_te, mean(plot_data[plot_data$Variable=="MdAPE_TEST","Value"]), round(dist_data$artifacts$Actual[10],3), round(dist_data$artifacts$ADL1LASSO[10],3)),
      "SMAPE_TRAIN"       = c(CI95_smape_tr, mean(plot_data[plot_data$Variable=="SMAPE_TRAIN","Value"]), round(dist_data$artifacts$Actual[11],3), round(dist_data$artifacts$ADL1LASSO[11],3)),
      "SMAPE_TEST"        = c(CI95_smape_te, mean(plot_data[plot_data$Variable=="SMAPE_TEST","Value"]), round(dist_data$artifacts$Actual[12],3), round(dist_data$artifacts$ADL1LASSO[12],3)),
      "SMdAPE_TRAIN"      = c(CI95_smdape_tr, mean(plot_data[plot_data$Variable=="SMdAPE_TRAIN","Value"]), round(dist_data$artifacts$Actual[13],3), round(dist_data$artifacts$ADL1LASSO[13],3)),
      "SMdAPE_TEST"       = c(CI95_smdape_te, mean(plot_data[plot_data$Variable=="SMdAPE_TEST","Value"]), round(dist_data$artifacts$Actual[14],3), round(dist_data$artifacts$ADL1LASSO[14],3))
    )
    
    artifacts_long <- artifacts %>%
      t() %>%
      as.data.frame() %>%
      mutate(Variable = rownames(.)) %>%
      rename(
        `5%` = V1,
        `50%` = V2,
        `95%` = V3,
        mean = V4,
        Actual = V5,
        ADL1LASSO = V6
      )
    var_levels <- unique(plot_data$Variable)
    var_map <- setNames(seq_along(var_levels), var_levels)
    plot_data <- plot_data %>%
      mutate(x = var_map[Variable])
    artifacts_long <- artifacts_long %>%
      mutate(x = var_map[Variable])
    datatable(
      artifacts_long %>% 
        select(-x),  # Remove the numeric x column
      options = list(pageLength = 5, scrollX = TRUE)
    )
  })
  
  output$tests_summary <- renderPrint({
    results <- reactive_tests_results()
    req(results)
    
    for (var in names(results)) {
      cat("\n=================================================================\n ")
      cat(paste("Variable:", var, "\n"))
      cat("-------------------------------------------------------------------\n")
      cat("Normality Test (Residuals):\n")
      print(results[[var]]$normality$shapiro_resid)
      cat("\nNormality Test (Prediction Errors):\n")
      print(results[[var]]$normality$shapiro_pred)
      cat("\nStationarity Test (Residuals):\n")
      print(results[[var]]$stationarity$adf_resid)
      cat("\nStationarity Test (Prediction Errors):\n")
      print(results[[var]]$stationarity$adf_pred)
      cat("\n")
    }
  })
  
  output$acf_pacf_plot <- renderPlotly({
    acf_data <- reactive_acf_pacf()
    req(acf_data)
    
    plots_list <- list()
    n_vars <- length(acf_data$residuals)
    res_df <- data.frame(acf_data$residuals)
    pred_df <- data.frame(acf_data$prediction_errors)
    res2acf  <- cbind(res_df[,1],res_df[,2],res_df[,3])
    pred2acf <- cbind(pred_df[,1],pred_df[,2],pred_df[,3])
    names(res2acf) <- names(pred2acf) <- c("MEB_AVG", "ACTUAL", "ADL1LASSO")
    for (i in 1:3) {
      acf_res <- acf(res2acf[,i], plot = FALSE)
      pacf_res <- pacf(res2acf[,i], plot = FALSE)
      
      acf_pred <- acf(pred2acf[,i], plot = FALSE)
      pacf_pred <- pacf(pred2acf[,i], plot = FALSE)
      
      plots_list[[length(plots_list) + 1]] <- ggplotly(autoplot(acf_res) + labs(y =paste0("ACF ", names(res2acf)[i],"(Train)")))
      plots_list[[length(plots_list) + 1]] <- ggplotly(autoplot(pacf_res) + labs(y =paste0("PACF ", names(res2acf)[i],"(Train)")))
      plots_list[[length(plots_list) + 1]] <- ggplotly(autoplot(acf_pred) + labs(y =paste0("ACF  ", names(pred2acf)[i],"Test")))
      plots_list[[length(plots_list) + 1]] <- ggplotly(autoplot(pacf_pred) + labs(y =paste0("PACF ", names(pred2acf)[i],"(Test)")))
    }
    
    # Arrange all plots in a grid
    subplot(plots_list, nrows = ceiling(length(plots_list)/4), margin = 0.02, shareX = TRUE, shareY = FALSE, titleY = TRUE) %>%
      layout(title = "ACF and PACF Plots for Residuals and Prediction Errors")  
    })
  
  output$plot_title <- renderText({
    req(input$plot_type, input$route)
    paste(input$plot_type, "for Detector:", input$route) 
  })
  
  output$main_plot <- renderPlotly({
    df <- reactive_plot_data()
    req(df)
    current_plot_type <- req(input$plot_type)
    
    p <- tryCatch({
      if (!("Predictor" %in% names(df))) {
        stop("Data frame missing 'Predictor' column.")
      }
      
      if (current_plot_type == "Coefficients") {
        dt_long <- reshape2::melt(
          df,
          id.vars = "Predictor",
          variable.name = "Series",
          value.name = "MEB"
        )
        dt_long <- dt_long %>%
          arrange(desc(MEB))
        
        dt_long$Predictor <- factor(dt_long$Predictor, levels = unique(dt_long$Predictor))
        p <- tryCatch({
          plot_ly(
            dt_long,
            x = ~Predictor,
            y = ~MEB,
            color = ~Series,
            colors = colorRampPalette(Polychrome::palette36.colors(36))(length(unique(dt_long$Series))),
            type = "bar",
            hoverinfo = "text",
            text = ~paste('Predictor:', Predictor, '<br>Series:', Series, '<br>Value:', round(MEB, 4))
          ) |>
            layout(
              title = "",
              xaxis = list(title = "Predictor Feature", type = "category", tickangle = 45),
              yaxis = list(title = "Coefficient Value"),
              showlegend = TRUE
            )
        }, error = function(e) {
          showNotification(paste("Plot Error:", e$message), type = "warning", duration = 5)
          plot_ly() |> layout(title = "Error Generating Plot", annotations = list(text = paste("Details:", e$message), showarrow = FALSE))
        })
      } else if (ncol(df) >= 2) {
        
        if (ncol(df) > 2) {
          melted <- reshape2::melt(df, id.vars = "Predictor", variable.name = "Series", value.name = "Value")
          melted$Series <- factor(melted$Series) 
          n_series <- length(unique(melted$Series))
          palette_dynamic <- colorRampPalette(palette36.colors(36))(n_series)
          
          plot_ly(
            melted,
            x = ~Predictor,
            y = ~Value,
            color = ~Series,
            colors = palette_dynamic,
            type = "scatter",
            mode = "lines",
            hoverinfo = 'text',
            text = ~paste('Time/Index:', Predictor, '<br>Series:', Series, '<br>Value:', round(Value, 4))
          ) |>
            layout(
              title = "",
              xaxis = list(title = "Time / Observation Index"),
              yaxis = list(title = "Value")
            )
          
        } else {
          plot_ly(
            df,
            x = ~Predictor,
            y = ~Value,
            type = "scatter",
            mode = "lines",
            hoverinfo = 'text',
            text = ~paste('Time/Index:', Predictor, '<br>Value:', round(Value, 4))
          ) |>
            layout(
              title = "",
              xaxis = list(title = "Time / Observation Index"),
              yaxis = list(title = "Value")
            )
        }
        
      } else {
        stop("Plotting error: Data frame has insufficient columns.")
      }
    }, error = function(e) {
      showNotification(paste("Plotting Error:", e$message), type = "warning", duration=5)
      plot_ly() |> layout(title = "Error Generating Plot", annotations = list(text = paste("Details:", e$message), showarrow=FALSE))
    })
    p
  })
  
  output$MEB_vs_ACTUAL_title <- renderText({
    req(input$plot_type, input$route)
    paste(input$plot_type, "for Detector:", input$route) 
  })
  
  output$comparison_plot <- renderPlotly({
    df <- plotData_meb(
      location = paste0(input$route, "_volume"),
      method = c(input$robust, (input$method == "MTE")),
      crit = input$criterion,
      show_ = "Coefficients"
    )
    
    act_df <- plotData(
      location = paste0(input$route, "_volume"),
      method = c(input$robust, (input$method == "MTE")),
      crit = input$criterion,
      meb_ = FALSE,
      show_ = "Coefficients"
    )[, input$route]
    
    ad_l1_lasso <- plotData(
      location = paste0(input$route, "_volume"),
      method = c(input$robust, (input$method == "MTE")),
      crit = input$criterion,
      meb_ = TRUE,
      show_ = "Coefficients"
    )
    
    req(df, act_df, ad_l1_lasso)
    
    mean_meb <- if (ncol(df) > 2) {
      rowMeans(df[, 1:(ncol(df) - 1)], na.rm = TRUE)
    } else {
      df[, 1]
    }
    
    posthoc_adl1lasso <- if (ncol(ad_l1_lasso) > 2) {
      rowMeans(ad_l1_lasso[, 1:(ncol(ad_l1_lasso) - 1)], na.rm = TRUE)
    } else {
      ad_l1_lasso[, 1]
    }
    
    df_combined <- data.frame(
      Mean_MEB = mean_meb,
      Actual = act_df,
      PostHoc_ADL1LASSO = posthoc_adl1lasso,
      Predictor = df$Predictor
    )
    
    names(df_combined) <- c("Mean_MEB", "Actual", "PostHoc_ADL1LASSO", "Predictor")
    
    tryCatch({
      if (!("Predictor" %in% names(df_combined))) {
        stop("Data frame missing 'Predictor' column.")
      }
      
      df_long <- reshape2::melt(
        df_combined,
        id.vars = "Predictor",
        variable.name = "Series",
        value.name = "Value"
      )
      
      df_long$Predictor <- factor(df_long$Predictor, levels = unique(df_long$Predictor))
      
      plot_ly(
        data = df_long,
        x = ~Predictor,
        y = ~Value,
        color = ~Series,
        colors = colorRampPalette(Polychrome::palette36.colors(36))(length(unique(df_long$Series))),
        type = "bar",
        hoverinfo = "text",
        text = ~paste('Predictor:', Predictor, '<br>Series:', Series, '<br>Value:', round(Value, 4))
      ) |>
        layout(
          title = "Mean MEB vs Actual vs PostHoc AD-L1 Lasso",
          xaxis = list(title = "Predictor", type = "category", tickangle = 45),
          yaxis = list(title = "Coefficient Value"),
          barmode = "group",
          showlegend = TRUE
        )
      
    }, error = function(e) {
      showNotification(paste("Plot Error:", e$message), type = "warning", duration = 5)
      plot_ly() |> layout(title = "Error Generating Plot", annotations = list(text = paste("Details:", e$message), showarrow = FALSE))
    })
  })
  
  output$locationWise_plot <- renderPlotly({
    df <- reactive_plotData_loc() 
    req(df) 
    tryCatch({
      if (!"Predictor" %in% names(df)) {
        stop("Data frame must have a 'Predictor' column.")
      }
      predictor_col <- "Predictor"
      value_cols <- names(df)[names(df) != predictor_col]
      if(length(value_cols) == 0) stop("No value columns found besides Predictor.")
      
      if (nrow(df) == 0) {
        stop("Data frame has no rows.")
      }
      dt_long <- reshape2::melt(
        df,
        id.vars = predictor_col,
        measure.vars = value_cols,
        variable.name = "Series", 
        value.name = "Value"      
      )
      if (is.factor(df[[predictor_col]])) {
        dt_long[[predictor_col]] <- factor(dt_long[[predictor_col]], levels = levels(df[[predictor_col]]))
      } else {
        dt_long[[predictor_col]] <- factor(dt_long[[predictor_col]], levels = unique(df[[predictor_col]]))
      }
      num_series <- length(value_cols)
      p <- plot_ly(
        data = dt_long,
        x = ~get(predictor_col),
        y = ~Value,
        color = ~Series,
        type = "bar",
        hoverinfo = "text",
        text = ~paste('Predictor:', get(predictor_col), '<br>Series:', Series, '<br>Value:', round(Value, 4))
      ) |>
        layout(
          title = paste("Coefficients across Locations:", input$data_source_loc), # Use correct input ID
          xaxis = list(
            title = predictor_col,
            type = "category",
            tickangle = -45
          ),
          yaxis = list(title = "Coefficient Value"),
          barmode = 'group',
          legend = list(title=list(text='Location')),
          showlegend = TRUE
        )
      
      return(p)
      
    }, error = function(e) {
      # --- Error Handling ---
      showNotification(paste("Location Plot Error:", e$message), type = "error", duration = 10)
      plot_ly() |> layout(title = "Error Generating Location Plot",
                          annotations = list(text = paste("Details:", e$message), showarrow = FALSE))
    })
  })
  
  output$vip_title <- renderText({
    req(input$route)
    paste("Variable Importance (VIP) Scores for Detector:", input$route)
  })
  
  output$vip_plot <- renderPlotly({
    dt <- reactive_vip_data()
    req(dt)
    
    p <- tryCatch({
      if (!all(c("Predictor", "Values") %in% names(dt))) stop(paste0("VIP data missing 'Predictor' or 'Value' column.\n The names are ",paste0(names(dt),collapse=", ")))
      
      dt <- dt[order(dt$Values, decreasing = TRUE), ]
      dt$Predictor <- factor(dt$Predictor, levels = dt$Predictor)
      
      plot_ly(
        data = dt,
        x = ~Predictor,
        y = ~Values,
        type = "bar",
        marker = list(color = "steelblue"),  
        hoverinfo = "text",
        text = ~paste('Feature:', Predictor, '<br>VIP Score:', round(Values, 4))
      ) |>
        layout(
          title = "",
          xaxis = list(title = "Features", type = "category"),
          yaxis = list(title = "VIP Score"),
          showlegend = FALSE  
        )
      
    }, error = function(e) {
      showNotification(paste("VIP Plot Error:", e$message), type = "warning", duration = 5)
      plot_ly() |> layout(title = "Error Generating VIP Plot", annotations = list(text = paste("Details:", e$message), showarrow = FALSE))
    })
    
    p
  })
  
  output$heatmap_plot <- renderPlotly({
    mat <- reactive_VIPs_tot() 
    plot_ly(
      x = mat[,9],
      y = colnames(mat[,-9]),
      z = t(mat),  # transpose to flip axes
      type = "heatmap",
      colors = colorRamp(c("navyblue", "white", "coral"))
    )
  })
  
  output$map_output <- renderLeaflet({
    leaflet(df_locations) |>
      addTiles(group = "OpenStreetMap") |>
      addProviderTiles(providers$CartoDB.Positron, group = "CartoDB (Light)") |>
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") |>
      addMarkers(
        lng = ~lng, lat = ~lat,
        popup = ~paste("<b>Detector:</b>", name),
        label = ~name,
        group = "Detectors"
      ) |>
      addLayersControl(
        baseGroups = c("OpenStreetMap", "CartoDB (Light)", "Esri Satellite"),
        overlayGroups = c("Detectors"),
        options = layersControlOptions(collapsed = FALSE)
      ) |>
      addMeasure(
        position = "bottomleft",
        primaryLengthUnit = "meters",
        primaryAreaUnit = "sqmeters",
        activeColor = "#3D535D",
        completedColor = "#7D4479") |>
      addScaleBar(position = "bottomright")
  })
  
  observeEvent(input$route, {
    req(input$route)
    selected_loc_name <- input$route
    loc_data <- df_locations[df_locations$name == selected_loc_name, ]
    
    if (nrow(loc_data) == 1) {
      leafletProxy("map_output", data = loc_data) |>
        flyTo(lng = loc_data$lng, lat = loc_data$lat, zoom = 16) |>
        clearGroup("Pulse") |>
        addPulseMarkers(
          lng = loc_data$lng,
          lat = loc_data$lat,
          label = loc_data$name,
          icon = makePulseIcon(heartbeat = 1.5, color = "red"),
          group = "Pulse"
        )
    } else {
      showNotification(paste("Location data not found for:", selected_loc_name), type="warning")
      leafletProxy("map_output") |> clearGroup("Pulse")
    }
  })
}

shinyApp(ui = ui, server = server)

# library(getip)
# runApp("shinyApp.R",host = "139.91.62.83",port = 1997) # IPv4 Address. 192.168.1.22. the uoc server is 147.52.205.205
