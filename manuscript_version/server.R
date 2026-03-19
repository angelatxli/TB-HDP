## Load Packages -----
library(shiny)
library(bslib)
library(shinyjs)
library(shinydashboard)
library(purrr)
library(markdown)

####
animals <- c("BALB/c Mouse", "Rat", "Dog")

species_params <- data.frame(
  Species = c("Human", "BALB/c Mouse", "Rat", "Dog"),
  BW = c(60, 0.02, 0.25, 10),
  LiverWeight = c(25.7, 80, 40, 25),  # liver weight (g/kg)
  MPPGL = c(40, 50, 45, 30), # microsomal protein per gram liver (mg/g liver)
  HPGL = c(120, 135, 120, 100), # million cells per gram liver
  Qh = c(20.7, 126, 77.4, 56.1) # hepatic blood flow (mL/min/kg)
  )

source("./helpers/allometric_scaling.R")
source("./helpers/ivive.R")
source("./helpers/mrgsolve.R")

server <- function(input, output, session) {
  
  observeEvent(input$tabs, {
    # List of all your notification IDs
    ids <- c("warn_lo", "warn_sp", "warn_pc", "pc_msg")
    
    # Remove all notifications whenever the user switches tabs
    # This keeps the 'About' and 'Citation' pages clean
    for (id in ids) {
      removeNotification(id)
    }
  })
  
# -------- LO -----
  
  ## Reset ----
  observeEvent(input$reset_lo, {
    # --- Compound Info ---
    updateTextInput(session,"drugname_lo", value = NA)

    # --- Simulation Params ---
    updateNumericInput(session, "ka_lo", value = NA)
    updateNumericInput(session, "dose_lo", value = NA)
    updateNumericInput(session, "ndoses_lo", value = NA)
    updateNumericInput(session, "inter_lo", value = NA)
    
    updateNumericInput(session, "fu_lo", value = NA)
    updateNumericInput(session, "heppk_lo", value = NA)
    updateNumericInput(session, "micpk_lo", value = NA)
    
    updateNumericInput(session, "preclin_cl1_lo", value = NA)
    updateNumericInput(session, "preclin_vss1_lo", value = NA)
    updateNumericInput(session, "preclin_cl2_lo", value = NA)
    updateNumericInput(session, "preclin_vss2_lo", value = NA)
    updateNumericInput(session, "preclin_cl3_lo", value = NA)
    updateNumericInput(session, "preclin_vss3_lo", value = NA)
    
    # --- PD Parameters ---
    updateNumericInput(session, "MIC", value = NA)

    plot3_data(NULL)
    plot6_data(NULL)

    showNotification("Parameters reset.", 
                     type = "default",
                     id = "reset_lo")
  })
  
  ## Sample params ----
  observeEvent(input$sample_params_lo, {
    # --- Compound Info ---
    updateTextInput(session, "drugname_lo", value = "Demo Compound")

    # --- Simulation Parameters ---
    updateNumericInput(session, "ka_lo", value = 0.5)
    updateNumericInput(session, "dose_lo", value = 200)
    updateNumericInput(session, "ndoses_lo", value = 14)
    updateNumericInput(session, "inter_lo", value = 24)
    
    # --- PK Parameters ---
    updateNumericInput(session, "fu_lo", value = 0.8)
    updateNumericInput(session, "heppk_lo", value = 12)
    updateNumericInput(session, "micpk_lo", value = 12)
    
    # --- Preclinical PK ---
    updateNumericInput(session, "preclin_cl1_lo", value = 20)
    updateNumericInput(session, "preclin_vss1_lo", value = 45)
    updateNumericInput(session, "preclin_cl2_lo", value = 40)
    updateNumericInput(session, "preclin_vss2_lo", value = 30)
    updateNumericInput(session, "preclin_cl3_lo", value = 50)
    updateNumericInput(session, "preclin_vss3_lo", value = 56)
    
    # --- PD Parameters ---
    updateNumericInput(session, "MIC", value = 0.001)
    
    # Send toast message to let the user know it worked
    showNotification("Sample parameters for Demo Compound loaded.", 
                     type = "default",
                     id = "reset_lo")
    })
  
  
  ## Reactivity ----
  # Create a reactive that bundles all inputs for the LO tab
  inputs_lo <- reactive({
    list(input$drugbackbone,
         input$pk_route_lo, input$ka_lo,
         input$dose_lo, input$ndoses_lo, input$inter_lo, 
         input$pkmethod_lo, 
         input$heppk_lo, input$fu_lo, input$micpk_lo,
         input$preclin_cl1_lo, input$preclin_cl2_lo, input$preclin_cl3_lo,
         input$preclin_vss1_lo, input$preclin_vss2_lo, input$preclin_vss3_lo,
         input$MIC)
    })
  
  # When any input changes, gray out the plot
  observeEvent(inputs_lo(), {
    shinyjs::addClass(id = "plot3_lo", class = "out-of-sync")
    shinyjs::addClass(id = "plot6_lo", class = "out-of-sync")
    shinyjs::addClass(id = "param_table_lo", class = "out-of-sync")
    shinyjs::addClass(id = "result_table_lo", class = "out-of-sync")
    
    if (input$tabs == "pred_lo") {
      showNotification("Input changed. Please update plots.", 
                       type = "default",
                       id = "input_lo")
      }
    })
  
  # When the button is clicked, remove gray effect
  observeEvent(input$click_lo, {
    shinyjs::removeClass(id = "plot3_lo", class = "out-of-sync")
    shinyjs::removeClass(id = "plot6_lo", class = "out-of-sync")
    shinyjs::removeClass(id = "param_table_lo", class = "out-of-sync")
    shinyjs::removeClass(id = "result_table_lo", class = "out-of-sync")
    })
  
  
  observeEvent(input$click_lo, {
    if (!is.na(input$fu_lo) && input$fu_lo < 0.001) {
      showNotification(
        "WARNING: Compound is highly protein bound (fu < 0.1%). 
        Unbound plasma PK/MIC may be unreliable; please continue to lesion studies.",
        type = "warning",
        id = "warn_lo",
        duration = 10
      )
      }
  })
  
  observe({
    # 1. Define required inputs and their user-friendly names
    reqs <- list(
      "Dose" = input$dose_lo,
      "Number of doses" = input$ndoses_lo,
      "Dosing interval" = input$inter_lo,
      "MIC" = input$MIC
    )
    
    # 2. Check for NA or zero values
    missing <- names(reqs)[sapply(reqs, function(x) is.null(x) || is.na(x) || x <= 0)]
    
    if (length(missing) > 0) {
      
      shinyjs::disable("click_lo") # Disable button
      
      # Create message
      output$lo_validation_ui <- renderUI({
        div(style = "background-color: #fcf8e3; border: 1px solid #faebcc; padding: 10px; border-radius: 4px; margin-bottom: 10px;",
            tags$b(style = "color: #8a6d3b;", icon("exclamation-triangle"), " Please enter:"),
            tags$ul(style = "color: #8a6d3b; margin-top: 5px; margin-bottom: 0;",
                    lapply(missing, function(m) tags$li(m))
            )
        )
      })
      shinyjs::addClass("click_lo", "btn-outline-secondary")
      shinyjs::removeClass("click_lo", "btn-primary")
    } 
    
    else {
      # Enable button
      shinyjs::enable("click_lo")
      output$lo_validation_ui <- renderUI({
        div(style = "background-color: #dff0d8; border: 1px solid #d6e9c6; padding: 10px; border-radius: 4px; margin-bottom: 10px;",
            tags$span(style = "color: #3c763d;", icon("check-circle"), " Ready to simulate.")
        )
      })
      
      # Clear tooltip and restore primary color
      shinyjs::addClass("click_lo", "btn-primary")
      shinyjs::removeClass("click_lo", "btn-outline-secondary")
    }
  })
  
  ## --- Clearance ----
  
  ### ---- Single species allometric scaling
  BW <- reactive({
    species_params %>% filter(Species == input$species1_lo) %>% pull(BW)
    })
  
  ### ---- Multi-species allometric scaling
  df_lo <- reactive({
    species <- c(input$species1_lo, 
                 input$species2_lo, 
                 input$species3_lo)
    
    cls <- c(input$preclin_cl1_lo, 
             input$preclin_cl2_lo, 
             input$preclin_cl3_lo)
    
    vs <- c(input$preclin_vss1_lo, 
            input$preclin_vss2_lo, 
            input$preclin_vss3_lo)
    
    weights_df <- species_params %>% filter(Species %in% species)
    weights <- weights_df$BW[match(species, weights_df$Species)]
    
    df_temp <- data.frame(CL = cls, V = vs, Species = species, BodyWeight = weights)
    df_temp$logBW <- log10(df_temp$BodyWeight)
    df_temp$logCL <- log10(df_temp$CL)
    df_temp$logV <- log10(df_temp$V)
    
    df_temp
    })
  
  model_Cl_results_lo <- reactive({
    model_cl_fit(df_lo())
    })
  
  model_V_results_lo <- reactive({
    model_v_fit(df_lo())
    })
  
  #### Plot --------
  plot_lo_obj <- reactive({
    res <- model_Cl_results_lo()
    plot_log_fit(df_lo(), res, "logCL", "log10 Clearance (L/h)")
    })
  
  output$plot_lo <- renderPlot({
    plot_lo_obj()
    })
  
  plot2_lo_obj <- reactive({
    res <- model_V_results_lo()
    plot_log_fit(df_lo(), res, "logV", "log10 Volume (L)")
    })
  
  output$plot2_lo <- renderPlot({
    plot2_lo_obj()
    })
  
  clin_CL_lo <- reactive({
    if (input$pkmethod_lo == "ivive_h") {
      hep_ivive(input$heppk_lo,input$fu_lo)
      #*input$fu_lo
      }
    else if (input$pkmethod_lo == "ivive_lm") {
      lm_ivive(input$micpk_lo,input$fu_lo)
      #*input$fu_lo
      }
    else if (input$pkmethod_lo == "alloscale") {
      res <- model_Cl_results_lo()
      res$clin_CL
      }
    else {
      NA  # in case neither method is selected
      }
    })
  
  clin_V_lo <- reactive({
    bw_val <- BW()
    input$preclin_vss1_lo*((55/bw_val)^1)
    })
  
  output$clin_V_lo <- renderText({
    val <- clin_V_lo()
    res <- model_V_results_lo()
    paste0(
      "Estimated Human Vss: ", round(val, 2)," L\n",
      "Estimated Human Vss: ", round(res$clin_V, 2), " L")
    })

  ## Plasma PK --------
  sim_result <- eventReactive(input$click_lo, {
    res <- clin_CL_lo()
    res2 <- model_V_results_lo()
    cmt_choice <- ifelse(input$pk_route_lo == "Oral", "GUT", "CP")
    
    run_pk_sim(
      cmt_choice, 
      input$ka_lo,
      res, 
      res2$clin_V, 
      input$dose_lo, 
      input$ndoses_lo, 
      input$inter_lo)
    })
  
  ### Plot --------
  plot3_data <- reactiveVal(NULL)
  
  observeEvent(input$click_lo, {
    res <- sim_result()
    plot3_data(plot_pk(res, input$MIC, NA, NA, NA))
    })
  
  output$plot3_lo <- renderPlot({
    validate(
      need(
        !is.na(input$dose_lo) && !is.na(input$ndoses_lo) && !is.na(input$inter_lo),
        "Please enter dose information."
        )
      )
    validate(
      need(
        !is.na(clin_CL_lo()),
        "Please enter PK parameters."
        )
      )
    plot3_data()
    })
  
  ## Plasma Coverage--------
  coverage_result <- eventReactive(input$click_lo, {

    res2 <- model_V_results_lo()
    mic_unbound_adj <- input$MIC / input$fu_lo
    
    run_coverage_sim(
        clin_CL  = clin_CL_lo(),
        clin_V   = res2$clin_V,
        dose_vec = c(0.1, 1, 5, 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000),
        ndoses   = input$ndoses_lo,
        inter    = input$inter_lo,
        MIC      = mic_unbound_adj,
        macIC90  = NA,
        casMBC90 = NA,
        ec50     = NA
      )
    
    })
  
  ### Plot --------
  plot6_data <- reactiveVal(NULL)
  
  observeEvent(input$click_lo, {
    coverage_df <- coverage_result()
    plot6_data(plot_coverage(coverage_df, input$MIC, NA, NA, NA))
    })
  
  output$plot6_lo <- renderPlot({
    validate(
      need(
        !is.na(input$dose_lo) && !is.na(input$ndoses_lo) && !is.na(input$inter_lo),
        "Please enter dose information."
        )
      )
    validate(
      need(
        !is.na(clin_CL_lo()),
        "Please enter PK parameters."
        )
      )
    validate(
      need(
        !is.na(input$MIC),
        "Please enter MIC."
        )
      )
    plot6_data()
  })
  
  
  ## Outputs ----
  
  ### Parameters ----
  param_table_obj_lo <- reactive({
    
    # 1. Human CL
    CL_h <- clin_CL_lo()
    
    # 2. Human unbound CL
    CL_h_u <- CL_h / input$fu_lo
    
    # 3. Human plasma Css 
    Css <- input$dose_lo/CL_h
    
    # 4. Human plasma Css
    Css_u <- input$dose_lo/CL_h_u
    
    df <- data.frame(
      Parameter = c("Predicted human total CL (L/h)", 
                    "Predicted human unbound CL (L/h)", 
                    "Predicted human total plasma Css (mg/L)",
                    "Predicted human unbound plasma Css (mg/L)"),
      Estimate = as.character(c(round(CL_h, 2),
                                round(CL_h_u, 2),
                                round(Css,2),
                                round(Css_u,2)))
      )
    
    colnames(df)[1] <- "Parameter"
    colnames(df)[2] <- "Value"
    
    df
    }
    )
  
  output$param_table_lo <- DT::renderDataTable({
    param_table_obj_lo()
    }, 
    options = list(
      paging = FALSE,
      searching = FALSE,
      info = FALSE,
      dom = 't', # show only table
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
    ),
    rownames = FALSE,
    class = 'cell-border stripe hover compact'
    )
  
  ### Dose ----
  result_table_obj_lo <- reactive({
    coverage_df <- coverage_result()
    plasma_MIC = NA
    
    if (!is.na(input$MIC)){
      plasma_MIC <- calc_ec50_interp(coverage_df, target = 50, coverage_col = "perc_TAM")
      }
    
    df <- data.frame(
      Parameter = c("Unbound plasma PK/MIC"
                    ),
      Estimate = as.character(c(plasma_MIC
                                ))
      )
    
    colnames(df)[1] <- "PK-PD Target"
    colnames(df)[2] <- "Projected ED50 (mg)"
    
    df
    }
    )
  
  output$result_table_lo <- DT::renderDataTable({
    result_table_obj_lo()
    }, 
    options = list(
      paging = FALSE,
      searching = FALSE,
      info = FALSE,
      dom = 't', # show only table
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
    ),
    rownames = FALSE,
    class = 'cell-border stripe hover compact'
  )
  
  
  # -------- SP --------
  
  ## Reset ----
  observeEvent(input$reset_sp, {
    # --- Compound Info ---
    updateTextInput(session,"drugname_sp", value = NA)
    updateSelectInput(session, "drugbackbone_sp", selected = " ")
    
    # --- Simulation Params ---
    updateNumericInput(session, "ka_sp", value = NA)
    updateNumericInput(session, "dose_sp", value = NA)
    updateNumericInput(session, "ndoses_sp", value = NA)
    updateNumericInput(session, "inter_sp", value = NA)
    
    updateNumericInput(session, "fu_sp", value = NA)
    updateNumericInput(session, "heppk_sp", value = NA)
    updateNumericInput(session, "micpk_sp", value = NA)
    
    updateNumericInput(session, "preclin_cl1_sp", value = NA)
    updateNumericInput(session, "preclin_vss1_sp", value = NA)
    updateNumericInput(session, "preclin_cl2_sp", value = NA)
    updateNumericInput(session, "preclin_vss2_sp", value = NA)
    updateNumericInput(session, "preclin_cl3_sp", value = NA)
    updateNumericInput(session, "preclin_vss3_sp", value = NA)
    
    # --- PD Parameters ---
    updateNumericInput(session, "PC_caseum", value = NA)
    updateNumericInput(session, "casMBC90_sp", value = NA)
    
    plot4_data(NULL)
    plot5_data(NULL)
    plot7_data(NULL)
    
    showNotification("Parameters reset.", 
                     type = "default",
                     id = "reset_sp")
  })
  
  
  ## Sample params ----
  observeEvent(input$sample_params_sp, {
    # --- Compound Info ---
    updateTextInput(session, "drugname_sp", value = "Demo Compound")
    updateSelectInput(session, "drugbackbone_sp", selected = " ")
    
    # --- Simulation Parameters ---
    updateNumericInput(session, "ka_sp", value = 0.5)
    updateNumericInput(session, "dose_sp", value = 200)
    updateNumericInput(session, "ndoses_sp", value = 14)
    updateNumericInput(session, "inter_sp", value = 24)
    
    # --- PK Parameters ---
    updateNumericInput(session, "fu_sp", value = 0.8)
    updateNumericInput(session, "heppk_sp", value = 15.5)
    updateNumericInput(session, "micpk_sp", value = 10.2)
    
    # --- Preclinical PK ---
    updateNumericInput(session, "preclin_cl1_sp", value = 20)
    updateNumericInput(session, "preclin_vss1_sp", value = 45)
    updateNumericInput(session, "preclin_cl2_sp", value = 40)
    updateNumericInput(session, "preclin_vss2_sp", value = 30)
    updateNumericInput(session, "preclin_cl3_sp", value = 50)
    updateNumericInput(session, "preclin_vss3_sp", value = 56)
    
    # --- PD Parameters ---
    updateNumericInput(session, "PC_caseum", value = 0.5)
    updateNumericInput(session, "casMBC90_sp", value = 0.005)
    
    # Send a little toast message to let the user know it worked
    showNotification("Sample parameters for Demo Compound loaded.", 
                     type = "default",
                     id = "reset_sp")
  })
  
  
  ## Reactivity ----
  inputs_sp <- reactive({
    list(input$pk_route_sp,input$ka_sp,
         input$dose_sp, input$ndoses_sp, input$inter_sp, 
         input$pkmethod_sp,
         input$heppk_sp, input$fu_sp, input$micpk_sp,
         input$preclin_cl1_sp, input$preclin_cl2_sp, input$preclin_cl3_sp,
         input$preclin_vss1_sp, input$preclin_vss2_sp, input$preclin_vss3_sp,
         input$PC_caseum,
         input$casMBC90_sp
         )
    })
  
  # When any input changes, gray out the plot
  observeEvent(inputs_sp(), {
    shinyjs::addClass(id = "plot4_sp", class = "out-of-sync")
    shinyjs::addClass(id = "plot7_sp", class = "out-of-sync")
    shinyjs::addClass(id = "param_table_sp", class = "out-of-sync")
    shinyjs::addClass(id = "result_table_sp", class = "out-of-sync")
    
    if (input$tabs == "pred_sp") {
      showNotification("Input changed. Please update plots.", 
                       type = "default",
                       id = "input_sp")
    }
    })
  
  # When the button is clicked, remove gray effect
  observeEvent(input$click_sp, {
    shinyjs::removeClass(id = "plot4_sp", class = "out-of-sync")
    shinyjs::removeClass(id = "plot7_sp", class = "out-of-sync")
    shinyjs::removeClass(id = "param_table_sp", class = "out-of-sync")
    shinyjs::removeClass(id = "result_table_sp", class = "out-of-sync")
  })
  
  observeEvent(input$click_sp, {
    # Calculate Css,avg for the lesion
    # Css = (AUC_lesion / tau)
    # Css_lesion = (Dose / CL) * PC / tau
    
    auc_lesion <- (input$dose_sp / clin_CL_sp()) * input$PC_caseum
    css_lesion <- auc_lesion / input$inter_sp
    ratio <- css_lesion / input$casMBC90_sp
    
    print(ratio)
    
    if (ratio < 1) {
      showNotification(
        "WARNING: Compound is not highly penetrating (caseum Css/casMBC90 <1). 
        Caseum PK/casMBC90 may be unreliable; please continue to murine studies",
        type = "warning",
        duration = 10
        )
      }
    })
  
  observe({
    # 1. Define required inputs and their user-friendly names
    reqs <- list(
      "Dose" = input$dose_sp,
      "Number of doses" = input$ndoses_sp,
      "Dosing interval" = input$inter_sp,
      "Plasma-to-caseum PC" = input$PC_caseum,
      "casMBC90" = input$casMBC90_sp
      )
    
    # 2. Check for NA or zero values
    missing <- names(reqs)[sapply(reqs, function(x) is.null(x) || is.na(x) || x <= 0)]
    
    if (length(missing) > 0) {
      
      shinyjs::disable("click_sp") # Disable button
      
      # Create message
      output$sp_validation_ui <- renderUI({
        div(style = "background-color: #fcf8e3; border: 1px solid #faebcc; padding: 10px; border-radius: 4px; margin-bottom: 10px;",
            tags$b(style = "color: #8a6d3b;", icon("exclamation-triangle"), " Please enter:"),
            tags$ul(style = "color: #8a6d3b; margin-top: 5px; margin-bottom: 0;",
                    lapply(missing, function(m) tags$li(m))
            )
        )
      })
      shinyjs::addClass("click_pc", "btn-outline-secondary")
      shinyjs::removeClass("click_pc", "btn-primary")
    } 
    
    else {
      # Enable button
      shinyjs::enable("click_sp")
      output$sp_validation_ui <- renderUI({
        div(style = "background-color: #dff0d8; border: 1px solid #d6e9c6; padding: 10px; border-radius: 4px; margin-bottom: 10px;",
            tags$span(style = "color: #3c763d;", icon("check-circle"), " Ready to simulate.")
        )
      })
      
      # Clear tooltip and restore primary color
      shinyjs::addClass("click_sp", "btn-primary")
      shinyjs::removeClass("click_sp", "btn-outline-secondary")
    }
  })
  
  
  ## --- Clearance ----
  ### ---- Single species allometric scaling
  BW <- reactive({
    species_params %>% filter(Species == input$species1_sp) %>% pull(BW)
    })
  
  ### ---- Multi-species allometric scaling
  df_sp <- reactive({
    species <- c(input$species1_sp, input$species2_sp, input$species3_sp)
    cls <- c(input$preclin_cl1_sp, input$preclin_cl2_sp, input$preclin_cl3_sp)
    vs <- c(input$preclin_vss1_sp, input$preclin_vss2_sp, input$preclin_vss3_sp)
    
    weights_df <- species_params %>% filter(Species %in% species)
    weights <- weights_df$BW[match(species, weights_df$Species)]
    
    df_temp <- data.frame(CL = cls, V = vs, Species = species, BodyWeight = weights)
    df_temp$logBW <- log10(df_temp$BodyWeight)
    df_temp$logCL <- log10(df_temp$CL)
    df_temp$logV <- log10(df_temp$V)
    
    df_temp
  })
  
  model_Cl_results_sp <- reactive({
    model_cl_fit(df_sp())
    })
  
  model_V_results_sp <- reactive({
    model_v_fit(df_sp())
    })
  
  plot_sp_obj <- reactive({
    res <- model_Cl_results_sp()
    plot_log_fit(df_sp(), res, "logCL", "log10 Clearance (L/h)")
    })
  
  output$plot_sp <- renderPlot({
    plot_sp_obj()
    })
  
  plot2_sp_obj <- reactive({
    res <- model_V_results_sp()
    plot_log_fit(df_sp(), res, "logV", "log10 Volume (L)")
    })
  
  output$plot2_sp <- renderPlot({
    plot2_sp_obj()
    })
  
  clin_CL_sp <- reactive({
    if (input$pkmethod_sp == "ivive_h") {
      hep_ivive(input$heppk_sp,input$fu_sp)
      } 
    else if (input$pkmethod_sp == "ivive_lm") {
      lm_ivive(input$micpk_sp,input$fu_sp)
      }
    else if (input$pkmethod_sp == "alloscale") {
      res <- model_Cl_results_sp()
      res$clin_CL
      }
    else {
      NA
    }
  })
  
  clin_V_sp <- reactive({
    bw_val <- BW()
    input$preclin_vss1_sp*((55/bw_val)^1)
    })
  
  ### ---- Print clinical CL ----
  output$clin_CL_sp <- renderText({
    val <- clin_CL_sp()
    paste0(
      "Estimated Human Plasma CL: ", round(val, 2), " L/hr\n"
      )
    })
  
  ## Caseum PK --------
  sim_caseum_result <- eventReactive(input$click_sp, {
    res <- clin_CL_sp()
    res2 <- model_V_results_sp()
    cmt_choice <- ifelse(input$pk_route_sp == "Oral", "GUT", "CP")
    
    run_lesion_pk_sim(
      cmt_choice,
      input$ka_sp,
      res, 
      res2$clin_V, 
      input$dose_sp, 
      input$ndoses_sp, 
      input$inter_sp, 
      input$PC_caseum, 
      KPL = 10
      )
    })
  
  ### Plot  --------
  plot4_data <- reactiveVal(NULL)
  
  observeEvent(input$click_sp, {
    res <- sim_caseum_result()
    plot4_data(plot_lesion_pk(res, NA, NA, input$casMBC90_sp, NA))
    })
  
  output$plot4_sp <- renderPlot({
    validate(
      need(
        !is.na(input$dose_sp) && !is.na(input$ndoses_sp) && !is.na(input$inter_sp),
        "Please enter dose information."
        )
      )
    validate(
      need(
        !is.na(clin_CL_sp()) && !is.na(clin_V_sp()),
        "Please enter PK parameters."
      )
    )
    validate(
      need(
        !is.na(input$PC_caseum),
        "Please enter plasma-to-caseum partition coefficient."
        )
      )
    plot4_data()
  })
  
  ## Caseum Coverage --------
  coverage_result_caseum <- eventReactive(input$click_sp, {
    res2 <- model_V_results_sp()
    
    run_coverage_sim_lesion(
      clin_CL  = clin_CL_sp(),
      clin_V   = res2$clin_V,
      dose_vec = c(0.1, 1, 5, 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000),
      ndoses   = input$ndoses_sp,
      inter    = input$inter_sp,
      PC       = input$PC_caseum,
      KPL      = 10,
      MIC      = NA,
      macIC90  = NA,
      casMBC90 = input$casMBC90_sp,
      ec50     = NA
      )
    })
  
  ### Plot  --------
  plot7_data <- reactiveVal(NULL)
  
  observeEvent(input$click_sp, {
    coverage_caseum_df <- coverage_result_caseum()
    plot7_data(plot_coverage(coverage_caseum_df, NA, NA, input$casMBC90_sp, NA))
    })
  
  output$plot7_sp <- renderPlot({
    validate(
      need(
        !is.na(input$dose_sp) && !is.na(input$ndoses_sp) && !is.na(input$inter_sp),
        "Please enter dose information."
        )
      )
    validate(
      need(
        !is.na(clin_CL_sp()) && !is.na(clin_V_sp()),
        "Please enter PK parameters."
        )
      )
    validate(
      need(
        !is.na(input$PC_caseum),
        "Please enter plasma-to-caseum partition coefficient."
      )
    )
    validate(
      need(
        !is.na(input$casMBC90_sp),
        "Please enter caseum MBC90."
        )
      )
    plot7_data()
  })
  
  
  ## Outputs ----
  
  ### Parameters ----
  param_table_obj_sp <- reactive({
    
    # 1. Human CL 
    CL_h <- clin_CL_sp()
    
    # 2. Human unbound CL
    CL_h_u <- CL_h / input$fu_sp
    
    # 3. Human caseum Css 
    h_tau <- input$inter_sp
    Css <- input$dose_sp/(CL_h*h_tau)*input$PC_caseum
    
    # 4. Human unbound caseum Css
    Css_u <- input$dose_sp/(CL_h_u*h_tau)*input$PC_caseum
    
    df <- data.frame(
      Parameter = c("Predicted human total CL (L/h)", 
                    "Predicted human unbound CL (L/h)", 
                    "Predicted human total caseum Css (mg/L)",
                    "Predicted human unbound caseum Css (mg/L)"),
      Estimate = as.character(c(round(CL_h, 2), 
                                round(CL_h_u, 2),
                                round(Css, 4),
                                round(Css_u, 4)))
      )
    
    colnames(df)[1] <- "Parameter"
    colnames(df)[2] <- "Value"
    df
    }
  )
  
  output$param_table_sp <- DT::renderDataTable({
    param_table_obj_sp()
    }, 
    options = list(
      paging = FALSE,
      searching = FALSE,
      info = FALSE,
      dom = 't', # show only table
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
    ),
    rownames = FALSE,
    class = 'cell-border stripe hover compact'
    )
  
  ### Dose ----
  result_table_obj_sp <- reactive({
    
    coverage_lesion_df <- coverage_result_caseum()
    
    lesion_mbc90 = NA
    
    if (!is.na(input$casMBC90_sp)){
      lesion_mbc90 <- calc_ec50_interp(coverage_lesion_df, target = 50, coverage_col = "perc_TACAS")
    }
    
    df <- data.frame(
      Parameter = c( 
                    "Caseum PK/casMBC90"),
      Estimate = as.character(c(
                                lesion_mbc90))
      )
    colnames(df)[1] <- "PK-PD Target"
    colnames(df)[2] <- "Projected ED50 (mg)"
    df
    }
  )
  
  output$result_table_sp <- DT::renderDataTable({
    result_table_obj_sp()
      }, 
    options = list(
      paging = FALSE,
      searching = FALSE,
      info = FALSE,
      dom = 't', # show only table
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
    rownames = FALSE,
    class = 'cell-border stripe hover compact'
    )
  
  
  
  # -------- PC --------
  
  ## Reset ----
  observeEvent(input$reset_pc, {
    # --- Compound Info ---
    updateTextInput(session,"drugname_pc", value = NA)
    updateSelectInput(session, "drugbackbone_pc", selected = " ")
    
    # --- Simulation Params ---
    updateNumericInput(session, "ka_pc", value = NA)
    updateNumericInput(session, "dose_m_pc", value = NA)
    updateNumericInput(session, "dose_pc", value = NA)
    updateNumericInput(session, "ndoses_pc", value = NA)
    updateNumericInput(session, "inter_pc", value = NA)
    
    updateNumericInput(session, "fu_pc", value = NA)
    updateNumericInput(session, "heppk_pc", value = NA)
    updateNumericInput(session, "micpk_pc", value = NA)
    
    updateNumericInput(session, "preclin_cl1_pc", value = NA)
    updateNumericInput(session, "preclin_vss1_pc", value = NA)
    updateNumericInput(session, "preclin_cl2_pc", value = NA)
    updateNumericInput(session, "preclin_vss2_pc", value = NA)
    updateNumericInput(session, "preclin_cl3_pc", value = NA)
    updateNumericInput(session, "preclin_vss3_pc", value = NA)
    
    # --- PD Parameters ---
    updateNumericInput(session, "PC_caseum_pc", value = NA)
    updateNumericInput(session, "PC_cell", value = NA)
    updateNumericInput(session, "MIC_pc", value = NA)
    updateNumericInput(session, "casMBC90_pc", value = NA)
    updateNumericInput(session, "ec50", value = NA)
    updateNumericInput(session, "Emax", value = NA)
    updateNumericInput(session, "hill_co", value = NA)
    
    
    plot4_pc_data(NULL)
    plot5_data(NULL)
    plot7_pc_data(NULL)
    plot8_data(NULL)

    showNotification("Parameters reset.", 
                     type = "default",
                     id = "reset_pc")
  })
  
  ## Sample params ----
  observeEvent(input$sample_params_pc, {
    # --- Compound Info ---
    updateTextInput(session, "drugname_pc", value = "Demo Compound")
    updateSelectInput(session, "drugbackbone_pc", selected = " ")
    
    # --- Simulation Parameters ---
    updateNumericInput(session, "ka_pc", value = 0.5)
    updateNumericInput(session, "dose_m_pc", value = 20)
    updateNumericInput(session, "dose_pc", value = 200)
    updateNumericInput(session, "ndoses_pc", value = 14)
    updateNumericInput(session, "inter_pc", value = 24)
    
    # --- PK Parameters ---
    updateNumericInput(session, "fu_pc", value = 0.25)
    updateNumericInput(session, "heppk_pc", value = 15.5)
    updateNumericInput(session, "micpk_pc", value = 10.2)
    
    # --- Preclinical PK ---
    updateNumericInput(session, "preclin_cl1_pc", value = 20)
    updateNumericInput(session, "preclin_vss1_pc", value = 45)
    updateNumericInput(session, "preclin_cl2_pc", value = 40)
    updateNumericInput(session, "preclin_vss2_pc", value = 30)
    updateNumericInput(session, "preclin_cl3_pc", value = 50)
    updateNumericInput(session, "preclin_vss3_pc", value = 56)
    
    # --- PD Parameters ---
    updateNumericInput(session, "PC_caseum_pc", value = 0.5)
    updateNumericInput(session, "PC_cell", value = 0.8)
    updateNumericInput(session, "casMBC90_pc", value = 0.005)
    updateNumericInput(session, "ec50", value = 4)
    updateNumericInput(session, "Emax", value = 0.2)
    updateNumericInput(session, "hill_co", value = 1)
    
    # Send a toast message to let the user know it worked
    showNotification("Sample parameters for Demo Compound loaded.", 
                     type = "default",
                     id = "reset_pc")
    })
  
  
  ## Reactivity ----
  inputs_pc <- reactive({
    list(input$drugname_pc ,input$drugbackbone_pc,
         input$pk_route_pc, input$ka_pc,
         input$dose_m_pc, input$dose_pc, input$ndoses_pc, input$inter_pc, 
         input$pkmethod_pc, 
         input$heppk_pc, input$fu_pc, input$micpk_pc,
         input$preclin_cl1_pc, input$preclin_cl2_pc, input$preclin_cl3_pc,
         input$preclin_vss1_pc, input$preclin_vss2_pc, input$preclin_vss3_pc,
         input$PC_caseum_pc, input$PC_cell,
         input$casMBC90_pc, 
         input$ec50, input$Emax, input$hill_co
         )
    })
  
  # When any input changes, gray out the plots and send notification
  observeEvent(inputs_pc(), {
    shinyjs::addClass(id = "plot4_pc", class = "out-of-sync")
    shinyjs::addClass(id = "plot5", class = "out-of-sync")
    shinyjs::addClass(id = "plot7_pc", class = "out-of-sync")
    shinyjs::addClass(id = "plot8", class = "out-of-sync")
    shinyjs::addClass(id = "param_table_pc", class = "out-of-sync")
    shinyjs::addClass(id = "result_table_pc", class = "out-of-sync")
    
    if (input$tabs == "pred_pc") {
      showNotification("Input changed. Please update plots.", 
                     type = "default",
                     id = "warn_pc")
      }
  })
  
  # When the button is clicked, remove gray effect
  observeEvent(input$click_pc, {
    shinyjs::removeClass(id = "plot4_pc", class = "out-of-sync")
    shinyjs::removeClass(id = "plot5", class = "out-of-sync")
    shinyjs::removeClass(id = "plot7_pc", class = "out-of-sync")
    shinyjs::removeClass(id = "plot8", class = "out-of-sync")
    shinyjs::removeClass(id = "param_table_pc", class = "out-of-sync")
    shinyjs::removeClass(id = "result_table_pc", class = "out-of-sync")
  })
  
  ### Recommendation ----
  observeEvent(input$click_pc, {
    
    # 5. Mouse unbound CL
    m_CL_u <- input$preclin_cl1_pc * 0.02 / input$fu_pc
    
    # 6. Mouse unbound plasma Css
    m_dose <- input$dose_m_pc * 0.02
    m_tau <- input$inter_pc
    css_p <- m_dose/(m_CL_u*m_tau)

    ratio_1 <- css_p / input$ec50
    
    auc_lesion <- (input$dose_pc / clin_CL_pc()) * input$PC_caseum_pc
    css_lesion <- auc_lesion / input$inter_pc
    ratio_2 <- css_lesion / input$casMBC90_pc
    
    if (ratio_1 < 3) {
      showNotification(
        HTML("Compound is not highly potent, recommend using <b>cellular lesion PK/EC50</b> to project dose"),
        type = "message",
        id = "pc_msg",
        duration = 10
        )
      }
    
    if (ratio_1 > 3) {
      if (ratio_2 < 1){
        showNotification(
          HTML("Compound is highly potent but not highly penetrating, recommend using <b>cellular lesion PK/EC90</b> to project dose."),
          type = "message",
          id = "pc_msg",
          duration = 10
          )
        }
      if (ratio_2 >= 1) {
        showNotification(
          HTML("Compound is highly potent and highly penetrating, recommend using <b>caseum PK/casMBC90</b> to project dose."),
          type = "message",
          id = "pc_msg",
          duration = 10
          )
      }
      }
    })
  
  observe({
    # 1. Define required inputs and their user-friendly names
    reqs <- list(
      "Mouse Dose" = input$dose_m_pc,
      "Human Dose" = input$dose_pc,
      "Number of doses" = input$ndoses_pc,
      "Dosing interval" = input$inter_pc,
      "Plasma-to-caseum PC" = input$PC_caseum_pc,
      "Plasma-to-cell PC" = input$PC_cell,
      "casMBC90" = input$casMBC90_pc,
      "EC50" = input$ec50,
      "Emax" = input$Emax,
      "Hill coefficient" = input$hill_co
      )
    
    # 2. Check for NA or zero values
    missing <- names(reqs)[sapply(reqs, function(x) is.null(x) || is.na(x) || x <= 0)]
    
    if (length(missing) > 0) {
      # Disable button
      shinyjs::disable("click_pc")
      
      # Create message
      output$pc_validation_ui <- renderUI({
        div(style = "background-color: #fcf8e3; border: 1px solid #faebcc; padding: 10px; border-radius: 4px; margin-bottom: 10px;",
            tags$b(style = "color: #8a6d3b;", icon("exclamation-triangle"), " Please enter:"),
            tags$ul(style = "color: #8a6d3b; margin-top: 5px; margin-bottom: 0;",
                    lapply(missing, function(m) tags$li(m))
            )
        )
      })
      shinyjs::addClass("click_pc", "btn-outline-secondary")
      shinyjs::removeClass("click_pc", "btn-primary")
    } 
    
    else {
      # Enable button
      shinyjs::enable("click_pc")
      output$pc_validation_ui <- renderUI({
        div(style = "background-color: #dff0d8; border: 1px solid #d6e9c6; padding: 10px; border-radius: 4px; margin-bottom: 10px;",
            tags$span(style = "color: #3c763d;", icon("check-circle"), " Ready to simulate.")
            )
        })
      
      # Clear tooltip and restore primary color
      shinyjs::addClass("click_pc", "btn-primary")
      shinyjs::removeClass("click_pc", "btn-outline-secondary")
    }
  })
  
  ## --- Clearance ----
  
  ### ---- Single species allometric scaling
  BW <- reactive({
    species_params %>% filter(Species == input$species1_pc) %>% pull(BW)
    })
  
  ### ---- Multi-species allometric scaling
  df <- reactive({
    species <- c(input$species1_pc, input$species2_pc, input$species3_pc)
    cls <- c(input$preclin_cl1_pc, input$preclin_cl2_pc, input$preclin_cl3_pc)
    vs <- c(input$preclin_vss1_pc, input$preclin_vss2_pc, input$preclin_vss3_pc)
    
    weights_df <- species_params %>% filter(Species %in% species)
    weights <- weights_df$BW[match(species, weights_df$Species)]
    
    df_temp <- data.frame(CL = cls, V = vs, Species = species, BodyWeight = weights)
    df_temp$logBW <- log10(df_temp$BodyWeight)
    df_temp$logCL <- log10(df_temp$CL)
    df_temp$logV <- log10(df_temp$V)
    
    df_temp
    })
  
  model_Cl_results <- reactive({
    model_cl_fit(df())
    })
  
  model_V_results <- reactive({
    model_v_fit(df())
    })
  
  plot_obj <- reactive({
    res <- model_Cl_results()
    plot_log_fit(df(), res, "logCL", "log10 Clearance (L/h)")
    })
  
  output$plot_pc <- renderPlot({
    plot_obj()
    })
  
  plot2_obj <- reactive({
    res <- model_V_results()
    plot_log_fit(df(), res, "logV", "log10 Volume (L)")
    })
  
  output$plot2_pc <- renderPlot({
    plot2_obj()
    })
  
  clin_CL_pc <- reactive({
    if (input$pkmethod_pc == "ivive_h") {
      hep_ivive(input$heppk_pc,input$fu_pc)
    }
    else if (input$pkmethod_pc == "ivive_lm") {
      lm_ivive(input$micpk_pc,input$fu_pc)
    }
    else if (input$pkmethod_pc == "alloscale") {
      res <- model_Cl_results()
      res$clin_CL
    }
    else {
      NA  # in case neither method is selected
    }
  })
  
  clin_V_pc <- reactive({
    bw_val <- BW()
    input$preclin_vss1_pc*((55/bw_val)^1)
    })
  
  output$clin_CL_pc <- renderText({
    val <- clin_CL_pc()
    paste0(
      "Estimated Human Plasma CL: ", round(val, 2), " L/hr\n"
      )
    })
  
  output$clin_V <- renderText({
    val <- clin_V_pc()
    res <- model_V_results()
    paste0(
      "Estimated Human Vss: ", round(val, 2)," L\n",
      "Estimated Human Vss: ", round(res$clin_V, 2), " L"
      )
    })
  
  
  ## Caseum PK ------
  sim_caseum_result_pc <- eventReactive(input$click_pc, {
    res <- clin_CL_pc()
    res2 <- model_V_results()
    cmt_choice <- ifelse(input$pk_route_pc == "Oral", "GUT", "CP")
    
    run_lesion_pk_sim(
      cmt_choice,
      input$ka_pc, 
      res, 
      res2$clin_V, 
      input$dose_pc, 
      input$ndoses_pc, 
      input$inter_pc, 
      input$PC_caseum_pc, 
      KPL = 10)
  })
  
  ### Plot  --------
  
  plot4_pc_data <- reactiveVal(NULL)
  
  observeEvent(input$click_pc, {
    res <- sim_caseum_result_pc()
    plot4_pc_data(plot_lesion_pk(res, NA, NA, input$casMBC90_pc, NA))
  })
  
  output$plot4_pc <- renderPlot({
    validate(
      need(
        !is.na(input$dose_pc) && !is.na(input$ndoses_pc) && !is.na(input$inter_pc),
        "Please enter dose information."
      )
    )
    validate(
      need(
        !is.na(clin_CL_pc()) && !is.na(clin_V_pc()),
        "Please enter PK parameters."
      )
    )
    validate(
      need(
        !is.na(input$PC_caseum_pc),
        "Please enter plasma-to-caseum partition coefficient."
      )
    )
    plot4_pc_data()
  })
  
  ## Cell lesion PK ------
  sim_lesion_result <- eventReactive(input$click_pc, {
    res <- clin_CL_pc()
    res2 <- model_V_results()
    cmt_choice <- ifelse(input$pk_route_pc == "Oral", "GUT", "CP")
    
    run_lesion_pk_sim(
      cmt_choice,
      input$ka_pc,
      res, 
      res2$clin_V, 
      input$dose_pc, 
      input$ndoses_pc, 
      input$inter_pc, 
      input$PC_cell, 
      KPL = 10
      )
    })
  
  ### Plot  --------
  plot5_data <- reactiveVal(NULL)
  
  observeEvent(input$click_pc, {
    res <- sim_lesion_result()
    plot5_data(plot_lesion_pk(res, NA, NA, NA, input$ec50))
    })
  
  output$plot5 <- renderPlot({
    validate(
      need(
        !is.na(input$dose_pc) && !is.na(input$ndoses_pc) && !is.na(input$inter_pc),
        "Please enter dose information."
        )
      )
    validate(
      need(
        !is.na(clin_CL_pc()) && !is.na(clin_V_pc()),
        "Please enter PK parameters."
        )
      )
    validate(
      need(
        !is.na(input$PC_cell),
        "Please enter plasma-to-lesion partition coefficient."
        )
      )
    plot5_data()
    })

  ### Coverage --------
  plot7_pc_data <- reactiveVal(NULL)
  
  observeEvent(input$click_pc, {
    coverage_caseum_df <- coverage_result_caseum_pc()
    coverage_lesion_df <- coverage_result_lesion()
    plot7_pc_data(plot_coverage(coverage_caseum_df, NA, NA, input$casMBC90_pc, NA))
    plot8_data(plot_coverage(coverage_lesion_df, NA, NA, NA, input$ec50))
    })
  
  ## ---- Caseum --------
  coverage_result_caseum_pc <- eventReactive(input$click_pc, {
    res2 <- model_V_results()
    run_coverage_sim_lesion(
      clin_CL  = clin_CL_pc(),
      clin_V   = res2$clin_V,
      dose_vec = c(0.1, 1, 5, 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000),
      ndoses   = input$ndoses_pc,
      inter    = input$inter_pc,
      PC       = input$PC_caseum_pc,
      KPL      = 10,
      MIC      = NA,
      macIC90  = NA,
      casMBC90 = input$casMBC90_pc,
      ec50     = NA
      )
    })
  
  ### Plot ------
  output$plot7_pc <- renderPlot({
    validate(
      need(
        !is.na(input$dose_pc) && !is.na(input$ndoses_pc) && !is.na(input$inter_pc),
        "Please enter dose information."
        )
      )
    validate(
      need(
        !is.na(clin_CL_pc()) && !is.na(clin_V_pc()),
        "Please enter PK parameters."
        )
      )
    validate(
      need(
        !is.na(input$PC_caseum_pc),
        "Please enter plasma-to-caseum partition coefficient."
        )
      )
    validate(
      need(
        !is.na(input$casMBC90_pc),
        "Please enter PD metric."
        )
      )
    plot7_pc_data()
  })
  
  
  ### Cell lesion ------
  coverage_result_lesion <- eventReactive(input$click_pc, {
    res2 <- model_V_results()
    run_coverage_sim_lesion(
      clin_CL  = clin_CL_pc(),
      clin_V   = res2$clin_V,
      dose_vec = c(0.1, 1, 5, 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000),
      ndoses   = input$ndoses_pc,
      inter    = input$inter_pc,
      PC       = input$PC_cell,
      KPL      = 10,
      MIC      = NA,
      macIC90  = NA,
      casMBC90 = NA,
      ec50     = input$ec50
      )
    })
  
  ### Plot --------
  plot8_data <- reactiveVal(NULL)
  
  output$plot8 <- renderPlot({
    validate(
      need(
        !is.na(input$dose_pc) && !is.na(input$ndoses_pc) && !is.na(input$inter_pc),
        "Please enter dose information."
        )
      )
    validate(
      need(
        !is.na(clin_CL_pc()) && !is.na(clin_V_pc()),
        "Please enter PK parameters."
        )
      )
    validate(
      need(
        !is.na(input$PC_cell),
        "Please enter plasma-to-lesion partition coefficient."
        )
      )
    validate(
      need(
        !is.na(input$ec50),
        "Please enter EC50."
        )
      )
    plot8_data()
    })
  
  
  ## Outputs ----
  
  ### Parameters ----
  param_table_obj_pc <- reactive({
    
    # 1. Human clearance
    CL_h <- clin_CL_pc()
    
    # 2. Human unbound CL
    CL_h_u <- CL_h / input$fu_pc
    
    # 3. Human caseum Css
    h_tau <- input$inter_pc
    Css_c <- (input$dose_pc / (CL_h*h_tau)) * input$PC_caseum_pc
    
    # 3. Human unbound caseum Css
    Css_c_u <- (input$dose_pc / (CL_h_u*h_tau)) * input$PC_caseum_pc
    
    # 5. Mouse unbound CL
    m_CL_u <- input$preclin_cl1_pc * 0.02 / input$fu_pc
    
    # 6. Mouse unbound plasma Css
    m_dose <- input$dose_m_pc * 0.02
    m_tau <- input$inter_pc
    Css_p <- m_dose/(m_CL_u*m_tau)
    
    df <- data.frame(
      "Parameter" = c("Predicted human total CL (L/h)", 
                    "Predicted human unbound CL (L/h)", 
                    "Predicted human total caseum Css (mg/L)", 
                    "Predicted human unbound caseum Css (mg/L)", 
                    "Mouse unbound CL (L/h)",
                    "Mouse unbound plasma Css (mg/L)"),
      "Value" = as.character(c(round(CL_h, 2), 
                                round(CL_h_u, 2), 
                                round(Css_c, 2), 
                                round(Css_c_u, 2), 
                                round(m_CL_u, 4),
                                round(Css_p, 4)))
      )
    df
    }
  )
  
  output$param_table_pc <- DT::renderDataTable({
    param_table_obj_pc()
    }, 
    options = list(
      paging = FALSE,
      searching = FALSE,
      info = FALSE,
      dom = 't', # show only table
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
    ),
    rownames = FALSE,
    class = 'cell-border stripe hover compact'
    )
  
  
  ### Dose projections ----
  
  result_table_obj <- reactive({
    coverage_df <- coverage_result_caseum_pc()
    coverage_lesion_df <- coverage_result_lesion()
    
    cas_mbc90 = NA
    lesion_ec50 = NA
    lesion_ec90 = NA
    
    
    if (!is.na(input$casMBC90_pc)){
      cas_mbc90 <- calc_ec50_interp(coverage_df, target = 50, coverage_col = "perc_TACAS")
      }
    
    if (!is.na(input$ec50)){
      lesion_ec50 <- calc_ec50_interp(coverage_lesion_df, target = 50, coverage_col = "perc_TAE")
      lesion_ec90 <- calc_ec50_interp(coverage_lesion_df, target = 50, coverage_col = "perc_TAEC90")
    }
    
    css <- (input$dose_m_pc / input$preclin_cl1_pc)
    ratio_1 <- css / input$ec50
    
    auc_lesion <- (input$dose_pc / clin_CL_pc()) * input$PC_caseum_pc
    css_lesion <- auc_lesion / input$inter_pc
    ratio_2 <- css_lesion / input$casMBC90_pc
    
    # Identify the recommended row name
    rec_row <- ""
    if (ratio_1 < 3) {
      rec_row <- "Cellular lesion PK/EC50"
    } else if (ratio_2 < 1) {
      rec_row <- "Cellular lesion PK/EC90"
    } else {
      rec_row <- "Caseum PK/casMBC90"
    }
    
    df <- data.frame(
      "PK-PD Target" = c("Caseum PK/casMBC90", 
                         "Cellular lesion PK/EC50", 
                         "Cellular lesion PK/EC90"),
      "Projected ED50 (mg)" = as.character(c(
                                             cas_mbc90, 
                                             lesion_ec50, 
                                             lesion_ec90)),
      "IsRecommended" = c(FALSE, FALSE, FALSE), # Placeholder
      check.names = FALSE
      )
    
    df$IsRecommended[df$`PK-PD Target` == rec_row] <- TRUE
    
    df
    }
  )
  
  output$result_table_pc <- DT::renderDataTable({
    DT::datatable(result_table_obj(), 
                  options = list(
                    paging = FALSE,
                    searching = FALSE,
                    info = FALSE,
                    dom = 't',
                    # Hide the 3rd column (IsRecommended) so the user doesn't see the TRUE/FALSE
                    columnDefs = list(
                      list(visible = FALSE, targets = 2), 
                      list(className = 'dt-center', targets = "_all")
                    )
                    ),
                  rownames = FALSE,
                  class = 'cell-border stripe hover compact') %>%
      DT::formatStyle(
        columns = c('PK-PD Target', 'Projected ED50 (mg)'), # Column to check
        valueColumns = 'IsRecommended', # Column that holds the logic
        backgroundColor = DT::styleEqual(TRUE, '#d4edda'), # Muted green background
      )
  })
  
  
  # ---- Download ----

  observe({
    # Identify which phase the user wants a report for
    phase <- input$report_phase
    
    # Define the required data objects for each phase
    # These must match the names of the reactives you use in params_list
    is_ready <- switch(phase,
                       "LO" = !is.null(plot3_data()) && !is.null(result_table_obj_lo()),
                       "SP" = !is.null(plot4_data()) && !is.null(result_table_obj_sp()),
                       "PC" = !is.null(plot5_data()) && !is.null(result_table_obj())
                       )
    
    if (isTruthy(is_ready)) {
      shinyjs::enable("downloadReport")
      shinyjs::addClass("downloadReport", "btn-primary")
      shinyjs::removeClass("downloadReport", "btn-outline-secondary")
      } 
    
    else {
      shinyjs::disable("downloadReport")
      shinyjs::removeClass("downloadReport", "btn-primary")
      shinyjs::addClass("downloadReport", "btn-outline-secondary")
      }
    })
  
  # 2. Control the Message (Display the Missing Data warning)
  output$download_status <- renderUI({
    phase <- input$report_phase
    
    is_ready <- switch(phase,
                       "LO" = !is.null(plot3_data()) && !is.null(result_table_obj_lo()),
                       "SP" = !is.null(plot4_data()) && !is.null(result_table_obj_sp()),
                       "PC" = !is.null(plot8_data()) && !is.null(result_table_obj())
    )
    
    if (!isTruthy(is_ready)) {
      div(style = "margin-top: 10px;",
          helpText(style = "color: #a94442; font-weight: bold;",
                   icon("info-circle"), 
                   paste("Please run the simulation to download."))
      )
    } else {
      NULL # Show nothing when data is ready
    }
  })
  
  
  
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste("HDP_report-", Sys.Date(), ".pdf", sep = "")
    },
    
    content = function(file) {
      withProgress(message = 'Generating PDF report...', value = 0, {
        
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      get_input_summary <- function(phase) {
        sfx <- paste0("_", tolower(phase))
        
        # Create a list of dataframes for each "Box"
        # This makes it easy to add/remove rows based on the phase
        
        # 1. Compound & Basic Sim
        common <- data.frame(
          Category = "Compound & Simulation",
          Parameter = c("Drug Name", 
                        "Dose (mg)", 
                        "Number of Doses", 
                        "Interval (h)"
                        ),
          Value = c(
            input[[paste0("drugname", sfx)]],
            input[[paste0("dose", sfx)]],
            input[[paste0("ndoses", sfx)]],
            input[[paste0("inter", sfx)]]
            )
          )
        
        
        plasma_pk <- data.frame(
          Category = "Plasma PK",
          Parameter = c("fu",
                        "Human HEP CLint (µL/min/10^6 cells)",
                        "Human MIC CLint (uL/min/mg)",
                        "Species 1", 
                        "Species 1 CL (L/h)", 
                        "Species 1 Vss (L)", 
                        "Species 2", 
                        "Species 2 CL (L/h)", 
                        "Species 2 Vss (L)",
                        "Species 3", 
                        "Species 3 CL (L/h)",
                        "Species 3 Vss (L)"
                        ),
          Value = c(
            input[[paste0("fu", sfx)]],
            input[[paste0("heppk", sfx)]],
            input[[paste0("micpk", sfx)]],
            input[[paste0("species1", sfx)]],
            input[[paste0("preclin_cl1", sfx)]],
            input[[paste0("preclin_vss1", sfx)]],
            input[[paste0("species2", sfx)]],
            input[[paste0("preclin_cl2", sfx)]],
            input[[paste0("preclin_vss2", sfx)]],
            input[[paste0("species3", sfx)]],
            input[[paste0("preclin_cl3", sfx)]],
            input[[paste0("preclin_vss3", sfx)]]
            )
          )
        
        # 2. Phase-Specific Metrics (The logic you requested)
        pd_metrics <- data.frame(Category = character(), Parameter = character(), Value = character())
        
        if (phase == "LO") {
          pd_metrics <- rbind(pd_metrics, data.frame(
            Category = "PD Targets", 
            Parameter = "MIC (mg/L)", 
            Value = as.character(input$MIC)
          ))
        } 
        
        else if (phase == "SP") {
          pd_metrics <- rbind(pd_metrics, data.frame(
            Category = "Site-of-Action PK", 
            Parameter = "Plasma-to-Caseum PC",
            Value = as.character(input$PC_caseum)
            ))
          
          pd_metrics <- rbind(pd_metrics, data.frame(
            Category = "PD Targets", 
            Parameter = "casMBC90 (mg/L)", 
            Value = as.character(input$casMBC90_sp)
          ))
        } 
        
        else if (phase == "PC") {
          
          pd_metrics <- rbind(pd_metrics, data.frame(
            Category = "Site-of-Action PK", 
            Parameter = c("Plasma-to-Caseum PC", 
                          "Plasma-to-Cell PC"),
            Value = as.character(c(input$PC_caseum_pc, 
                                   input$PC_cell))
            ))
          
          pd_metrics <- rbind(pd_metrics, data.frame(
            Category = "PD Targets", 
            Parameter = c("casMBC90 (mg/L)", 
                          "Infection Model", 
                          "EC50 (mg/L)", 
                          "Emax"),
            Value = as.character(c(input$casMBC90_pc, 
                                   input$infmod, 
                                   input$ec50, 
                                   input$Emax))
          ))
        }
        
        # Combine all sections
        final_df <- rbind(common, plasma_pk, pd_metrics)           # Combine into single df
        final_df$Category[duplicated(final_df$Category)] <- ""     # Remove repeated catagories
        
        return(final_df)
      }
      
      params_list <- switch(input$report_phase,
                           "LO" = list(
                             exposure_plot = plot3_data(), # Plasma
                             coverage_plot = plot6_data(), # Plasma
                             
                             input_table = get_input_summary("LO"),
                             param_table   = param_table_obj_lo(),
                             result_table  = result_table_obj_lo()
                             ),
                           
                           "SP" = list(
                             exposure_caseum = plot4_data(), # Caseum
                             coverage_caseum = plot7_data(), # Caseum
                             
                             input_table = get_input_summary("SP"),
                             param_table   = param_table_obj_sp(),
                             result_table  = result_table_obj_sp()
                             ),
                           
                           "PC" = list(
                             exposure_caseum  = plot4_pc_data(),
                             exposure_cellular = plot5_data(), 
                             coverage_caseum  = plot7_pc_data(), 
                             coverage_cellular = plot8_data(),
                             
                             input_table = get_input_summary("PC"),
                             param_table      = param_table_obj_pc(),
                             result_table     = result_table_obj()
                             )
                           )
      
      setProgress(0.5)                                                    # Move bar to 50%
      
      rmarkdown::render(
        input = tempReport,
        output_format = rmarkdown::pdf_document(), 
        output_file = file,
        params = params_list,
        envir = new.env(parent = globalenv())
        )
      
      setProgress(1)                                                     # Done
      }
    )
    }
  )
  
  
  
}     # End of server