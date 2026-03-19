## ---- Load Packages ----
library(shiny)
library(shinyjs)
library(bslib)
library(shinydashboard)
library(shinyWidgets)

animals <- c("BALB/c Mouse", "Rat", "Dog")

ui <- dashboardPage(
  skin = "blue", 
  dashboardHeader(title = "TB HDP"),
   
  dashboardSidebar(
    hr(),
    sidebarMenu(
      id = "tabs",
      menuItem("About", tabName = "about",icon = icon("graduation-cap")),
      menuItem("Predict", icon = icon("area-chart"),
               menuSubItem("Lead Optimization", tabName = "pred_lo", icon = icon("angle-right")),
               menuSubItem("Selection Phase", tabName = "pred_sp", icon = icon("angle-right")),
               menuSubItem("Preclinical Candidate", tabName = "pred_pc", icon = icon("angle-right"))
               ),
      menuItem("Report", tabName = "rep", icon = icon("area-chart")),
      menuItem("Acknowledgements", tabName = "Ack",icon = icon("handshake")),
      menuItem("Citation", tabName = "Cite", icon = icon("link"))
      )
    ),
  
  dashboardBody(
    shinyjs::useShinyjs(),
    tabItems(
      
      # About Section------
      tabItem(tabName = "about",
              box(width = 5, collapsible = FALSE, solidHeader = TRUE,
                  title = "About",
                  includeMarkdown("www/about.md")
              ),
              box(width = 7, 
                  br(),
                  div(style = "text-align: center;", 
                      tags$img(
                        src = "Decision_Tree.png",
                        style = "max-width: 100%; height: 55vh;"
                      )
                  ),
                  br())
              ),
      
      # 1. LO ------
      tabItem(tabName = "pred_lo",
              column(class = "param-container",
                     
                     fluidRow(
                       column(
                         6,
                         div(
                           style = "text-align:center;",
                           actionButton("reset_lo", "Reset", width = "100%")
                           )
                         ),
                       column(
                         6,
                         div(
                           style = "text-align:center;",
                           actionButton("sample_params_lo", "Sample parameters", width = "100%")
                           )
                         )
                       ),
                
                br(),
                
                ### Inputs ------
                box(title = "Compound", width = 12, collapsible = TRUE, collapsed = FALSE,
                    fluidRow(column(12,textInput("drugname_lo", "Name", "Demo Compound")))
                    ),
                
                box(title = "Simulation parameters", width = 12, collapsible = TRUE, collapsed = TRUE,
                    fluidRow(
                      column(12, 
                             shinyWidgets::radioGroupButtons(
                               inputId = "pk_route_lo",
                               label = "PK Route",
                               choices = c("IV", "Oral"),
                               selected = "IV",
                               justified = TRUE,
                               status = "primary"
                             )
                      )
                    ),
                    
                    # 2. Conditional input for ka (only shows if Oral is selected)
                    conditionalPanel(
                      condition = "input.pk_route_lo == 'Oral'",
                      fluidRow(
                        column(6),
                        column(6, 
                               numericInput("ka_lo", HTML("Absorption rate k<sub>a</sub> (1/h)"), 
                                            value = 0.5, min = 0, step = 0.1)
                               )
                        )
                      ),
                    
                    hr(), # Visual separator
                    
                    fluidRow(column(12, numericInput("dose_lo", "Dose (mg)", value = 200, min = 0))),
                    fluidRow(column(12, numericInput("ndoses_lo", "Number of doses", value = 14, min = 1, max = 20))),
                    fluidRow(column(12, div(style = "white-space: nowrap;", numericInput("inter_lo", "Dosing interval (h)", value = 24, min = 1))))
                ),
                
                box(title = HTML("<i>In vitro</i> PK"), width = 12, collapsible = TRUE, collapsed = TRUE,
                    fluidRow(column(6, numericInput("fu_lo", "Human fu (0-1)", value = 0.8, min = 0, max = 1, step = 0.01))),
                    fluidRow(column(12, div(style = "white-space: nowrap;",numericInput("heppk_lo", HTML("Human Hepatocyte CL<sub>int</sub> (µL/min/10<sup>6</sup> cells)"), value = 12, min = 0)))),
                    fluidRow(column(12, numericInput("micpk_lo", HTML("Human Microsomal CL<sub>int</sub> (µL/min/mg)"), value = 12, min = 0)))
                ),
                
                box(title = HTML("<i>In vivo</i> PK"), width = 12, collapsible = TRUE, collapsed = TRUE,
                    fluidRow(
                      column(4, selectInput("species1_lo", "Preclinical species", choices = animals,selected = "Mouse")),
                      column(4, numericInput("preclin_cl1_lo", "Plasma CL (L/hr/kg)", value = 20, min = 0)),
                      column(4, numericInput("preclin_vss1_lo", "Plasma Vss (L)", value = 45, min = 0))
                      ),
                    
                    fluidRow(
                      column(4, selectInput("species2_lo", "Preclinical species", choices = animals,selected = "Rat")),
                      column(4, numericInput("preclin_cl2_lo", "Plasma CL (L/hr/kg)", value = 40, min = 0)),
                      column(4, numericInput("preclin_vss2_lo", "Plasma Vss (L)", value = 30, min = 0))           
                      ),
                    
                    fluidRow(
                      column(4, selectInput("species3_lo", "Preclinical species", choices = animals,selected = "Dog")),
                      column(4, numericInput("preclin_cl3_lo", "Plasma CL (L/hr/kg)", value = 50, min = 0)),
                      column(4, numericInput("preclin_vss3_lo", "Plasma Vss (L)", value = 56, min = 0))            
                      )
                  ),
                
                box(title = HTML("<i>In vitro</i> PD"), width = 12, collapsible = TRUE, collapsed = FALSE,
                    fluidRow(column(12, numericInput("MIC", "MIC (mg/L)", value = 0.001, min = 0)))
                    ),
                
                column(12, radioButtons("pkmethod_lo", "Clearance Prediction Method",
                                       c("IVIVE - Hepatocyte" = "ivive_h", 
                                         "IVIVE - Liver microsome" = "ivive_lm", 
                                         "Allometry" = "alloscale")),
                       
                       conditionalPanel(condition = "input.pkmethod_lo == 'alloscale'",
                                        br(),
                                        fluidRow(column(12, plotOutput("plot_lo", height = 300))),
                                        fluidRow(column(12, plotOutput("plot2_lo", height = 300)))
                                        )
                       ),
                width = 4
                
              ),
              
              
              ### Plots ------
              column(4, class = "plot-container",
                     
                     uiOutput("lo_validation_ui"),
                     
                     div(
                       style = "text-align:center;",
                       actionButton("click_lo", "Predict Dose", width = "50%")
                       ),
                     
                     br(),
                     
                     tabsetPanel(
                       tabPanel("Plasma",
                                fluidRow(column(12, box(title = 'Plasma Exposure', width = 12, solidHeader = TRUE, status = "primary", class = "box-custom", plotOutput("plot3_lo", height = "25vh")))),
                                fluidRow(column(12, box(title = 'Plasma Coverage', width = 12, solidHeader = TRUE, status = "primary", class = "box-custom", plotOutput("plot6_lo", height = "25vh"))))
                                ),
                       )
                     ),
              
              ### Outputs -----
              column(4,
                     
                     #### Parameters -----
                     fluidRow(column(12,
                                     br(),
                                     br(),
                                     box(width = 12, DT::dataTableOutput("param_table_lo"))
                                     )
                              ),
                     #### Dose -----
                     fluidRow(column(12,
                                     box(width = 12, DT::dataTableOutput("result_table_lo"))
                                     )
                              )
                     )
              ), 
      
      
      ## 2. SP ------
      tabItem(tabName = "pred_sp",
              column(class = "param-container",
                     
                     fluidRow(
                       column(
                         6,
                         div(
                           style = "text-align:center;",
                           actionButton("reset_sp", "Reset", width = "100%")
                         )
                       ),
                       column(
                         6,
                         div(
                           style = "text-align:center;",
                           actionButton("sample_params_sp", "Sample parameters", width = "100%")
                         )
                       )
                     ),
                     
                     br(),
                     
                     ### Inputs ------
                     box(title = "Compound", width = 12, collapsible = TRUE, collapsed = FALSE,
                         fluidRow(column(12,textInput("drugname_sp", "Name", "Demo Compound")))
                         ),
                     
                     box(title = "Simulation parameters", width = 12, collapsible = TRUE, collapsed = TRUE,
                         
                         fluidRow(
                           column(12, 
                                  shinyWidgets::radioGroupButtons(
                                    inputId = "pk_route_sp",
                                    label = "PK Route",
                                    choices = c("IV", "Oral"),
                                    selected = "IV",
                                    justified = TRUE,
                                    status = "primary"
                                  )
                           )
                         ),
                         
                         # 2. Conditional input for ka (only shows if Oral is selected)
                         conditionalPanel(
                           condition = "input.pk_route_sp == 'Oral'",
                           fluidRow(
                             column(6),
                             column(6, 
                                    numericInput("ka_sp", HTML("Absorption rate k<sub>a</sub> (1/h)"), 
                                                 value = 0.5, min = 0, step = 0.1)
                             )
                           )
                         ),
                         
                         
                         fluidRow(column(12, numericInput("dose_sp", "Dose (mg)", value = 200, min = 0))),
                         fluidRow(column(12, numericInput("ndoses_sp", "Number of doses", value = 14, min = 1, max = 20))),
                         fluidRow(column(12, div(style = "white-space: nowrap;", numericInput("inter_sp", "Dosing interval (h)", value = 24, min = 1))))
                     ),
                     
                     box(title = HTML("<i>In vitro</i> PK"), width = 12, collapsible = TRUE, collapsed = TRUE,
                         fluidRow(column(6, numericInput("fu_sp", "Human fu (0-1)", value = 0.8, min = 0, max = 1, step = 0.01))),
                         fluidRow(column(12, div(style = "white-space: nowrap;",numericInput("heppk_sp", HTML("Human Hepatocyte CL<sub>int</sub> (µL/min/10<sup>6</sup> cells)"), value = 12, min = 0)))),
                         fluidRow(column(12, numericInput("micpk_sp", HTML("Human Microsomal CL<sub>int</sub> (µL/min/mg)"), value = 12, min = 0)))
                     ),
                     
                     box(title = HTML("<i>In vivo</i> plasma PK"), width = 12, collapsible = TRUE, collapsed = TRUE,
                         fluidRow(
                           column(4, selectInput("species1_sp", "Preclinical species", choices = animals,selected = "Mouse")),
                           column(4, numericInput("preclin_cl1_sp", "Plasma CL (L/hr/kg)", value = 20, min = 0)),
                           column(4, numericInput("preclin_vss1_sp", "Plasma Vss (L)", value = 45, min = 0))
                         ),
                         
                         fluidRow(
                           column(4, selectInput("species2_sp", "Preclinical species", choices = animals,selected = "Rat")),
                           column(4, numericInput("preclin_cl2_sp", "Plasma CL (L/hr/kg)", value = 40, min = 0)),
                           column(4, numericInput("preclin_vss2_sp", "Plasma Vss (L)", value = 30, min = 0))           
                         ),
                         
                         fluidRow(
                           column(4, selectInput("species3_sp", "Preclinical species", choices = animals,selected = "Dog")),
                           column(4, numericInput("preclin_cl3_sp", "Plasma CL (L/hr/kg)", value = 50, min = 0)),
                           column(4, numericInput("preclin_vss3_sp", "Plasma Vss (L)", value = 56, min = 0))            
                         )
                     ),
                     
                     box(title = HTML("<i>In vivo </i> site-of-action PK"), width = 12, collapsible = TRUE, collapsed = TRUE,
                         fluidRow(column(12, numericInput("PC_caseum", "Plasma-to-caseum partition coefficient", value = 0.5, min = 0, step = 0.01)))
                     ),
                     
                     box(title = HTML("<i>Ex vivo</i> PD"), width = 12, collapsible = TRUE, collapsed = FALSE,
                         fluidRow(column(12, numericInput("casMBC90_sp", HTML("casMBC<sub>90</sub> (mg/L)"), value = 0.005, min = 0)))
                     ),
                     
                     column(12, radioButtons("pkmethod_sp", "Clearance Prediction Method",
                                             c("IVIVE - Hepatocyte" = "ivive_h", 
                                               "IVIVE - Liver microsome" = "ivive_lm", 
                                               "Allometry" = "alloscale")),
                            
                            conditionalPanel(condition = "input.pkmethod_sp == 'alloscale'",
                                             br(),
                                             fluidRow(column(12, plotOutput("plot_sp", height = 300))),
                                             fluidRow(column(12, plotOutput("plot2_sp", height = 300)))
                                             )
                            ),
                     width = 4
                     ),
              
              
              ### Plots ------
              column(4, class = "plot-container",
                     
                     uiOutput("sp_validation_ui"),
                     
                     div(
                       style = "text-align:center;",
                       actionButton("click_sp", "Predict Dose", width = "50%")
                     ),
                     br(),
                     
                     tabsetPanel(
                       id = "main_tabs",           
                       selected = "Caseum", 
                      
                       tabPanel("Caseum", value = "Caseum",
                                fluidRow(column(12, box(title = 'Caseum Exposure', width = 12, solidHeader = TRUE, status = "primary", class = "box-custom", plotOutput("plot4_sp", height = "25vh")))),
                                fluidRow(column(12, box(title = 'Caseum Coverage', width = 12, solidHeader = TRUE, status = "primary", class = "box-custom", plotOutput("plot7_sp", height = "25vh"))))
                                )
                       )
                     ),
                     
                     ### Outputs -----
                     column(4,
                            
                            #### Parameters -----
                            fluidRow(column(12,
                                            br(),
                                            br(),
                                            box(width = 12, DT::dataTableOutput("param_table_sp"))
                                            )
                                     ),
                            #### Dose -----
                            fluidRow(column(12,
                                            box(width = 12, DT::dataTableOutput("result_table_sp"))
                                            )
                                     )
                            )
              ), 
      
      ## 3. PC ------
      tabItem(tabName = "pred_pc",
              
              column(class = "param-container",
                     
                     fluidRow(
                       column(6,
                         div(
                           style = "text-align:center;",
                           actionButton("reset_pc", "Reset", width = "100%")
                           )
                         ),
                       column(6,
                         div(
                           style = "text-align:center;",
                           actionButton("sample_params_pc", "Sample parameters", width = "100%")
                           )
                         )
                       ),
                     
                     br(),
                     
                     ### Inputs ------
                     box(title = "Compound", width = 12, collapsible = TRUE, collapsed = FALSE,
                         fluidRow(column(12,textInput("drugname_pc", "Name", "Demo Compound"))),
                         fluidRow(column(12,selectInput("drugbackbone_pc", "Backbone regimen", list(" ", "BPa", "PaL"))))
                     ),
                     
                     box(title = "Simulation parameters", width = 12, collapsible = TRUE, collapsed = TRUE,
                         
                         fluidRow(
                           column(12, 
                                  shinyWidgets::radioGroupButtons(
                                    inputId = "pk_route_pc",
                                    label = "PK Route",
                                    choices = c("IV", "Oral"),
                                    selected = "IV",
                                    justified = TRUE,
                                    status = "primary"
                                  )
                           )
                         ),
                         
                         # 2. Conditional input for ka (only shows if Oral is selected)
                         conditionalPanel(
                           condition = "input.pk_route_pc == 'Oral'",
                           fluidRow(
                             column(6),
                             column(6, 
                                    numericInput("ka_pc", HTML("Absorption rate k<sub>a</sub> (1/h)"), 
                                                 value = 0.5, min = 0, step = 0.1)
                             )
                           )
                         ),
                         
                         
                         fluidRow(column(12, numericInput("dose_m_pc", "Mouse dose (mg/kg)", value = 20, min = 0))),
                         fluidRow(column(12, numericInput("dose_pc", "Human dose (mg)", value = 200, min = 0))),
                         fluidRow(column(12, numericInput("ndoses_pc", "Number of doses", value = 14, min = 1, max = 20))),
                         fluidRow(column(12, div(style = "white-space: nowrap;", numericInput("inter_pc", "Dosing interval (h)", value = 24, min = 1))))
                     ),
                     
                     box(title = HTML("<i>In vitro</i> PK"), width = 12, collapsible = TRUE, collapsed = TRUE,
                         fluidRow(column(6, numericInput("fu_pc", "Human fu (0-1)", value = 0.8, min = 0, max = 1, step = 0.01))),
                         fluidRow(column(12, div(style = "white-space: nowrap;",numericInput("heppk_pc", HTML("Human Hepatocyte CL<sub>int</sub> (µL/min/10<sup>6</sup> cells)"), value = 12, min = 0)))),
                         fluidRow(column(12, numericInput("micpk_pc", HTML("Human Microsomal CL<sub>int</sub> (µL/min/mg)"), value = 12, min = 0)))
                     ),
                     
                     box(title = HTML("<i>In vivo</i> plasma PK"), width = 12, collapsible = TRUE, collapsed = TRUE,
                         fluidRow(
                           column(4, selectInput("species1_pc", "Preclinical species", choices = animals,selected = "Mouse")),
                           column(4, numericInput("preclin_cl1_pc", "Plasma CL (L/hr/kg)", value = 20, min = 0)),
                           column(4, numericInput("preclin_vss1_pc", "Plasma Vss (L)", value = 45, min = 0))
                         ),
                         
                         fluidRow(
                           column(4, selectInput("species2_pc", "Preclinical species", choices = animals,selected = "Rat")),
                           column(4, numericInput("preclin_cl2_pc", "Plasma CL (L/hr/kg)", value = 40, min = 0)),
                           column(4, numericInput("preclin_vss2_pc", "Plasma Vss (L)", value = 30, min = 0))           
                         ),
                         
                         fluidRow(
                           column(4, selectInput("species3_pc", "Preclinical species", choices = animals,selected = "Dog")),
                           column(4, numericInput("preclin_cl3_pc", "Plasma CL (L/hr)", value = 50, min = 0)),
                           column(4, numericInput("preclin_vss3_pc", "Plasma Vss (L)", value = 56, min = 0))            
                         )
                     ),
                     
                     box(title = HTML("<i>In vivo</i> site-of-action PK"), width = 12, collapsible = TRUE, collapsed = TRUE,
                         fluidRow(column(12, numericInput("PC_caseum_pc", "Plasma-to-caseum partition coefficient", value = 0.5, min = 0, step = 0.01))),
                         fluidRow(column(12, numericInput("PC_cell", "Plasma-to-cellular partition coefficient", value = 0.8, min = 0, step = 0.01)))
                     ),
                     
                     
                     box(title = HTML("<i>Ex vivo</i> PD"), width = 12, collapsible = TRUE, collapsed = FALSE,
                         fluidRow(column(12, numericInput("casMBC90_pc", HTML("casMBC<sub>90</sub> (mg/L)"), value = 0.005, min = 0)))
                     ),
                     
                     box(title = HTML("<i>In vivo</i> PD"), width = 12, collapsible = TRUE, collapsed = FALSE,
                         fluidRow(
                           column(6, selectInput("infmod", "Infection Model", list("Acute", "Subacute", "Chronic"))),
                           column(6, numericInput("ec50", HTML("EC<sub>50</sub> (mg/L)"), value = 4, min = 0))),
                         fluidRow(column(12,numericInput("Emax", HTML("E<sub>max</sub>"), value = 0.2, min = 0, max = 5))),
                         fluidRow(column(12,numericInput("hill_co", HTML("Hill Coefficient"), value = 1, min = 0.1, max = 5)))
                     ),
                     

                     column(12, radioButtons("pkmethod_pc", "Clearance Prediction Method",
                                             c("IVIVE - Hepatocyte" = "ivive_h", 
                                               "IVIVE - Liver microsome" = "ivive_lm", 
                                               "Allometry" = "alloscale")),
                            
                            conditionalPanel(condition = "input.pkmethod_pc == 'alloscale'",
                                             br(),
                                             fluidRow(column(12, plotOutput("plot_pc", height = 300))),
                                             fluidRow(column(12, plotOutput("plot2_pc", height = 300)))
                                             )
                            ),
                     width = 4
                     
              ),
              
              
              ### Plots ------
              column(4, class = "plot-container",
                     
                     uiOutput("pc_validation_ui"),
                     
                     div(
                       style = "text-align:center;",
                       actionButton("click_pc", "Predict Dose", width = "50%", class = "btn-primary")
                       ),
                     br(),
                     
                     
                     
                     tabsetPanel(
                       id = "main_tabs",           
                       selected = "Cellular Lesion", 
                       
                       tabPanel("Cellular Lesion", value = "Cellular Lesion",
                              fluidRow(column(12, box(title = 'Cellular Lesion Exposure', width = 12, solidHeader = TRUE, status = "primary", class = "box-custom", plotOutput("plot5", height = "25vh")))),
                              fluidRow(column(12, box(title = 'Cellular Lesion Coverage', width = 12, solidHeader = TRUE, status = "primary", class = "box-custom", plotOutput("plot8", height = "25vh"))))
                              ),
                       tabPanel("Caseum", value = "Caseum",
                                fluidRow(column(12, box(title = 'Caseum Exposure', width = 12, solidHeader = TRUE, status = "primary", class = "box-custom", plotOutput("plot4_pc", height = "25vh")))),
                                fluidRow(column(12, box(title = 'Caseum Coverage', width = 12, solidHeader = TRUE, status = "primary", class = "box-custom", plotOutput("plot7_pc", height = "25vh"))))
                              )
                     )),
              
              ### Outputs -----
              column(4,
                     
                     #### Parameters -----
                     fluidRow(column(12,
                                     br(),
                                     br(),
                                     box(width = 12, DT::dataTableOutput("param_table_pc"))
                                     )
                              ),
                     #### Dose -----
                     fluidRow(column(12,
                                     box(width = 12, DT::dataTableOutput("result_table_pc"))
                                     )
                              ),
                     )
              ), 

      
      
      ## Report Section ------          
      tabItem(tabName = "rep",
              selectInput("report_phase", "Select Phase to Export:",
                          choices = c("Lead Optimization" = "LO", 
                                      "Selection Phase" = "SP", 
                                      "Preclinical Candidate" = "PC")),
                uiOutput("download_status"),
                downloadButton("downloadReport", "Download Report")
              ),
      
      ## Acknowledgements Section------
      
      tabItem(tabName = "Ack",
              box(width = 12, title = "Acknowledgements", collapsible = FALSE, solidHeader = TRUE,
                  "[To be completed]",
                  br(),
                  br(),
                  "This work was supported by the Gates Foundation.",
                  br(),
                  div(style = "text-align: center;", tags$img(src = "logo.png", height = "150px"))
              )
      ),
      
      
      ## Citation Section------
      
      tabItem(tabName = "Cite",
              box( width = 12, collapsible = FALSE, solidHeader = TRUE, title = "Citation",  
                   "[To be completed]",
                   br(),
                   br(),
                   "Full methodology and results of this project was published in ...",
                   br(),
                   br(),
                   "Github code..."
              )
      )
    )
  ),

  
  tags$head(
    tags$style(HTML("
    
    .out-of-sync {
      opacity: 0.2; 
      filter: grayscale(100%); 
      transition: opacity 0.3s; 
      }
    
    .box-header .box-title {
      font-size: 1.5rem;       /* use rem for relative sizing */
      white-space: normal;     /* allow wrapping */
      }
    
    .param-container .box-header {
    padding-top: 7px !important;
    padding-bottom: 7px !important;
    }
    
    .param-container .control-label {
    font-size: 1.2rem !important; 
    margin-bottom: 2px;
    }
    
    .selectize-dropdown-content .option {
    min-height: 1.5em;
    }
    
    /* 1. Reduce gap with side bar */
    .plot-container .row {
    margin-left: -10px !important;
    margin-right: -10px !important;
    }
    
    
    /* 2. Reduce gap between side-by-side boxes */
    .plot-container .row > [class*='col-'] {
    padding-left: 0px !important;   /* Default is 15px */
    padding-right: 0px !important;
    }
    
    
    /* 3. Keep headers slim */
    .plot-container .box-header {
    padding-top: 5px !important;
    padding-bottom: 5px !important;
    }
    
    .plot-container .box-body {
    padding: 7px !important;
    }
    
    #reset_lo, #reset_sp, #reset_pc, #sample_params_lo, #sample_params_sp, #sample_params_pc {
     font-family: 'Arial', sans-serif;  
     font-size: 1.1rem;
     font-weight: bold
     }
     
    #click_lo , #click_sp , #click_pc {
     background-color: #DDE8D5;
     color: #408A0B; 
     font-family: 'Arial', sans-serif;  
     font-size: 1.5rem;
     font-weight: bold;
     border-color: #408A0B;
     }
    
    /* Force 'message' type notifications to be green */
      .shiny-notification-message {
        background-color: #28a745 !important;
        color: white !important;
        border: 1px solid #1e7e34 !important;
        opacity: 0.8 !important;
      }
      
    .wrapper {
      background-color: #f5f5f5 !important;
    }
    
    .content-wrapper, .right-side {
      background-color: #f5f5f5 !important;
      min-height: 100vh !important; /* force full height */
    }
    
    .box.box-solid.box-custom > .box-header {
      background-color: #800080 !important;
      color: white !important;
    }
    
    .box.box-solid.box-custom {
      border-color: #800080 !important;
    }
    
    #result_table table.dataTable {
    font-size: 1.2rem !important;
    }
    
    "))
  )
)