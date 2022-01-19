library(shiny)
library(here)
library(tidyverse)
library(DT)
library(shinythemes)

### 1. Dataset

DNM_df <- read.delim(here("denovo_app/data", "DNMs.hg19.indelFilt.rptFilt.MAF001.singletons.RegDBv2.TURF.BrainSpecific.CADD.VEP.phastCons.SIFT.PolyPhen.DHS_fetal_brain_enh.DHS_fetal_brain_prom.1500bp_prom.autosomes.bed"))

### 2. Server
server <- function(input, output, session) {
  
  # 1. UI dropdown selection
  output$ui_select_child <- renderUI({
    selectInput(inputId = "select_child", label = "Child", choices = c("all", "proband","sibling"))
  })
  
  output$ui_select_regdb <- renderUI({
    selectInput(inputId = "select_regdb", label = "RegulomeDB", choices = c("all", "all 2s", "all 3s", "2a", "2b", "2c", "3a", "3b", "4", "5", "6", "7"))
  })
  
  output$ui_select_VEP <- renderUI({
    selectInput(inputId = "select_vep", label = "VEP", choices = c("all", as.character(unique(DNM_df$VEP))) )
  })
  
  output$ui_select_turf <- renderUI({
    selectInput(inputId = "select_turf", label = "TURF", choices = c("all", .6, .7, .8, .9))
  })
  
  output$ui_select_enhancer <- renderUI({
    selectInput(inputId = "select_enhancer", label = "Fetal brain enhancer", choices = c("all", "yes", "no"))
  })
  
  # 2. Reactive data set
  df_data <- reactive({
    
    # 1. Read UI selections
    ifelse(input$select_child != "all", child_selected <- input$select_child, child_selected <- DNM_df$child)
    
    # THIS DOESNT SEEM TO BE ACCURATELY PULLING ALL THE 2S AND 3S. THE COUNTS TO DONT ADD UP
    ifelse(input$select_regdb == "all", regdb_selected <- DNM_df$regDB2.0,
           ifelse(input$select_regdb == "all 2s", regdb_selected <- c("2a","2b","2c"),
                  ifelse(input$select_regdb == "all 3s", regdb_selected <- c("3a","3b"),
                         regdb_selected <- input$select_regdb
                  )
           )
    )    
    
    ifelse(input$select_vep != "all",
           vep_selected <- input$select_vep,
           vep_selected <- DNM_df$VEP)  
    
    ifelse(input$select_turf != "all",
           turf_selected <- input$select_turf,
           turf_selected <- DNM_df$TURF)
    
    ifelse(input$select_enhancer != "all",
           ifelse(input$select_enhancer == "yes", enhancer_selected <- "DHS_fetal_brain_enh", 
                  enhancer_selected <- "."),
           enhancer_selected <- DNM_df$fetal_brain_enh_dhs
    )
    
    # 2. Filter data
    filt_DNM_df <- DNM_df %>% filter(child == child_selected &
                                       regDB2.0 == regdb_selected &
                                       VEP == vep_selected & 
                                       TURF >= turf_selected &
                                       fetal_brain_enh_dhs == enhancer_selected)
    
    # 3. Return filtered table
    filt_DNM_df   
  })
  
  # 3. Output filtered data table (uses DT package)
  output$dt_table <- renderDataTable({
    datatable(df_data())
  })    
  
  # 4. Create mutation count table
  
  # Get total number of mutations in probands and siblings
  total_proband <- nrow(filter(DNM_df, child == "proband")) 
  total_sibling <- nrow(filter(DNM_df, child == "sibling"))
  
  # Reactive count of probands and siblings with and without mutations in the filtered table
  filtered_proband_count <- reactive({
    nrow(filter(df_data(), child == "proband"))
  })
  
  filtered_no_proband_count <- reactive({
    total_proband - filtered_proband_count()
  })
  
  filtered_sibling_count <- reactive({
    nrow(filter(df_data(), child == "sibling"))
  })
  
  filtered_no_sibling_count <- reactive({
    total_sibling - filtered_sibling_count()
  })
  
  # Create reactive matrix object with reactive counts
  twobytwo <- reactive({
    matrix(c(filtered_proband_count(), filtered_sibling_count(), filtered_no_proband_count(), filtered_no_sibling_count()), nrow = 2,
           dimnames= list(
             child=c('proband','sibling'),
             feature=c("meeting criteria","not meeting criteria")))
  }) 
  
  # Output counts matrix object
  output$counts <-renderTable({
    twobytwo()
  })
  
  # Output p-value from the fisher's exact test
  output$pval <- renderText({
    result <- fisher.test(twobytwo(), alternative = "greater")$p.value
    paste('Test for enrichment in probands of mutations meeting criteria. The p-value is:', result)
  })
  
}


### 3. UI
main_page <- tabPanel(
  title = "De novo table",
  
  sidebarLayout(
    sidebarPanel(
      h4("Select desired filters for mutations to be shown in the table"),
      br(),
      title = "Inputs",
      
      uiOutput("ui_select_child"),
      uiOutput("ui_select_regdb"),
      uiOutput("ui_select_VEP"),
      uiOutput("ui_select_turf"),
      uiOutput("ui_select_enhancer")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "De Novo Mutation table",
          DTOutput('dt_table') #use DT's DTOutput() to avoid conflict with shiny::dataTableOutput()
        ),
        tabPanel(
          title = "Table of counts",
          #plotOutput("plot_1")
          tableOutput("counts"),
          textOutput("pval")
        )
      )
    )
  )
)   
      
#   fluidPage(
#     fluidRow(
#       # Column widths should add up to 12 within a fluidRow() container
#       column(4,
#              uiOutput("ui_select_child")
#       ),
#       column(4,
#              uiOutput("ui_select_regdb")
#       ),
#       column(4,
#              uiOutput("ui_select_VEP")
#       )
#     ),  
#     fluidRow(
#       column(4,
#              uiOutput("ui_select_turf")
#       ),
#       column(4,
#              uiOutput("ui_select_enhancer")
#       )
#     ),
# 
#     fluidRow(dataTableOutput("dt_table"))
#   )
# )
# 
# counts_page <- tabPanel(
#   title = "Counts and Enrichment Testing",
#   tableOutput("counts"),
#   
#   textOutput("pval")
# )

about_page <- tabPanel(
  title = "About",
  titlePanel("About"),
  "Created with R Shiny",
  br(),
  "2021"
)
  
ui <- navbarPage( 
  title = 'De Novo Analyser',
  theme = shinytheme('united'),
  main_page,
  about_page
)

### 4. Run app
shinyApp(ui = ui, server = server)