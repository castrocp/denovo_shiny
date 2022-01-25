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
  
  ### Function to generate the text to be shown in the selection menu. Percentile along with corresponding threshold value.
  paste_select_option <- function(percentile, value){
    output_option <- paste(percentile, paste("( \u2265 ", value, ")", sep = ""), sep = " ") #unicode, greater than or equal to    
    return(output_option)
  }
  
  ### Function to calculate quantile thresholds
  calc_quantile <- function(col, percentage){
    quant <- round(quantile(col, percentage, names = FALSE, na.rm = TRUE), digits = 3) #set names to false to return quantile value without percentage label
    #na.rm = TRUE is for the cases in which factors were converted to numeric type and NAs were introduced
    return(quant)
  }
  
  output$ui_select_child <- renderUI({
    selectInput(inputId = "select_child", label = "Child", choices = c("all", "proband","sibling"))
  })
  
  output$ui_select_regdb <- renderUI({
    selectInput(inputId = "select_regdb", label = "RegulomeDB", choices = c("all", "all 2s", "all 3s", "2a", "2b", "2c", "3a", "3b", "4", "5", "6", "7"))
  })
  
  output$ui_select_turf <- renderUI({
    selectInput(inputId = "select_turf", label = "TURF", choices = c("all", .6, .7, .8, .9))
  })
  
  ## Calculate brain-specific functional score quantiles
  brainsp_top_1_percent <- calc_quantile(DNM_df$brainSp_score, .99)
  brainsp_top_5_percent <- calc_quantile(DNM_df$brainSp_score, .95)   
  brainsp_top_10_percent <- calc_quantile(DNM_df$brainSp_score, .90)   
  brainsp_top_15_percent <- calc_quantile(DNM_df$brainSp_score, .85)   
  brainsp_top_20_percent <- calc_quantile(DNM_df$brainSp_score, .80)   
  brainsp_top_25_percent <- calc_quantile(DNM_df$brainSp_score, .75)   
  
  pasted_brainsp_top_1_percent <- paste_select_option("top 1%", brainsp_top_1_percent)
  pasted_brainsp_top_5_percent <- paste_select_option("top 5%", brainsp_top_5_percent)
  pasted_brainsp_top_10_percent <- paste_select_option("top 10%", brainsp_top_10_percent)
  pasted_brainsp_top_15_percent <- paste_select_option("top 15%", brainsp_top_15_percent)
  pasted_brainsp_top_20_percent <- paste_select_option("top 20%", brainsp_top_20_percent)
  pasted_brainsp_top_25_percent <- paste_select_option("top 25%", brainsp_top_25_percent)
  
  
  output$ui_select_brainscore <- renderUI({
    selectInput(inputId = "select_brainscore", label = "Brain-specific functional score", 
                choices = c("all",
                            c(paste_select_option("top 1%", brainsp_top_1_percent), paste_select_option("top 5%", brainsp_top_5_percent),
                              paste_select_option("top 10%", brainsp_top_10_percent), paste_select_option("top 15%", brainsp_top_15_percent),
                              paste_select_option("top 20%", brainsp_top_20_percent), paste_select_option("top 25%", brainsp_top_25_percent)
                            )
                          )
    )
  }) #"top 5% (.414)", "top 10% (.296)", "top 15% (.234)", "top 20% (.204)", "top 25% (.182)"))

  
  output$ui_select_CADD <- renderUI({
    selectInput(inputId = "select_cadd", label = "CADD score", 
                choices = c("all", "top 1% (21.1)", "top 5% (12.7)", "top 10% (9.4)", "top 15% (7.5)", "top 20% (6.2)", "top 25% (5.2)"))
  })
  
  output$ui_select_VEP <- renderUI({
    selectInput(inputId = "select_vep", label = "VEP", choices = c("all", as.character(unique(DNM_df$VEP))) )
  })
  
  ## Calculate phastCons score quantiles
  ## Convert from factor type to numeric
  DNM_df$phastCons <- as.numeric(as.character(DNM_df$phastCons))
  
  phastcons_top_1_percent <- calc_quantile(DNM_df$phastCons, .99)
  phastcons_top_5_percent <- calc_quantile(DNM_df$phastCons, .95)   
  phastcons_top_10_percent <- calc_quantile(DNM_df$phastCons, .90)   
  phastcons_top_15_percent <- calc_quantile(DNM_df$phastCons, .85)   
  phastcons_top_20_percent <- calc_quantile(DNM_df$phastCons, .80)   
  phastcons_top_25_percent <- calc_quantile(DNM_df$phastCons, .75)   
  
  pasted_phastcons_top_1_percent <- paste_select_option("top 1%", phastcons_top_1_percent)
  pasted_phastcons_top_5_percent <- paste_select_option("top 5%", phastcons_top_5_percent)
  pasted_phastcons_top_10_percent <- paste_select_option("top 10%", phastcons_top_10_percent)
  pasted_phastcons_top_15_percent <- paste_select_option("top 15%", phastcons_top_15_percent)
  pasted_phastcons_top_20_percent <- paste_select_option("top 20%", phastcons_top_20_percent)
  pasted_phastcons_top_25_percent <- paste_select_option("top 25%", phastcons_top_25_percent)
  
  output$ui_select_phastcons <- renderUI({
    selectInput(inputId = "select_phastcons", label = "phastCons score",
                choices = c("all",
                            c(pasted_phastcons_top_1_percent, pasted_phastcons_top_5_percent,
                              pasted_phastcons_top_10_percent, pasted_phastcons_top_15_percent,
                              pasted_phastcons_top_20_percent, pasted_phastcons_top_25_percent
                            )
                          )
    )
  })#     99%    95%    90%    85%    80%    75% 
    #    817.52 719.60 664.00 630.00 602.00 576.00 
  
  output$ui_select_enhancer <- renderUI({
    selectInput(inputId = "select_enhancer", label = "Fetal brain enhancer", choices = c("all", "yes", "no"))
  })
  
  output$ui_select_promoter <- renderUI({
    selectInput(inputId = "select_promoter", label = "Fetal brain promoter", choices = c("all", "yes", "no"))
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
    
    ### Brain-specific score selection
    ifelse(input$select_brainscore == "all", brainscore_selected <- DNM_df$brainSp_score,
          ifelse(input$select_brainscore == pasted_brainsp_top_1_percent, brainscore_selected <- brainsp_top_1_percent, 
                ifelse(input$select_brainscore == pasted_brainsp_top_5_percent, brainscore_selected <- brainsp_top_5_percent, 
                      ifelse(input$select_brainscore == pasted_brainsp_top_10_percent, brainscore_selected <- brainsp_top_10_percent,
                            ifelse(input$select_brainscore == pasted_brainsp_top_15_percent, brainscore_selected <- brainsp_top_15_percent, 
                                  ifelse(input$select_brainscore == pasted_brainsp_top_20_percent, brainscore_selected <- brainsp_top_20_percent, 
                                        brainscore_selected <- brainsp_top_25_percent)
                            )
                      )
                )
          )
    )
           
    
    ifelse(input$select_cadd == "all", cadd_selected <- DNM_df$CADD,
           ifelse(input$select_cadd == "top 1% (21.1)", cadd_selected <- 21.1, 
                  ifelse(input$select_cadd == "top 5% (12.7)", cadd_selected <- 12.7, 
                         ifelse(input$select_cadd == "top 10% (9.4)", cadd_selected <- 9.4,
                                ifelse(input$select_cadd == "top 15% (7.5)", cadd_selected <- 7.5, 
                                       ifelse(input$select_cadd == "top 20% (6.2)", cadd_selected <- 6.2, 
                                              cadd_selected <- 5.2) #top 25% cutoff (5.2)
                                )
                         )
                  )
           )
    )
    
    ### phastCons score selection
    ifelse(input$select_phastcons == "all", phastcons_selected <- "all",
           ifelse(input$select_phastcons == pasted_phastcons_top_1_percent, phastcons_selected <- phastcons_top_1_percent, 
                  ifelse(input$select_phastcons == pasted_phastcons_top_5_percent, phastcons_selected <- phastcons_top_5_percent, 
                         ifelse(input$select_phastcons == pasted_phastcons_top_10_percent, phastcons_selected <- phastcons_top_10_percent,
                                ifelse(input$select_phastcons == pasted_phastcons_top_15_percent, phastcons_selected <- phastcons_top_15_percent, 
                                       ifelse(input$select_phastcons == pasted_phastcons_top_20_percent, phastcons_selected <- phastcons_top_20_percent, 
                                              phastcons_selected <- phastcons_top_25_percent)
                                )
                         )
                  )
           )
    )
    
    ifelse(input$select_enhancer != "all",
           ifelse(input$select_enhancer == "yes", enhancer_selected <- "DHS_fetal_brain_enh", 
                  enhancer_selected <- "."),
           enhancer_selected <- DNM_df$fetal_brain_enh_dhs
    )
    
    ifelse(input$select_promoter != "all",
           ifelse(input$select_promoter == "yes", promoter_selected <- "DHS_fetal_brain_prom", 
                  promoter_selected <- "."),
           promoter_selected <- DNM_df$fetal_brain_prom_dhs
    )
    
    # 2. Filter data
    filt_DNM_df <- DNM_df %>% filter(child == child_selected &
                                       regDB2.0 == regdb_selected &
                                       VEP == vep_selected & 
                                       TURF >= turf_selected &
                                       brainSp_score >= brainscore_selected &
                                       #(if (phastcons_selected=="all") (phastCons==T | is.na(phastCons)) else phastCons >= phastcons_selected) &
                                       (case_when(phastcons_selected == "all" ~ phastCons %in% DNM_df$phastCons,
                                                  TRUE ~ phastCons >= phastcons_selected)) &
                                       CADD >= cadd_selected &
                                       fetal_brain_enh_dhs == enhancer_selected &
                                       fetal_brain_prom_dhs == promoter_selected)
    
    # 3. Return filtered table
    filt_DNM_df   
  })
  
  # 3. Output filtered data table (uses DT package)
  output$dt_table <- renderDataTable({
    datatable(df_data(), options = list(pageLength = 15)) #options = list(lengthMenu = c(5, 30, 50), pageLength = 5)
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
  rnames <- c("probands", "siblings")
  cnames <- c("meeting criteria","not meeting criteria")
  twobytwo <- reactive({
    matrix(c(filtered_proband_count(), filtered_sibling_count(), filtered_no_proband_count(), filtered_no_sibling_count()), nrow = 2,
           dimnames= list(rnames, cnames)
    )
  }) 
  
  # Output counts matrix object
  output$counts <-renderTable({
    twobytwo()
  }, rownames = TRUE)
  
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
      width = 3,
      h4("Select desired filters for mutations to be shown in the table"),
      br(),
      title = "Inputs",
      
      uiOutput("ui_select_child"),
      uiOutput("ui_select_regdb"),
      uiOutput("ui_select_turf"),
      uiOutput("ui_select_brainscore"),
      uiOutput("ui_select_CADD"),
      uiOutput("ui_select_VEP"),
      uiOutput("ui_select_phastcons"),
      uiOutput("ui_select_enhancer"),
      uiOutput("ui_select_promoter")
      
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
  theme = shinytheme('cerulean'),
  main_page,
  about_page
)

### 4. Run app
shinyApp(ui = ui, server = server)