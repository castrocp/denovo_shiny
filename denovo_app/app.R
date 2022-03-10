library(shiny)
library(here)
library(tidyverse)
library(DT)
library(shinythemes)

### 1. Dataset

DNM_df <- read.delim(here("denovo_app/data", "DNMs.hg19.indelFilt.rptFilt.MAF001.singletons.RegDBv2.TURF.BrainSpecific.CADD.VEP.phastCons.SIFT.PolyPhen.DHS_fetal_brain_enh.DHS_fetal_brain_prom.1500bp_prom.autosomes.bed"))

### 2. Server
server <- function(input, output, session) {
  
  ########## 1. UI dropdown selection ##########
  
  ### Function to generate the text to be shown in the selection menu. Percentile along with corresponding threshold value.
  paste_select_option <- function(percentile, value){
    output_option <- paste(percentile, paste("( \u2265 ", value, ")", sep = ""), sep = " ") #unicode, greater than or equal to    
    return(output_option)
  }
  
  ### Function to calculate quantile thresholds
  calc_quantile <- function(col, percentile){
    quant <- round(quantile(col, percentile, names = FALSE, na.rm = TRUE), digits = 3) #set names to false to return quantile value without percentile label
    #na.rm = TRUE is for the cases in which factors were converted to numeric type and NAs were introduced
    return(quant)
  }
  
  ### Child selection
  output$ui_select_child <- renderUI({
    selectInput(inputId = "select_child", label = "Child", choices = c("all", "proband","sibling"))
  })
  
  ### RegulomeDB score selection
  output$ui_select_regdb <- renderUI({
    selectInput(inputId = "select_regdb", label = "RegulomeDB", choices = c("all", "all 2s", "all 3s", "2a", "2b", "2c", "3a", "3b", "4", "5", "6", "7"))
  })
  
  ### TURF score selection and quantile calculation
  turf_top_1_percent <- calc_quantile(DNM_df$TURF, .99)
  turf_top_5_percent <- calc_quantile(DNM_df$TURF, .95)
  turf_top_10_percent <- calc_quantile(DNM_df$TURF, .90)
  turf_top_15_percent <- calc_quantile(DNM_df$TURF, .85)
  turf_top_20_percent <- calc_quantile(DNM_df$TURF, .80)
  turf_top_25_percent <- calc_quantile(DNM_df$TURF, .75)
  
  pasted_turf_top_1_percent <- paste_select_option("top 1%", turf_top_1_percent)
  pasted_turf_top_5_percent <- paste_select_option("top 5%", turf_top_5_percent)
  pasted_turf_top_10_percent <- paste_select_option("top 10%", turf_top_10_percent)
  pasted_turf_top_15_percent <- paste_select_option("top 15%", turf_top_15_percent)
  pasted_turf_top_20_percent <- paste_select_option("top 20%", turf_top_20_percent)
  pasted_turf_top_25_percent <- paste_select_option("top 25%", turf_top_25_percent)
  
  output$ui_select_turf <- renderUI({
    selectInput(inputId = "select_turf", label = "TURF",
                choices = c("all",
                            c(pasted_turf_top_1_percent, pasted_turf_top_5_percent,
                              pasted_turf_top_10_percent, pasted_turf_top_15_percent,
                              pasted_turf_top_20_percent, pasted_turf_top_25_percent)
                           )
               )
                           
  })
  
  ### Brain-specific functional score selection and quantile calculation
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
                            c(pasted_brainsp_top_1_percent, pasted_brainsp_top_5_percent,
                              pasted_brainsp_top_10_percent, pasted_brainsp_top_15_percent,
                              pasted_brainsp_top_20_percent, pasted_brainsp_top_25_percent) #"top 5% (.414)", "top 10% (.296)", "top 15% (.234)", "top 20% (.204)", "top 25% (.182)"))
                           )
                )
  }) 

  ### CADD score selection and quantile calculation
  cadd_top_1_percent <- calc_quantile(DNM_df$CADD, .99)
  cadd_top_5_percent <- calc_quantile(DNM_df$CADD, .95)
  cadd_top_10_percent <- calc_quantile(DNM_df$CADD, .90)
  cadd_top_15_percent <- calc_quantile(DNM_df$CADD, .85)
  cadd_top_20_percent <- calc_quantile(DNM_df$CADD, .80)
  cadd_top_25_percent <- calc_quantile(DNM_df$CADD, .75)
  
  pasted_cadd_top_1_percent <- paste_select_option("top 1%", cadd_top_1_percent)
  pasted_cadd_top_5_percent <- paste_select_option("top 5%", cadd_top_5_percent)
  pasted_cadd_top_10_percent <- paste_select_option("top 10%", cadd_top_10_percent)
  pasted_cadd_top_15_percent <- paste_select_option("top 15%", cadd_top_15_percent)
  pasted_cadd_top_20_percent <- paste_select_option("top 20%", cadd_top_20_percent)
  pasted_cadd_top_25_percent <- paste_select_option("top 25%", cadd_top_25_percent)
  
  output$ui_select_CADD <- renderUI({
    selectInput(inputId = "select_cadd", label = "CADD score", 
                choices = c("all",
                            c(pasted_cadd_top_1_percent, pasted_cadd_top_5_percent,
                              pasted_cadd_top_10_percent, pasted_cadd_top_15_percent,
                              pasted_cadd_top_20_percent, pasted_cadd_top_25_percent) #"top 1% (21.1)", "top 5% (12.7)", "top 10% (9.4)", "top 15% (7.5)", "top 20% (6.2)", "top 25% (5.2)")
                )
    )
  })
  
  ### VEP category selection 
  output$ui_select_VEP <- renderUI({
    selectInput(inputId = "select_vep", label = "VEP", choices = c("all", as.character(unique(DNM_df$VEP))) )
  })
  
  ### phastCons score selection and quantile calculation
  DNM_df$phastCons <- as.numeric(as.character(DNM_df$phastCons))   #Convert from factor type to numeric because of NAs

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
                              pasted_phastcons_top_10_percent, pasted_phastcons_top_15_percent,  #     99%    95%    90%    85%    80%    75% 
                              pasted_phastcons_top_20_percent, pasted_phastcons_top_25_percent)  #    817.52 719.60 664.00 630.00 602.00 576.00 
                )
    )
  })
  
  ### SIFT selection UI
  
  
  
  ### Enhancer selection
  output$ui_select_enhancer <- renderUI({
    selectInput(inputId = "select_enhancer", label = "Fetal brain enhancer", choices = c("all", "yes", "no"))
  })
  
  ### Promoter selection
  output$ui_select_promoter <- renderUI({
    selectInput(inputId = "select_promoter", label = "Fetal brain promoter", choices = c("all", "yes", "no"))
  })
  
  ###  promoter selection
  output$ui_select_promoter_by_proximity <- renderUI({
    selectInput(inputId = "select_promoter_by_proximity", label = "Promoter by TSS proximity", choices = c("all", "yes", "no"))
  })
  
  
  ####### 2. Create reactive data set based on user selections
  df_data <- reactive({
    
    # 1. Read UI selections
    
    ### Select child
    ifelse(input$select_child != "all", child_selected <- input$select_child, child_selected <- DNM_df$child)
    
    ### Select RegulomeDB score
    ifelse(input$select_regdb == "all", regdb_selected <- DNM_df$regDB2.0,
           ifelse(input$select_regdb == "all 2s", regdb_selected <- c("2a","2b","2c"),
                  ifelse(input$select_regdb == "all 3s", regdb_selected <- c("3a","3b"),
                         regdb_selected <- input$select_regdb
                  )
           )
    )    
    
    ### TURF selection
    ifelse(input$select_turf == "all", turf_selected <- DNM_df$TURF,
           ifelse(input$select_turf == pasted_turf_top_1_percent, turf_selected <- turf_top_1_percent, 
                  ifelse(input$select_turf == pasted_turf_top_5_percent, turf_selected <- turf_top_5_percent, 
                         ifelse(input$select_turf == pasted_turf_top_10_percent, turf_selected <- turf_top_10_percent,
                                ifelse(input$select_turf == pasted_turf_top_15_percent, turf_selected <- turf_top_15_percent, 
                                       ifelse(input$select_turf == pasted_turf_top_20_percent, turf_selected <- turf_top_20_percent, 
                                              turf_selected <- turf_top_25_percent)
                                )
                         )
                  )
           )
    )
    
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
    
    ### CADD score selection       
    ifelse(input$select_cadd == "all", cadd_selected <- DNM_df$CADD,
           ifelse(input$select_cadd == pasted_cadd_top_1_percent, cadd_selected <- cadd_top_1_percent,
                   ifelse(input$select_cadd == pasted_cadd_top_5_percent, cadd_selected <- cadd_top_5_percent,
                           ifelse(input$select_cadd == pasted_cadd_top_10_percent, cadd_selected <- cadd_top_10_percent,
                                   ifelse(input$select_cadd == pasted_cadd_top_15_percent, cadd_selected <- cadd_top_15_percent,
                                           ifelse(input$select_cadd == pasted_cadd_top_20_percent, cadd_selected <- cadd_top_20_percent,
                                                   cadd_selected <- cadd_top_25_percent)
                                   )
                           )
                   )
           )
    )
    
    ### VEP category selection
    ifelse(input$select_vep != "all",
           vep_selected <- input$select_vep,
           vep_selected <- DNM_df$VEP) 
    
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
    
    ### Enhancer selection
    ifelse(input$select_enhancer != "all",
           ifelse(input$select_enhancer == "yes", enhancer_selected <- "DHS_fetal_brain_enh", 
                  enhancer_selected <- "."),
           enhancer_selected <- DNM_df$fetal_brain_enh_dhs
    )
    
    ### Promoter selection
    ifelse(input$select_promoter != "all",
           ifelse(input$select_promoter == "yes", promoter_selected <- "DHS_fetal_brain_prom", 
                  promoter_selected <- "."),
           promoter_selected <- DNM_df$fetal_brain_prom_dhs
    )
    
    ### Promoter-by-proximity selection
    ifelse(input$select_promoter_by_proximity != "all",
           ifelse(input$select_promoter_by_proximity == "no", promoter_by_proximity_selected <- ".", 
                  promoter_by_proximity_selected <- "yes"),
           promoter_by_proximity_selected <- DNM_df$promoter_1500bp
    )
    
    
    # 2. Filter data
    filt_DNM_df <- DNM_df %>% filter(child %in% child_selected &
                                       regDB2.0 %in% regdb_selected &
                                       VEP == vep_selected & 
                                       TURF >= turf_selected &
                                       brainSp_score >= brainscore_selected &
                                       (case_when(phastcons_selected == "all" ~ phastCons %in% DNM_df$phastCons, #return all NA rows also
                                                  TRUE ~ phastCons >= phastcons_selected)) &
                                       CADD >= cadd_selected &
                                       fetal_brain_enh_dhs == enhancer_selected &
                                       fetal_brain_prom_dhs == promoter_selected &
                                       (case_when(promoter_by_proximity_selected == "yes" ~ promoter_1500bp != ".",
                                                  TRUE ~ promoter_1500bp == promoter_by_proximity_selected))
                                     )
                                       #promoter_1500bp == promoter_by_proximity_selected)
    
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
  
  # Calculate proportions
  pro_prop <- reactive({
    filtered_proband_count() / total_proband
  })
  
  sib_prop <- reactive({
    filtered_sibling_count() / total_sibling
  })
  
  # Output proportions
  output$proband_prop <- renderText({
    pro_prop()
  })
  output$sibling_prop <- renderText({
    sib_prop()
  })
  
  # Output p-value from the fisher's exact test
  output$pval <- renderText({
    result <- fisher.test(twobytwo(), alternative = "greater")$p.value
    paste('The unadjusted p-value is:', result)
  })
  
}


### 3. UI
main_page <- tabPanel(
  title = "De novo table",
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Filter mutations by feature"),
      title = "Inputs",
      
      uiOutput("ui_select_child"),
      uiOutput("ui_select_regdb"),
      uiOutput("ui_select_turf"),
      uiOutput("ui_select_brainscore"),
      uiOutput("ui_select_CADD"),
      uiOutput("ui_select_VEP"),
      uiOutput("ui_select_phastcons"),
      uiOutput("ui_select_enhancer"),
      uiOutput("ui_select_promoter"),
      uiOutput("ui_select_promoter_by_proximity")
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "De Novo Mutation table",
          DTOutput('dt_table') #use DT's DTOutput() to avoid conflict with shiny::dataTableOutput()
        ),
        tabPanel(
          title = "Comparison of proband vs. sibling",
          h3("The table below displays mutation counts across probands and siblings,
             based on mutations meeting the selected-filter criteria"),
          br(),
          tableOutput("counts"),
          h3("Proportions for mutations of interest"),
          strong("Proportion of selected mutations from total in probands: "), textOutput("proband_prop"),
          br(),
          strong("Proportion of selected mutations from total in siblings: "), textOutput("sibling_prop"),
          h3("Test for enrichment of mutations represented in the above table"),
          "Using Fisher's exact test and the counts from the table above.",
          br(),
          "Test for enrichment of mutations in probands:",
          br(),
          br(),
          textOutput("pval")
        )
      )
    )
  )
)   
      

about_page <- tabPanel(
  title = "About",
  titlePanel("About"),
  "Created with R Shiny",
  br(),
  "2021"
)
  
ui <- navbarPage( 
  title = 'De Novo Browser',
  theme = shinytheme('cerulean'),
  main_page,
  about_page
)

### 4. Run app
shinyApp(ui = ui, server = server)