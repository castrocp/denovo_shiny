library(shiny)
library(shinythemes)
library(DT)
library(here)
library(tidyverse)

DNM_df <- read.delim(here("denovo_app/data", "DNMs.hg19.indelFilt.rptFilt.MAF001.singletons.RegDBv2.TURF.BrainSpecific.CADD.VEP.phastCons.SIFT.PolyPhen.DHS_fetal_brain_enh.DHS_fetal_brain_prom.1500bp_prom.autosomes.bed"))
DNM_df$fetal_brain_prom_dhs <- unlist(lapply(DNM_df$fetal_brain_prom_dhs, function(x) {return(x=="DHS_fetal_brain_prom")}))
DNM_df$fetal_brain_enh_dhs <- unlist(lapply(DNM_df$fetal_brain_enh_dhs, function(x) {return(x=="DHS_fetal_brain_enh")}))

# Define UI -----------
# Shiny uses a family of functions that turn R objects into output for the user interface
# Each function creates a specific type of output: plotOutput, tableOutput, textOuput, etc.
# Each `--Output` function takes as its first argument the name of the element, which is created on the server side.

not_sel <- "Not Selected" # to show when a selection has not yet been made in the drop-down menu

main_page <- tabPanel(
  title = "Analysis",
  titlePanel("Testing Fisher's Exact Test"),
  
  sidebarLayout(
    sidebarPanel(
      title = "Inputs",
    
      # The `input` argument for the server section function is a list-like object that stores the values of widgets
      # The value is saved under the name used for `inputID` and later accessed by the server as `input$inputID`
      # Widget inputs are also referred to as reactive values
      
      #selectInput(inputId = "sel_feature", label = "feature", choices = names(DNM_df)),
      #selectInput(inputId = "sel_feature", label = "Feature", choices = "child", selected = "child"),
      #selectInput("num_var_2", "Additional Genomic Feature", choices = c(not_sel, "TURF > .8", "RegulomeDB 2s")),
      
      uiOutput("dropdownlist"),
      
      br(),
      actionButton("add_feature", "Add Feature", icon = icon("plus")),
      br(),
      br(),
      actionButton("run_button", "Run Test", icon = icon("play")),
      uiOutput('addfeature')
      
      
  ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "De Novo Mutation table",
          DTOutput('df') #use DT's DTOutput() to avoid conflict with shiny::dataTableOutput()
        ),
        tabPanel(
          title = "Table of counts",
          #plotOutput("plot_1")
          tableOutput("counts"),
          textOutput("pval")
        ),
        tabPanel(
          title = "test linking of features to dataframe",
          
          textOutput("proband_count"),
          textOutput("sibling_count"),
          textOutput("proband_count2"),
          textOutput("test")
          #textOutput("additional_selected_feature")
          
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
  title = 'De Novo Analyser',
  theme = shinytheme('united'),
  main_page,
  about_page
)

  

# Define server logic -----------
# You provide instructions on where to display object in the UI section.
# You build the objects under the server function.
# Each object is placed in a list-like object named `output`
# Each object should contain the output of one of Shiny's `render--` functions. ex. renderPlot, renderTable, etc.
# The element name should match the name that's used in the UI section.

server <- function(input, output) {
  
  # Fisher's exact test
  fet <- function(x){
    #total_proband <- 30
    #total_sibling <- 30
    #r1c1 <- 15 #proband with
    #r1c2 <- total_proband - r1c1 #proband without
    #r2c1 <- 13 #sibling with
    #r2c2 <- total_sibling - r2c1 #sibling without
    
    #twobytwo <- matrix(c(r1c1,r1c2,r2c1,r2c2), nrow = 2)
    
    fet_output <- fisher.test(x)
    return(fet_output)
  }
  
  
  # Create a reactive dataframe based on selected input
  #data_input <- reactive({ # brackets are there to allow multiple rows of code. Otherwise can omit.
  #  req(input$sel_feature) # req makes sure we have the necessary input
  #  filter(DNM_df, input$sel_feature=="proband")
  #})

  
  # Get total number of mutations in probands and siblings
  total_proband <- nrow(filter(DNM_df, child == "proband"))
  total_sibling <- nrow(filter(DNM_df, child == "sibling"))
  
  output$proband_count2 <- renderText({
    total_proband
    #nrow(filter(DNM_df, get(input$sel_feature) == "13123"))
  })
  
  # Can get the number of entries in dataframe dependent on column values
  output$proband_count <- renderText({
    #total_proband
    nrow(filter(DNM_df, get(input$sel_feature) == "DHS_fetal_brain_prom"))
  })
  
  output$sibling_count <- renderText({
    #total_sibling
    nrow(filter(DNM_df, get(input$sel_feature) != "DHS_fetal_brain_prom"))
  })
  
  # Another possible method for accessing values using selected column name  
    #switch(input$sel_feature,
    #       chrom = {nrow(filter(DNM_df, famID == "13850"))}, 
    #       start = {"b"},
    #       {"default"})


  

  # Populate 2x2 matrix with counts of mutations corresponding to selected feature
  get.table <- function(input){
    
    ## Might make sense to use a dictionary to populate the conditional statement
    ## dependent on which feature the user selects
    
    ## ex. if the feature is fetal_brain_prom_dhs, then you're looking for
    ## the "DHS_fetal_brain_prom" annotation in the data table to count.
    ## So the pair would be fetal_brain_prom_dhs and "DHS_fetal_brain_prom"
    
    ##if USER-SELECTED FEATURE IS IN DICTIONARY:
    ##  RETURN THE VALUE CORRESPONDING TO THAT SELECTION
    
    #if user selects regDB2.0, run:  get(input$sel_feature) == "2a"
    #if user selects TURF, run: get(input$sel_feature) >= .89 (top 1% score)
    
    
    ## USE BOOLEAN AND "WHICH" FOR COUNTING FEATURES
    
    # mutations in proband for selected feature
    #pro_with <-  nrow(filter(DNM_df,
    #                         child == "proband"
    #                         & get(input$sel_feature) == "DHS_fetal_brain_prom"))
    pro_with <- length(which(DNM_df[,input$sel_feature] == TRUE & DNM_df$child == "proband"))
    
    # mutations in proband without selected feature present
    pro_wo <- total_proband - pro_with
    
    # mutations in sibling for selected feature
    #sib_with <- nrow(filter(DNM_df,
    #                        child == "sibling"
    #                        & get(input$sel_feature) == "DHS_fetal_brain_prom")) # mutations in sibling for selected feature
    sib_with <- length(which(DNM_df[,input$sel_feature] == TRUE & DNM_df$child == "sibling"))
    
    # mutations in siblings without selected feature present
    sib_wo <- total_sibling - sib_with
    
    twobytwo <- matrix(
      c(pro_with, sib_with, pro_wo, sib_wo), nrow = 2,
      dimnames= list(
        child=c('proband','sibling'),
        feature=c("with feature","without feature")))
  }


    
    
  
    

  


  
  
  
  
  #output$plot <- renderPlot({
    #barplot(data_input())
    #})
  
  #observeEvent(data_input(),{ # function to execute code according to a changing input field
    #choices <- c(not_sel,names(data_input())) # choices is a list of the columns in the data.table
    #updateSelectInput(inputId = "num_var_1", choices = choices)
    #updateSelectInput(inputId = "num_var_2", choices = choices)
    #updateSelectInput(inputId = "fact_var", choices = choices)
  #})
  
  output$df = renderDT(DNM_df, rownames=FALSE)
  
  # Show the 2x2 table of counts
  output$counts <- renderTable({
    get.table(input)
    #x <- get.table(input)
    #y <- as.data.frame(addmargins(x))
    #y[,c(1,2,3)] <- apply(y[,c(1,2,3)],1, as.integer)
    #y
  })
  
  # Show the p-value from the fisher's exact test
  output$pval <- renderText({
    x <- get.table(input)
    result <- fet(x)
    paste('p-value:', result) #"result" should reference the pvalue which is returned by the fet function
  })
  
  

# Use "()" after the reactive element variable name to get the value of the element
# Otherwise, it tries to interact with the reactive element itself and throws an error when trying to print the value  
  observeEvent(input$run_button,{
    feature1 <- reactive(input$prom_count)
    output$selected_feature <- renderText({
      paste('The feature that was selected is:', feature1())
    })
  })
  
  observeEvent(input$run_button,{
    feature2 <- reactive(input$num_var_1)
    output$additional_selected_feature <- renderText({
      paste('The second feature that was selected is:', feature2())
    })
  })
  
  
  
  dropdownslist <- reactiveValues()
  dropdownslist[[as.character(1)]] <- selectInput(inputId = "sel_feature_1", label = "feature 1", choices = names(DNM_df))
  
  dropdown_count <- reactiveVal(1)
    
  observeEvent(input$add_feature,{
    dropdown_count(dropdown_count() + 1)
    dropdownslist[[as.character(dropdown_count())]] <- selectInput(inputId = paste("sel_feature_", dropdown_count(), sep = ""), label = paste("feature", dropdown_count(), sep = " "), choices = names(DNM_df))
    
output$dropdownlist <- renderUI({dropdownslist})
    
    #print(dropdownslist)
    #print(dropdown_count())
  })
               
  output$dropdownlist <- renderUI({dropdownslist})
  
}
  

# Run the app -----------
shinyApp(ui = ui, server = server)

# reference tutorial from:
# https://towardsdatascience.com/how-to-build-a-data-analysis-app-in-r-shiny-143bee9338f7