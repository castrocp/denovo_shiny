library(shiny)
library(shinythemes)
library(DT)
library(here)

# using tutorial from:
# https://towardsdatascience.com/how-to-build-a-data-analysis-app-in-r-shiny-143bee9338f7

DNM_df <- read.delim(here("denovo_app/data", "DNMs.hg19.indelFilt.rptFilt.MAF001.singletons.RegDBv2.TURF.BrainSpecific.CADD.VEP.phastCons.SIFT.PolyPhen.DHS_fetal_brain_enh.DHS_fetal_brain_prom.1500bp_prom.autosomes.bed"))

not_sel <- "Not Selected" # to show when a selection has not yet been made in the drop-down menu
###
main_page <- tabPanel(
  title = "Analysis",
  titlePanel("Testing Fisher's Exact Test"),
  
  sidebarLayout(
    sidebarPanel(
      title = "Inputs",
      #fileInput("csv_input", "Select CSV File to Import", accept = ".csv"),
      selectInput(inputId = "prom_count", label = "feature", choices = names(DNM_df)),
      selectInput("num_var_1", "Genomic Feature", choices = c(not_sel, "fetal brain promoter", "promoter 1500bp upstream", "fetal brain enhancer")),
      selectInput("num_var_2", "Additional Genomic Feature", choices = c(not_sel, "TURF > .8", "RegulomeDB 2s")),

      br(),
      actionButton("add_feature", "Add Feature", icon = icon("plus")),
      br(),
      actionButton("run_button", "Run Test", icon = icon("play"))
      
      
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
          
          textOutput("selected_feature"),
          textOutput("additional_selected_feature")
          
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

# Define UI -----------
ui <- navbarPage( 
    title = 'De Novo Analyser',
    theme = shinytheme('united'),
    main_page,
    about_page
)

    

  

# Define server logic -----------

server <- function(input, output) {
  

  
  #Create the values to fill the 2x2 matrix in with
  total_proband <- 30
  total_sibling <- 30
  
  # I'll want to fill these numbers in
  get.table <- function(input){
    r1c1 <- 15 #proband with
    r1c2 <- total_proband - r1c1 #proband without
    r2c1 <- 13 #sibling with
    r2c2 <- total_sibling - r2c1 #sibling without
    
    # Create 2x2 matrix
    twobytwo <- matrix(
      c(r1c1,r1c2,r2c1,r2c2), nrow = 2,
      dimnames = list(
        child=c('proband','sibling'),
        feature=c('present','not present') #whatever feature we're looking for DNMs in
      )
    )
  }
  
  # Fisher's exact test.   DO I HAVE TO BUILD THE MATRIX WITHIN HERE OR CAN I JUST USE THE ONE BUILT ABOVE?
  fet <- function(x){
    total_proband <- 30
    total_sibling <- 30
    r1c1 <- 15 #proband with
    r1c2 <- total_proband - r1c1 #proband without
    r2c1 <- 13 #sibling with
    r2c2 <- total_sibling - r2c1 #sibling without
    
    twobytwo <- matrix(c(r1c1,r1c2,r2c1,r2c2), nrow = 2)
    
    fet_output <- fisher.test(twobytwo)
    return(fet_output)
  }
  
  
  
  data_input <- reactive({ # brackets are there to allow multiple rows of code. Otherwise can omit.
    req(input$csv_input) # req makes sure we have the necessary input
    fread(input$csv_input$datapath)
  })
  
  observeEvent(data_input(),{ # function to execute code according to a changing input field
    choices <- c(not_sel,names(data_input())) # choices is a list of the columns in the data.table
    updateSelectInput(inputId = "num_var_1", choices = choices)
    updateSelectInput(inputId = "num_var_2", choices = choices)
    #updateSelectInput(inputId = "fact_var", choices = choices)
  })
  
  output$df = renderDT(DNM_df, rownames=FALSE)
  
  # Show the 2x2 table of counts
  output$counts <- renderTable({
    x <- get.table(input)
    y <- as.data.frame(addmargins(x))
    y[,c(1,2,3)] <- apply(y[,c(1,2,3)],1, as.integer)
    y
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
  
  
  
  
}
  
#  output$selected_var <- renderText({
#    paste("You have selected", input$var)
#  })
  
#  output$min_max <- renderText({
#    paste("You chose a range that goes from", input$range[1], "to", input$range[2])
#  })
  
#  output$dnm_table = DT::renderDataTable({
#    read.delim("/data/DNMs.hg19.indelFilt.rptFilt.MAF001.RegDBv2.singletons.CADD.VEP.phastCons.SIFT.PolyPhen.fetal_brain_enhancer.DHS_fetal_brain_enh.DHS_fetal_brain_prom.turf_score.organScore.brainSp_score")
#  })



# Run the app -----------
shinyApp(ui = ui, server = server)