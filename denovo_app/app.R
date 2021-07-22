library(shiny)
library(DT)

#dnm_df <- read.delim("/Users/christopher/Documents/Autism Project/Denovo_Shiny/denovo_shiny/denovo_app/data/DNMs.hg19.indelFilt.rptFilt.MAF001.RegDBv2.singletons.CADD.VEP.phastCons.SIFT.PolyPhen.fetal_brain_enhancer.DHS_fetal_brain_enh.DHS_fetal_brain_prom.turf_score.organScore.brainSp_score.bed")



# Define UI -----------
ui <- fluidPage( 
  titlePanel("A Shiny app for viewing de novo mutations"),
  
  navbarPage(
    title = 'DataTable Options',
    tabPanel('This is a link',     DT::dataTableOutput('ex1')),
    tabPanel('Hi',        DT::dataTableOutput('ex2')),
    tabPanel('Something',      DT::dataTableOutput('ex3')),
    tabPanel('Some Stuff',       DT::dataTableOutput('ex4')),
    tabPanel('Something else',  DT::dataTableOutput('ex5'))
  ),
  
  DT::dataTableOutput("dnm_table"),
  
  
  sidebarLayout(
    sidebarPanel(
      h2("This is the sidebar panel"),
      p("I can put something here that allows the user to perform statistical tests"),
    
    helpText("Create demographic maps with 
               information from the 2010 US Census."),
    
    selectInput("var", 
                label = "Choose a variable to display",
                choices = list("Percent White", 
                               "Percent Black",
                               "Percent Hispanic", 
                               "Percent Asian"),
                selected = "Percent White"),
    
    sliderInput("range", 
                label = "Range of interest:",
                min = 0, max = 100, value = c(0, 100))
    ),
  
    mainPanel(
      h1("This is where I can put the data table containing the list of mutations"),
      p("The table will be sortable by the different mutation attributes"),
      textOutput("selected_var"),
      textOutput("min_max")
    )
  )
) #fluidPage creates a display that automatically adjusts to the dimensions of the user's browser window.


# Define server logic -----------
server <- function(input, output) {
  
  output$selected_var <- renderText({
    paste("You have selected", input$var)
  })
  
  output$min_max <- renderText({
    paste("You chose a range that goes from", input$range[1], "to", input$range[2])
  })
  
  output$dnm_table = DT::renderDataTable({
    read.delim("/data/DNMs.hg19.indelFilt.rptFilt.MAF001.RegDBv2.singletons.CADD.VEP.phastCons.SIFT.PolyPhen.fetal_brain_enhancer.DHS_fetal_brain_enh.DHS_fetal_brain_prom.turf_score.organScore.brainSp_score")
  })
}


# Run the app -----------
shinyApp(ui = ui, server = server)