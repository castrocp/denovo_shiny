library(shiny)

# Define UI -----------
ui <- fluidPage( 
  titlePanel("A Shiny app for viewing de novo mutations"),
  
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
      p("The table will be sortable by the different mutation attributes")
    )
  )
) #fluidPage creates a display that automatically adjusts to the dimensions of the user's browser window.


# Define server logic -----------
server <- function(input, output) {
  
}


# Run the app -----------
shinyApp(ui = ui, server = server)