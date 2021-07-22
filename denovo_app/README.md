A Shiny app is contained in a single script called `app.R`. The script lives in a directory and can be run with `runApp("directory_name")`.  
`app.R` has three components:  
- a user interface (`ui`) object  
- a server function  
- a call to the `shinyApp` function  
  
The directory `app.R` appears in will become the working directory of the Shiny app.  
Shiny will run the code at the start of `app.R`, before the `server` function, only once during the life of the app.  
Shiny will run code placed in the `server` function multiple times.  


  

