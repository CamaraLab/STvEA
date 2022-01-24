library(shiny)
source("ui.R")
source("server.R")

shinyApp(ui = build.ui(), server = server)
