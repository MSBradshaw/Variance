library(shiny)
library(ggplot2)
library(tibble)
library(readr)
library(dplyr)
library(plotrix)

ui <- fluidPage(
  fileInput("file1", "Choose CSV File",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv",".pileup"), multiple = TRUE
  ),
  selectInput("image_type", "Display Type",
              c('normal','variance','consensus'), selected = 2
              ),
  numericInput("front", "Trim from front:", 600, min = 1, max = 10000),
  numericInput("end", "Trim from back:", 600, min = 1, max = 10000),
      htmlOutput("image_area")
  
)
server <- function(input, output) {
  source('parsepileup.R')
  
  output$image_area <- renderUI({
    if( is.na(input$file1) ){
      return(NULL)
    }
    names <- c('8664_2','8664_3','8664_4','8664_5','8664_6')
    names <- c('8663_A','8663_O','8663_T','8664_2','8664_3','8664_4','8664_5','8664_6','8665_N','8668_A','8668_B','8668_G','8678','8795_TA','8796','8796_TA','8797','8800','8801','8802_TC','8803_GA','8805','8806','8810','8935_B','8935_P')
   # names <- c('MCM7','RPB1','RPB2')
    paths <- c('final_dirs/8664_2/pileup.pileup',
               'final_dirs/8664_3/pileup.pileup',
               'final_dirs/8664_4/pileup.pileup',
               'final_dirs/8664_5/pileup.pileup',
               'final_dirs/8664_6/pileup.pileup')
    #input$front,input$end
    return( HTML( build_image( input$file1$datapath,input$image_type,'',input$front,input$end,names ) ) )
    })
}
shinyApp(ui = ui, server = server)
