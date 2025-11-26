library(shiny)
library(bslib)
library(tidyverse)
library(GEOquery)

merged_shiny = read_tsv("merged_shiny.tsv")

source("../R/99_proj_func.R")

ui <- page_sidebar(
  title = 'Boxplot showing difference in gene expression between non-control and control.',
  sidebar = sidebar(
    sliderInput(
      inputId = "genes",
      label = "Number of genes:",
      min = 2,
      max = 100,
      value = 10,
      step = 2
    ),
    radioButtons(
      "geneCat",
      "Do you want your genes to be..",
      choices = list("Randomly picked" = 'random',
                     "Positive difference" = 'max', 
                     "Negative difference" = 'min', 
                     "Positive and negative difference" = 'maxmin'),
      selected = 'random'
    ),
    radioButtons(
      "plotOrtable",
      "Do you want to show?",
      choices = list('Plot', 'Table'),
      selected = 'Plot'
    )
  ),
  conditionalPanel(
    condition = "input.plotOrtable == 'Plot'",
    plotOutput(outputId = "genePlot"),
    textOutput('descriptiveText')
  ),
  conditionalPanel(
    condition = "input.plotOrtable == 'Table'",
    tableOutput(outputId = "geneTable")
  )
)



server <- function(input, output) {
  output$genePlot <- renderPlot({
    
    merged_final_0 = gene_category(merged_final_1,input$geneCat, input$genes)
    
    ggplot(merged_final_0) + 
      geom_boxplot(mapping = aes(x = Gene_ID,
                                 y = mean_gene_expression,
                                 color = diff)) +
      scale_color_gradient2(low = "red",
                            mid = "white",
                            high = "blue",
                            midpoint = 0) +
      theme(axis.text.x=element_blank())
  })
  output$geneTable = renderTable({
    gene_category(merged_final_0,input$geneCat, input$genes)
    
  })
  output$descriptiveText <- renderText({
    "The more visible the boxplot is, the higher the difference. Blue means an positive difference, so there is a higher gene expression in non-control patients. Red means an negative difference, so there is a lower gene expression in non-control patients."
  })
}

shinyApp(ui = ui, server = server)