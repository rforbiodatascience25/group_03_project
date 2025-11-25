library(shiny)
library(bslib)
library(tidyverse)
library(GEOquery)

merged_final = read_tsv('merged_final.tsv')

merged_final_0 = merged_final

#

merged_final_0 = select(merged_final_0, Gene_ID, sample, expression)

#
gse <- getGEO("GSE54514", GSEMatrix = TRUE)

# If multiple platforms, take the first one
gse <- gse[[1]]
pheno_data <- pData(gse)
pheno_data_tibble = tibble(pheno_data)

pheno_data_control <- pheno_data_tibble |> 
  mutate(control = str_detect(title, '(Control)'))  |> 
  mutate(control = factor(control, levels = c(TRUE, FALSE)))  |> 
  group_by(control)  |> 
  select(geo_accession, control)

merged_final_0 <- left_join(x = merged_final_0,
                            y = pheno_data_control, 
                            by = join_by(sample == geo_accession))

#Get list of distinct gene_ID
Gene_ID_distinct <- merged_final_0 |> distinct(Gene_ID)

#Add control and gene_id as factors and groups
merged_final_0 <- merged_final_0 |> 
  mutate(control = factor(control, levels = c(TRUE, FALSE)),
         Gene_ID = factor(Gene_ID,Gene_ID_distinct$Gene_ID))  |> 
  group_by(control,Gene_ID)

#Get the mean gene expression for each gene
merged_final_mean <- summarize(merged_final_0, mean = mean(expression))

#Pivot wide, turning control into two columns with values from mean
merged_final_mean = pivot_wider(merged_final_mean, names_from = control,
                                values_from = mean)

#Add column with difference in gene expression between control and non-control
#For each gene
merged_final_mean = mutate(merged_final_mean,
                           diff = `FALSE`-`TRUE`)

#Join df_gene and df_gene_mean
merged_final_0 <- left_join(x = merged_final_0,
                            y = merged_final_mean,
                            by = join_by(Gene_ID))

merged_final_0 = ungroup(merged_final_0)

merged_final_0 = select(merged_final_0, Gene_ID, `TRUE`, `FALSE`, diff)
merged_final_0 = pivot_longer(merged_final_0, col = c(`TRUE`, `FALSE`), names_to = 'control', values_to = 'mean_gene_expression')
merged_final_0 = distinct(merged_final_0)

merged_final_1 = merged_final_0
merged_final_1

gene_category <- function(merged_final_1, category, genes) {
  if (category == 'random'){
    Gene_slice <- merged_final_1 |> distinct(Gene_ID) |> slice_sample(n = genes)
    Gene_ID_slice <- Gene_slice$Gene_ID
    merged_final_0 <- merged_final_1 |>
      filter(Gene_ID %in% Gene_ID_slice)
    merged_final_0 <- merged_final_0 |> arrange(desc(diff))
    return (merged_final_0)
  }
  else if (category == 'max'){
    Gene_slice <- merged_final_1 |> slice_max(diff, n = genes*2)
    Gene_ID_slice <- Gene_slice$Gene_ID
    merged_final_0 <- merged_final_1 |>
      filter(Gene_ID %in% Gene_ID_slice)
    merged_final_0 <- merged_final_0 |> arrange(desc(diff))
    return (merged_final_0)
    
  }
  else if (category == 'min'){
    Gene_slice <- merged_final_1 |> slice_min(diff, n = genes*2)
    print(Gene_slice)
    Gene_ID_slice <- Gene_slice$Gene_ID
    merged_final_0 <- merged_final_1 |>
      filter(Gene_ID %in% Gene_ID_slice)
    merged_final_0 <- merged_final_0 |> arrange(desc(diff))
    return (merged_final_0)
    
  }
  else if (category == 'maxmin'){
    Gene_slice_max <- merged_final_1 |> slice_max(diff, n = genes)
    Gene_slice_min <- merged_final_1 |> slice_min(diff, n = genes)
    Gene_ID_slice <- c(Gene_slice_max$Gene_ID,Gene_slice_min$Gene_ID)
    merged_final_0 <- merged_final_1 |>
      filter(Gene_ID %in% Gene_ID_slice)
    merged_final_0 <- merged_final_0 |> arrange(desc(diff))
    return (merged_final_0)
    
  }
}

# Define UI for app that draws a histogram ----
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

# Define server logic required to draw a histogram ----
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

