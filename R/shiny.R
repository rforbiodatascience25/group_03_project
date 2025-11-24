library(shiny)
library(bslib)
library(tidyverse)
library(GEOquery)

#
gse <- getGEO("GSE54514", GSEMatrix = TRUE)

# If multiple platforms, take the first one
gse <- gse[[1]]

# Inspect
exprs_data <- exprs(gse)
pheno_data <- pData(gse)
feature_data <- fData(gse)

#

#Since its a matrix we make it df
df_ed <- as.data.frame(exprs_data)

#And we keep the rownames 
df_ed$gene_ID <- rownames(df_ed)
df_ed <- df_ed[, c("gene_ID", colnames(df_ed)[colnames(df_ed) != "gene_ID"])]

df_ed_tibble = tibble(df_ed)

df_ed_long <- df_ed_tibble |> pivot_longer(cols = GSM1317896:GSM1318058, 
                                           names_to = 'Gene_name', 
                                           values_to = 'Gene_expression')

#

pheno_data <- pData(gse)
pheno_data_tibble = tibble(pheno_data)

pheno_data_control <- pheno_data_tibble |> 
  mutate(control = str_detect(title, '(Control)'))  |> 
  mutate(control = factor(control, levels = c(TRUE, FALSE)))  |> 
  group_by(control)  |> 
  select(geo_accession, control)

#

#combine  df_ed_long and pheno_data_control
df_gene <- left_join(x = df_ed_long,
                     y = pheno_data_control, 
                     by = join_by(Gene_name == geo_accession))

#Get list of distinct gene_ID
gene_ID_distinct <- df_gene |> distinct(gene_ID)

#Add control and gene_id as factors and groups
df_gene <- df_gene |> 
  mutate(control = factor(control, levels = c(TRUE, FALSE)),
         gene_ID = factor(gene_ID,gene_ID_distinct$gene_ID))  |> 
  group_by(control,gene_ID)

#Get the mean gene expression for each gene
df_gene_mean <- summarize(df_gene, mean = mean(Gene_expression))

#Pivot wide, turning control into two columns with values from mean
df_gene_mean = pivot_wider(df_gene_mean, names_from = control,
                           values_from = mean)

#Add column with difference in gene expression between control and non-control
#For each gene
df_gene_mean = mutate(df_gene_mean,
                      diff = `FALSE`-`TRUE`)

#Join df_gene and df_gene_mean
df_gene <- left_join(x = df_gene,
                     y = df_gene_mean,
                     by = join_by(gene_ID))

df_gene = ungroup(df_gene)

df_gene = select(df_gene, gene_ID, `TRUE`, `FALSE`, diff)
df_gene = pivot_longer(df_gene, col = c(`TRUE`, `FALSE`), names_to = 'control', values_to = 'mean_gene_expression')
df_gene = distinct(df_gene)

df_gene_0 = df_gene

gene_category <- function(df_gene_0, category, genes) {
  if (category == 'random'){
    Gene_slice <- df_gene_0 |> distinct(gene_ID) |> slice_sample(n = genes)
    gene_ID_slice <- Gene_slice$gene_ID
    df_gene <- df_gene_0 |>
      filter(gene_ID %in% gene_ID_slice)
    df_gene <- df_gene |> arrange(desc(diff))
    return (df_gene)
  }
  else if (category == 'max'){
    Gene_slice <- df_gene_0 |> slice_max(diff, n = genes*2)
    gene_ID_slice <- Gene_slice$gene_ID
    df_gene <- df_gene_0 |>
      filter(gene_ID %in% gene_ID_slice)
    df_gene <- df_gene |> arrange(desc(diff))
    return (df_gene)
    
  }
  else if (category == 'min'){
    Gene_slice <- df_gene_0 |> slice_min(diff, n = genes*2)
    gene_ID_slice <- Gene_slice$gene_ID
    df_gene <- df_gene_0 |>
      filter(gene_ID %in% gene_ID_slice)
    df_gene <- df_gene |> arrange(desc(diff))
    return (df_gene)
    
  }
  else if (category == 'maxmim'){
    Gene_slice_max <- df_gene_0 |> slice_max(diff, n = genes)
    Gene_slice_min <- df_gene_0 |> slice_min(diff, n = genes)
    c(Gene_slice_max$gene_ID,Gene_slice_min$gene_ID)
    gene_ID_slice <- c(Gene_slice_max$gene_ID,Gene_slice_min$gene_ID)
    df_gene <- df_gene_0 |>
      filter(gene_ID %in% gene_ID_slice)
    df_gene <- df_gene |> arrange(desc(diff))
    return (df_gene)
    
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
                     "Positive and negative difference" = 'maxmim'),
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
    
    df_gene = gene_category(df_gene_0,input$geneCat, input$genes)
    
    ggplot(df_gene) + 
      geom_boxplot(mapping = aes(x = gene_ID,
                                 y = mean_gene_expression,
                                 color = diff)) +
      scale_color_gradient2(low = "red",
                            mid = "white",
                            high = "blue",
                            midpoint = 0) +
      theme(axis.text.x=element_blank())
  })
  output$geneTable = renderTable({
    gene_category(df_gene_0,input$geneCat, input$genes)
    
  })
  output$descriptiveText <- renderText({
    "The more visible the boxplot is, the higher the difference. Blue means an positive difference, so there is a higher gene expression in non-control patients. Red means an negative difference, so there is a lower gene expression in non-control patients."
  })
  
  
}

shinyApp(ui = ui, server = server)

