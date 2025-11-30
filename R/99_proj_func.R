#All repeated code should be turned into functions and put here

#Import function gene_category
gene_category <- function(df, category, genes) {
  #df = a dataframe with a "Gene_ID" column and a 'difference_mean_gene_expression' column
  ##Gene_ID = a string unique for each row
  ##difference_mean_gene_expression = a float
  #category = one of the four categories, 'random', 'max', 'min' or 'maxmin'
  #genes = a int
  
  #return df with random unique Gene_ID equal to genes
  if (category == "random"){
    Gene_slice <- df |>
      distinct(Gene_ID) |> 
      slice_sample(n = genes)
    
    Gene_ID_slice <- Gene_slice$Gene_ID
    
    df <- df |> 
      filter(Gene_ID %in% Gene_ID_slice) |> 
      arrange(desc(difference_mean_gene_expression))
    
    return (df)
  }
  
  #return df with the unique Gene_ID, equal to genes, 
  #with the maximum difference in mean gene expression
  else if (category == "max"){
    Gene_slice <- df |> 
      slice_max(difference_mean_gene_expression, 
                n = genes*2)
    
    Gene_ID_slice <- Gene_slice$Gene_ID
    
    df <- df |>
      filter(Gene_ID %in% Gene_ID_slice) |> 
      arrange(desc(difference_mean_gene_expression))
    
    return (df)
  }
  
  #return df with the unique Gene_ID, equal to genes, 
  #with the minimum difference in mean gene expression
  else if (category == "min"){
    Gene_slice <- df |> 
      slice_min(difference_mean_gene_expression,
                n = genes*2)
    
    Gene_ID_slice <- Gene_slice$Gene_ID
    
    df <- df |>
      filter(Gene_ID %in% Gene_ID_slice) |> 
      arrange(desc(difference_mean_gene_expression))
    
    return (df)
    
  }
  
  #return df with the unique Gene_ID, equal to genes, 
  #with almost no difference in mean gene expression
  else if (category == "mid"){
    df_mid <- df |>
      mutate(difference_mean_gene_expression = abs(difference_mean_gene_expression))
    
    Gene_slice <- df_mid |> 
      slice_min(difference_mean_gene_expression,
                n = genes*2)
    
    Gene_ID_slice <- Gene_slice$Gene_ID
    
    df <- df |>
      filter(Gene_ID %in% Gene_ID_slice) |> 
      arrange(desc(difference_mean_gene_expression))
    
    return (df)
  }
}