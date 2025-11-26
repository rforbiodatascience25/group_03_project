#All repeated code should be turned into functions and put here

#Return a sliced version of the dataframe df
gene_category <- function(df, category, genes) {
  #df = a dataframe with a "Gene_ID" column and a 'diff' column
  ##Gene_ID = a string unique for each row
  ##diff = a float
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
      arrange(desc(diff))
    
    return (df)
  }
  
  #return df with the unique Gene_ID, equal to genes, with the maximum diff
  else if (category == "max"){
    Gene_slice <- df |> 
      slice_max(diff, 
                n = genes*2)
    
    Gene_ID_slice <- Gene_slice$Gene_ID
    
    df <- df |>
      filter(Gene_ID %in% Gene_ID_slice) |> 
      arrange(desc(diff))
    
    return (df)
  }
  
  #return df with the unique Gene_ID, equal to genes, with the minimum diff
  else if (category == "min"){
    Gene_slice <- df |> 
      slice_min(diff,
                n = genes*2)
    
    Gene_ID_slice <- Gene_slice$Gene_ID
    
    df <- df |>
      filter(Gene_ID %in% Gene_ID_slice) |> 
      arrange(desc(diff))
    
    return (df)
    
  }
  
  #return df with the unique Gene_ID, equal to genes
  #half with the maximum diff and the other half with the minimum diff
  else if (category == "maxmin"){
    Gene_slice_max <- df |> 
      slice_max(diff, 
                n = genes)
    
    Gene_slice_min <- df |> 
      slice_min(diff, 
                n = genes)
    
    Gene_ID_slice <- c(Gene_slice_max$Gene_ID,
                       Gene_slice_min$Gene_ID)
    
    df <- df |>
      filter(Gene_ID %in% Gene_ID_slice) |>
      arrange(desc(diff))
    
    return (df)
  }
}