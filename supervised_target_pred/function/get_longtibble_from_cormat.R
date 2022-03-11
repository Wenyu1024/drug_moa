# The function is used to transform pairwise similairity matrix
# into a long tible with either raw similarity or rank transformed value
library(tidyverse)
min_max_normalization <- function(x){(x-min(x))/ (max(x)-min(x)) }

get_longtibble_from_cormat <- 
  function(mat, metric="rank_percentage_pairwise" ){
    
    if (metric == "similiarity"){
      print("Similarity returned")  
      # quantile normalization can be done later if needed.
    }
    
    if (metric == "rank_percentage_pairwise"){
      print("Rank percentage pairwise returned")
      diag(mat) <- -1000 # so the drug pair AA will not impact.
      mat = apply(
        mat, 
        2, 
        FUN = function(x){rank(x,ties.method = "max")/length(x)}
        # FUN = function(x){rank(-x,ties.method = "max")}
      )
      mat = (mat + t(mat))/2
    }
    
    mat[upper.tri(mat,diag = T)] <- -1000
    mat %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "drug1") %>% 
      pivot_longer(cols = -drug1, names_to="drug2", values_to= "value") %>% 
      filter(value != -1000) %>% 
      drop_na() %>% 
      mutate(value = min_max_normalization(value)) %>% 
      mutate(value = round(value, 3))  
  }
