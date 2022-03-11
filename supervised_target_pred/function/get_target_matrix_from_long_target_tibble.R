library(tidyverse) 
get_target_mat <- function(target_tibble,return_format= "mat"){
  if (!("binding_score" %in% colnames(target_tibble))){  
    target_tibble <- target_tibble %>% mutate(binding_score = 1)}
  
  target_tibble <- target_tibble %>% 
    distinct() %>% 
    group_by(drug, target) %>%
    summarize(binding_score= median(binding_score,na.rm = T)) %>% 
    ungroup() %>% 
    pivot_wider(id_cols = drug,
                names_from= target, 
                values_from=binding_score) %>% 
    arrange(drug) %>%     
    replace(is.na(.), 0) 
    # select(-drug)
  target_mat= as.matrix(target_tibble[,-1])
  row.names(target_mat) <- target_tibble$drug
  target_mat <- target_mat[,colSums(target_mat)!=0]
  if (return_format== "mat"){  return(target_mat) }
  if (return_format== "widedf"){  return(target_tibble) }
}

