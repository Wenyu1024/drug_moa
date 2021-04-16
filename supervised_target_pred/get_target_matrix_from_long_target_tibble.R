library(tidyverse) 
get_target_mat <- function(target_tibble){
  target_mat <- target_tibble %>% 
    distinct() %>% 
    group_by(drug, target_gene) %>%
    summarize(binding_score= median(binding_score,na.rm = T)) %>% 
    ungroup() %>% 
    pivot_wider(id_cols = drug,
                names_from= target_gene, 
                values_from=binding_score) %>% 
    arrange(drug) %>%     
    replace(is.na(.), 0) %>% 
    select(-drug)
  target_mat= as.matrix(target_mat)
  return(target_mat)
}

