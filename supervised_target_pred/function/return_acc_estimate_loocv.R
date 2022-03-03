library(furrr)
library(tidyverse)
library(tidymodels)


return_acc_estimate_cv <- function(target_tibble,predictors_tibble=NULL,cor_mat= NULL, similiarity= "spearman",acc_metric= "AUC"){
  set.seed(0000)
  source('no_tunning_weighted_averaging.R')
  source('get_target_matrix_from_long_target_tibble.R')
  
  if (("binding_score" %in% colnames(target_tibble))){  
    acc_metric= "spearman_cor"}
  
  if (is.null(cor_mat) ){
    overlapping_drug= sort(intersect(target_tibble$drug, predictors_tibble$drug))
    if (similiarity == "tahimoto"){
      predictors_mat <- predictors_tibble %>% 
        filter(drug %in% overlapping_drug) %>% 
        arrange(drug) %>% 
        select(-drug) 
      cor_mat = 1- as.matrix(vegan::vegdist(x = (predictors_mat),method = "jaccard",upper = F))
      
    } else{
      predictors_mat <- predictors_tibble %>% 
        filter(drug %in% overlapping_drug) %>% 
        arrange(drug) %>% 
        select(-drug) %>% 
        mutate_all(.funs = scale)
      cor_mat = cor(t(predictors_mat),method = similiarity) %>% replace(is.na(.), 0)
      cor_mat[cor_mat< 0 ] <- 0 
      colnames(cor_mat) <- row.names(cor_mat) <- overlapping_drug
    }
  }
  
  if (is.null(predictors_tibble)) {
    overlapping_drug= sort(intersect(target_tibble$drug, colnames(cor_mat)))
    idx <- match(overlapping_drug, colnames(cor_mat))
    cor_mat <- cor_mat[idx,idx]
  }
    
  target_mat <- get_target_mat(target_tibble = 
                                     target_tibble %>% 
                                     filter(drug %in% overlapping_drug) %>% 
                                     arrange(drug))

  fold_obj <- loo_cv(tibble(idx= 1:nrow(target_mat))) %>% 
    mutate(test_idx = map_dbl(splits, .f= function(split) {split %>% testing() %>% pull(idx)})) %>%     
    mutate(acc= future_map_dbl(
      .x = .data$test_idx,
      .f = ~no_tunning_weighted_averaging(target_mat = target_mat, cor_mat= cor_mat, test_idx= .x, acc_metric = acc_metric)))
  
  res <- fold_obj %>% 
    select(test_idx, acc) %>% 
    mutate(test_idx= unlist(test_idx)) %>% 
    arrange(test_idx) %>% 
    mutate(drug= row.names(target_mat))
  return(res)
  }


