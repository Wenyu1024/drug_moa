# LOOCV
# predictors_tibble = feature_imp_ridge_ctrp_ces1 %>% rename(drug=broad_cpd_id)
# target_tibble= ctrp_target_tibble %>%
#   ungroup() %>%
#   select(broad_cpd_id, gene_symbol_of_protein_target) %>%
#   rename(drug=broad_cpd_id, target_gene=gene_symbol_of_protein_target) %>% 
#   mutate(binding_score= 1)


return_acc_estimate_loocv <- function(target_tibble,predictors_tibble=NULL,cor_mat= NULL){
  if (is.null(cor_mat) ){
    overlapping_drug= intersect(target_tibble$drug, predictors_tibble$drug)
    predictors_mat <- predictors_tibble %>% 
      filter(drug %in% overlapping_drug) %>% 
      arrange(drug) %>% 
      select(-drug) %>% 
      mutate_all(.funs = scale)
    cor_mat = cor(t(predictors_mat)) %>% replace(is.na(.), 0)
    cor_mat[cor_mat< 0 ] <- 0 
  }
  
  if (is.null(predictors_tibble)) {
    overlapping_drug= intersect(target_tibble$drug, colnames(cor_mat))
  }
    
  target_mat <- get_target_mat(target_tibble = 
                                     target_tibble %>% 
                                     filter(drug %in% overlapping_drug) %>% 
                                     arrange(drug))

  fold_obj <- loo_cv(tibble(idx= 1:nrow(target_mat))) %>% 
    mutate(train_idx = map(splits, .f= function(split) {split %>% training() %>% pull(idx)})) %>% 
    mutate(test_idx = map(splits, .f= function(split) {split %>% testing() %>% pull(idx)})) %>%     
    mutate(acc= map_dbl(.x = .data$test_idx,
                        .f = ~no_tunning_weighted_averaging(target_mat = target_mat, cor_mat= cor_mat, test_idx= .x)))
  
  res <- fold_obj %>% select(test_idx, acc)
  return(res)
  }

