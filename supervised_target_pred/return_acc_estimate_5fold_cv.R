return_acc_estimate_cv <- function(target_tibble,predictors_tibble=NULL,cor_mat= NULL, similiarity= "spearman",acc_metric= "AUC"){
  set.seed(0000)
  source('/projappl/project_2003466/drug_moa/supervised_target_pred/no_tunning_weighted_averaging.R')
  source('/projappl/project_2003466/drug_moa/supervised_target_pred/get_target_matrix_from_long_target_tibble.R')
  
  if (("binding_score" %in% colnames(target_tibble))){  
    acc_metric= "spearman_cor"}
  
  if (is.null(cor_mat) ){
    overlapping_drug= intersect(target_tibble$drug, predictors_tibble$drug)

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
  
  fold_obj <- vfold_cv(data = tibble(idx= 1:nrow(target_mat)),v = 5,
                       repeats = 3) %>%
    mutate(train_idx = future_map(splits, .f= function(split) {split %>% training() %>% pull(idx)})) %>%
    mutate(test_idx = future_map(splits, .f= function(split) {split %>% testing() %>% pull(idx)})) %>%
    mutate(acc= future_map_dbl(.x = .data$test_idx,
                        .f = ~no_tunning_weighted_averaging(target_mat = target_mat, cor_mat= cor_mat, test_idx= .x, acc_metric = acc_metric)))

  # w <- no_tunning_weighted_averaging(target_mat = target_mat, cor_mat= cor_mat, test_idx= 1:20)

  res <- fold_obj %>% select(id, id2, acc)
  return(res)
}




