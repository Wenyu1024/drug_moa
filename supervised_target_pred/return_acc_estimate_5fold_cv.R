return_acc_estimate_foldcv <- function(target_tibble,predictors_tibble=NULL,cor_mat= NULL){
  if (is.null(cor_mat) ){
    overlapping_drug= intersect(target_tibble$drug, predictors_tibble$drug)
    predictors_mat <- predictors_tibble %>% 
      filter(drug %in% overlapping_drug) %>% 
      arrange(drug) %>% 
      select(-drug) %>% 
      mutate_all(.funs = scale)
    cor_mat = cor(t(predictors_mat),method = "spearman") %>% replace(is.na(.), 0)
    cor_mat[cor_mat< 0 ] <- 0 
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
    mutate(train_idx = map(splits, .f= function(split) {split %>% training() %>% pull(idx)})) %>% 
    mutate(test_idx = map(splits, .f= function(split) {split %>% testing() %>% pull(idx)})) %>%     
    mutate(acc= map_dbl(.x = .data$test_idx,
                        .f = ~no_tunning_weighted_averaging(target_mat = target_mat, cor_mat= cor_mat, test_idx= .x)))
  
  res <- fold_obj %>% select(id, id2, acc)
  return(res)
}




