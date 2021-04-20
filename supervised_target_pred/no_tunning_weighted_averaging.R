library(tidyverse)
library(pROC)

no_tunning_weighted_averaging <- function(target_mat, cor_mat, test_idx){
  test_idx <- unlist(test_idx)
  if (length(test_idx)>1){
      cor_test <- cor_mat[test_idx,-test_idx] # 
      # cor_test[cor_test < (mean(cor_test)+1*sd(cor_test)) ] <- 0
      #if any drug fail to be similar to any training drugs we just take 1 as the weighting
      cor_test[apply(cor_test, 1, sum) == 0, ] <- 1
      
      pred1 = (cor_test) %*% as.matrix(target_mat[-test_idx,]) 
      label= target_mat[test_idx,]
      # add a minus sign to predictions because roc function automatically set controls > cases
      res <- map_dbl(.x= 1:length(test_idx), 
                     .f = function(idx){
                       roc_auc_vec(truth = as.factor(label[idx,]), estimate = pred1[idx,],event_level = "second",estimator = "binary") 
                     })
      res <- mean(res)
      return(res)
  }
  
  
  if (length(test_idx)==1){ 
    cor_test <- cor_mat[test_idx,-test_idx]
    # cor_test[cor_test < (mean(cor_test)+1*sd(cor_test)) ] <- 0
    #if any drug fail to be similar to any training drugs we just take 1 as the weighting
    if (sum(cor_test)==0) {cor_test= rep(1,nrow(cor_mat) -1)} 
    cor_mat_small = as.matrix(cor_test)
    pred1 = t(cor_mat_small) %*% as.matrix(target_mat[-test_idx,]) 
    pred1 <- as.vector(pred1)
    acc <- vector(mode = "numeric",length = length(test_idx))
    label= as.vector(target_mat[test_idx,])
    drug_auc <- roc_auc_vec(truth = as.factor(label), pred1,event_level = "second",estimator = "binary") 
    return(drug_auc)
    # return(pred1)
    # return(label)
    # res= tibble(response = label, predictor = pred1)
    # return(res)
  }
}