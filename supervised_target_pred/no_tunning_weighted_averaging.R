library(tidyverse)
library(tidymodels)
# library(pROC)

# acc_metric can be either "AUC" or "spearman_cor"
# if acc_metric= spearman_cor, target_mat has to be continues
# if acc_metric= AUC, target_mat has to be binary.

no_tunning_weighted_averaging <- function(target_mat, cor_mat, test_idx,acc_metric= "not applicable", pred_new=F,pred_all=F){
  if (pred_new){
    cor_test <- cor_mat[test_idx,-test_idx] # 
    cor_test[apply(cor_test, 1, sum) == 0, ] <- 1
    pred1 = (cor_test) %*% as.matrix(target_mat)
    row.names(pred1) <- test_idx
    res <- pred1
    
  }else {
    test_idx <- unlist(test_idx) 
    if (length(test_idx)>1){
      cor_test <- cor_mat[test_idx,-test_idx] # 
      # cor_test[cor_test < (mean(cor_test)+1*sd(cor_test)) ] <- 0
      #if any drug fail to be similar to any training drugs we just take 1 as the weighting
      cor_test[apply(cor_test, 1, sum) == 0, ] <- 1
      
      pred1 = (cor_test) %*% as.matrix(target_mat[-test_idx,]) 
      label= target_mat[test_idx,]
      # add a minus sign to predictions because roc function automatically set controls > cases
      if (acc_metric== "AUC") {
        res <- map_dbl(.x= 1:length(test_idx), 
                       .f = function(idx){
                         df <- tibble(truth = as.factor(label[idx,]), estimate = pred1[idx,]) 
                         roc_auc_vec(df$truth,df$estimate, event_level = "second", estimator = "binary",na_rm = T )
                       })
        res <- mean(res)
      }
      
      if (acc_metric== "spearman_cor") {
        res <- map_dbl(.x= 1:length(test_idx), 
                       .f = function(idx){
                         df <- tibble(truth = label[idx,], estimate = pred1[idx,]) 
                         cor(df$truth,df$estimate, use = "complete.obs",method = "spearman" )
                       })
        res <- mean(res)
      }
      
      # if (acc_metric== "NULL") {
      #   res <- pred1
      # }
      
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
      
      
      if (acc_metric== "AUC") { 
        res <- roc_auc_vec(truth = as.factor(label), pred1,event_level = "second",estimator = "binary") 
      }
      if (acc_metric== "spearman_cor") { 
        res <- cor(label, pred1,use = "complete.obs",method = "spearman") 
      }
      
      if ((acc_metric== "not applicable")&pred_all) {
        res <- pred1
      }
    }
  }
  return(res) 
}