#function for loocv

library(pROC)
library(tidyverse)
no_tunning_weighted_averaging <- function(target_mat, cor_mat, test_idx){
  cor_test <- cor_mat[test_idx,-test_idx]
  # cor_test[cor_test < (mean(cor_test)+1*sd(cor_test)) ] <- 0
  #if any drug fail to be similar to any training drugs we just take 1 as the weighting
  if (sum(cor_test)==0) {cor_test= rep(1,nrow(cor_mat) -1)} 
  cor_mat_small = as.matrix(cor_test)
  pred1 = t(cor_mat_small) %*% as.matrix(target_mat[-test_idx,]) 
  pred1 <- as.vector(pred1)
  acc <- vector(mode = "numeric",length = length(test_idx))
  label= as.vector(target_mat[test_idx,])
  # add a minus sign to predictions because roc function automatically set controls > cases
  roc_obj <- pROC::roc(response = label, predictor = -pred1) 
  drug_auc <- pROC::auc(roc_obj )  
  return(drug_auc)
}
