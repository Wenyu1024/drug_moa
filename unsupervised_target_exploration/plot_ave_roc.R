library(tidyverse)
# library(tidymodels)
library(ROCR)
plot_ave_roc <- function(annotated_pred_df, title, addlegend= T,label_level= 0){
  sig_to_inhibitortargt_pred <- function(annotated_sig_df) {
    # for con sen sig we do relu, for con exp sig we do negative relu
    data1 <- 
      annotated_sig_df %>%
      filter(imptype== "ess-sig") %>% 
      mutate(imp = case_when(imp>0 ~ imp, imp<0 ~ 0)) 
    
    data2 <- 
      annotated_sig_df %>%
      filter(imptype== "exp-sig") %>% 
      mutate(imp = case_when(imp<0 ~ abs(imp), imp>0 ~ 0))   
    
    annotate_pred_df <- bind_rows(data1,data2)
    return(annotate_pred_df)
  }
  
  sig_to_downstream_pred <- function(annotated_sig_df) {
    # for downstream direction does not matter
    annotate_pred_df <- annotated_sig_df %>% 
      mutate(imp = abs(imp)) 
    return(annotate_pred_df)
  }
  
  #prepare prediction
  if (label_level ==0 ) {annotated_pred_df <-  sig_to_inhibitortargt_pred(annotated_pred_df)}  
  if (label_level != 0) {annotated_pred_df <-  sig_to_downstream_pred(annotated_pred_df)}  
  # prepare label
  annotated_pred_df <- annotated_pred_df %>%
    # filter(anno_type1!="Unknown") %>%
    mutate(anno_type2 = case_when(
      anno_type!="Target" ~anno_type,
      TRUE ~ "0")) %>%    
    mutate(anno_type2= as.numeric(anno_type2) ) %>% 
    drop_na() %>% 
    mutate(anno_binary= factor(x = (anno_type2 < (label_level+1)),levels = c(TRUE,FALSE)))
  
  tmp0 <- annotated_pred_df %>%
    group_by(drug, imptype) %>% 
    summarise(AUC= yardstick::roc_auc_vec(truth = anno_binary,estimate = imp)) %>%
    ungroup()
  
  tmp <- annotated_pred_df %>%
    mutate(labels=anno_binary,predictions= imp) %>% 
    select(drug, imptype, labels,predictions) %>% 
    nest_by(drug, imptype)  %>% 
    inner_join(tmp0) %>% drop_na() # filter out drugs where AUC calculation is not available.
  
  ## here I used transpose from purrr package to change the index order of this multi-level list
  ROCR.consensig <- tmp %>% filter(imptype== "ess-sig") %>% pull(data) %>% transpose()
  ROCR.conexpsig <- tmp %>% filter(imptype== "exp-sig") %>% pull(data) %>% transpose()
  
  pred <- prediction(ROCR.consensig$predictions, ROCR.consensig$labels)
  perf <- performance(pred,'tpr','fpr')
  pred1 <- prediction(ROCR.conexpsig$predictions, ROCR.conexpsig$labels)
  perf1 <- performance(pred1,'tpr','fpr')
  plot(perf, avg="vertical", lwd=3, col="red",spread.estimate="stderror",plotCI.lwd=2,xlab= "FPR", ylab= "TPR",main=NULL)
  plot(perf1, avg="vertical", lwd=3, col="blue",spread.estimate="stderror",plotCI.lwd=2,add=T,xlab= "FPR", ylab= "TPR",main=NULL)
  if(addlegend== T){legend("bottomright", legend = c("ess-sig", "exp-sig"), lty = 1, lwd = 2, col = c("red", "blue"))}
}