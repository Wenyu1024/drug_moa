library(tidyverse)
# note only the drugs with both ConSensig and ConExpsig is involved
CTRP_binary_PPIold_data <- read_csv( "~/cluster_scratch/prior/CTRP_binary_PPIold_data.csv") #217 drugs
GDSC_binary_PPIold_data <- read_csv( "~/cluster_scratch/prior/GDSC_binary_PPIold_data.csv") #n= 57
PRISM_binary_PPIold_data <- read_csv( "~/cluster_scratch/prior/PRISM_binary_PPIold_data.csv") # n= 918 drugs

set.seed(0)
sig_to_inhibitortargt_pred <- function(annotated_sig_df) {
  # for con sen sig we do relu, for con exp sig we do negative relu
  data1 <- 
    annotated_sig_df %>%
    filter(imptype== "ConSen-Sig") %>% 
    mutate(imp = case_when(imp>0 ~ imp, imp<0 ~ 0)) 
  
  data2 <- 
    annotated_sig_df %>%
    filter(imptype== "ConExp-Sig") %>% 
    mutate(imp = case_when(imp<0 ~ abs(imp)
                           # , imp>0 ~ 0
                           ))   
  
  annotate_pred_df <- bind_rows(data1,data2)
  return(annotate_pred_df)
}

sig_to_downstream_pred <- function(annotated_sig_df) {
  # for downstream direction does not matter
  annotate_pred_df <- annotated_sig_df %>% 
    mutate(imp = abs(imp)) 
  return(annotate_pred_df)
}

get_auc <- 
  function(dataset, label_level= 0){
    if (label_level ==0 ) {dataset <-  sig_to_inhibitortargt_pred(dataset)}  
    if (label_level != 0) {dataset <-  sig_to_downstream_pred(dataset)}  
    dataset %>% 
      # filter(anno_type1!="Unknown") %>% 
      mutate(anno_type2 = case_when(
        anno_type!="Target" ~anno_type,
        TRUE ~ "0")) %>%    
      mutate(anno_type2= as.numeric(anno_type2) ) %>% 
      drop_na() %>% 
      mutate(anno_binary= factor(x = (anno_type2 < (label_level+1)),levels = c(TRUE,FALSE))) %>%
      group_by(drug, imptype) %>% 
      summarise(AUC= yardstick::roc_auc_vec(truth = anno_binary,estimate = imp)) %>% 
      ungroup() 
  }

ctrp_unsupervised_auc <- get_auc(CTRP_binary_PPIold_data, label_level= 0) 
gdsc_unsupervised_auc <- get_auc(GDSC_binary_PPIold_data, label_level= 0) 
prism_unsupervised_auc <- get_auc(PRISM_binary_PPIold_data, label_level= 0) 
save(list = c("ctrp_unsupervised_auc", "gdsc_unsupervised_auc", "prism_unsupervised_auc"),
     file = "~/cluster_scratch/unsupervised_target_pred/unsupervised_target_pred_method_comparison.RData")
#note the analysis above is to compare between ConSensig and ConExpSig




