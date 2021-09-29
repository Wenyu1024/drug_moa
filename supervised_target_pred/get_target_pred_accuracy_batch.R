#Target prediction accuracy for different sen dataset/datapoint, estimation method and target source 
library(tidyverse)
library(tidymodels)
get_target_pred_accuracy_batch <- 
  function(dataset, id_list= NULL, estimation_method= "fold_cv",target_source){
    # setwd("~/cluster_wrk/drug_moa/supervised_target_pred")
    setwd("/projappl/project_2003466/drug_moa/supervised_target_pred/")
    
    if (estimation_method== "fold_cv") {
      source('return_acc_estimate_5fold_cv.R')
    } 
    if (estimation_method== "loocv") {
      source('return_acc_estimate_loocv.R')
    } 
    if (!is.null(id_list)) {target_source <- target_source %>% filter(drug %in% id_list )}
    
    set.seed(0000)
    
    acc_maccs<- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("feature_maccs_", dataset,collapse = "") ))  ,
      similiarity = "tahimoto")
    
    acc_extended<- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("feature_extended_", dataset,collapse = "") ))  ,
      similiarity = "tahimoto")
    
    acc_consensus <- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("drug_consensus_", dataset,collapse = "") )) ) 
    
    acc_comb1<- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("feature_imp_ridge_", dataset, "_comb1" ,collapse = "") ))  ) 
    
    acc_df <- acc_extended %>% mutate(method = "structure_ECFP") %>%
      bind_rows(acc_maccs %>% mutate(method = "structure_MACCS")) %>%
      bind_rows(acc_consensus %>% mutate(method = "Consensus_exp_perturb")) %>%
      bind_rows(acc_comb1 %>% mutate(method = "senbased_genesig_comb1")) 
    return(acc_df)
  }