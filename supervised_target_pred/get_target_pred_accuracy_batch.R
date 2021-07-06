#Target prediction accuracy for different sen dataset/datapoint, estimation method and target source 
get_target_pred_accuracy_batch <- 
  function(dataset, id_list= NULL, estimation_method= "fold_cv",target_source){
    source('~/cluster_wrk/drug_moa/supervised_target_pred/no_tunning_weighted_averaging.R')
    source('~/cluster_wrk/drug_moa/supervised_target_pred/get_target_matrix_from_long_target_tibble.R')
    
    if (estimation_method== "fold_cv") {
      source('~/cluster_wrk/drug_moa/supervised_target_pred/return_acc_estimate_5fold_cv.R')
    } 
    if (estimation_method== "loocv") {
      source('~/cluster_wrk/drug_moa/supervised_target_pred/return_acc_estimate_loocv.R')
    } 
    if (!is.null(id_list)) {target_source <- target_source %>% filter(drug %in% id_list )}
    
    set.seed(0000)
    acc_ces1 <- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("feature_imp_ridge_", dataset, "_ces1" ,collapse = "") )) )
    
    acc_ceres<- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("feature_imp_ridge_", dataset, "_ceres" ,collapse = "") )) ) 
    
    acc_demeter2<- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("feature_imp_ridge_", dataset, "_demeter2" ,collapse = "") ))  ) 
    
    acc_exp<- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("feature_imp_ridge_", dataset, "_exp" ,collapse = "") ))  )
    
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
    
    acc_comb2<- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("feature_imp_ridge_", dataset, "_comb2" ,collapse = "") ))  ) 
    
    acc_comb3<- return_acc_estimate_cv(
      target_tibble = target_source,
      predictors_tibble = eval(as.symbol(
        str_c("feature_imp_ridge_", dataset, "_comb3" ,collapse = "") ))  ) 

        
    
    acc_df <- acc_ces1 %>% mutate(method = "ces1") %>% 
      bind_rows(acc_ceres %>% mutate(method = "ceres")) %>% 
      bind_rows(acc_demeter2 %>% mutate(method = "demeter2")) %>% 
      bind_rows(acc_exp %>% mutate(method = "exp")) %>% 
      bind_rows(acc_extended %>% mutate(method = "structure_ECFP")) %>%
      bind_rows(acc_maccs %>% mutate(method = "structure_MACCS")) %>%
      bind_rows(acc_consensus %>% mutate(method = "Consensus_exp_perturb")) %>% 
      bind_rows(acc_comb1 %>% mutate(method = "senbased_genesig_comb1")) %>%
      bind_rows(acc_comb2 %>% mutate(method = "senbased_genesig_comb2")) %>% 
      bind_rows(acc_comb3 %>% mutate(method = "senbased_genesig_comb3"))
    
    return(acc_df)
  }