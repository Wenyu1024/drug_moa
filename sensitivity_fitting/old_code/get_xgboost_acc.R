library(furrr)
library(tidyverse)
library(tidymodels)
# source('~/drug_target/analysis/get_fixed_xgboost_no_tunning.R')
# sen_df refers to the cellular sensitivity data of a drug


estimate_performance_xgboost <- function(sen_df, ess){
  df <- sen_df  %>%   
    as_tibble() %>% 
    inner_join( 
      y= janitor::clean_names(ess)%>% 
        rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
      by= c("DepMap_ID"="dep_map_id")) %>% 
    select(-master_ccl_id, -apparent_ec50_umol, -DepMap_ID) 
  
  
  df1 <-  df %>% 
    nested_cv(outside=vfold_cv(v = 5), inside= vfold_cv(v = 5)) %>% 
    mutate(train_data= map(.x = splits, .f = training)) %>% 
    mutate(test_data= map(.x = splits, .f = ~testing(.x))) %>% 
    mutate(xg_wflow_final_fit = 
             map2(.x = train_data,.y = inner_resamples,
                  .f = ~get_fixed_xgboost(training_data = .x, training_vfold = .y))) %>%  
    mutate(pred_test = map2(.x = xg_wflow_final_fit,
                            .y = test_data,
                            .f = ~predict(.x, new_data=.y))) %>% 
    mutate(spearman_cor =  map2_dbl(.x= pred_test, .y= test_data,
                                    .f= function(x,y){
                                      cor(x = x %>% pull(.pred), y = y %>% pull(area_under_curve), method = "spearman")}))
  res= df1 %>% pull(spearman_cor)
  # return(df1)
  return(res)
  
}




# load("~/drug_target/data/ctrpv2_data.RData")
# load("~/drug_target/data/new_ess_modeling.RData")
# load("/home/cloud-user/drug_target/analysis/xg_input.RData")

data <- ctrpv2_data %>%
  ungroup() 
  # dplyr::slice(1:72)
plan(multicore,workers= 36)
# plan(sequential)
a <- Sys.time()
res2 <- data %>% 
  mutate(ces1_perf= future_map(.x = sensitivity,
                               .f = ~estimate_performance_xgboost(sen_df = .x,ess = ces1),
                               .options = furrr_options(seed = 0000))) 
  # mutate(ces2_perf= future_map(sensitivity,
  #                              .f = ~estimate_performance_xgboost(sen_df = .x,ess = ces2),
  #                              .options = furrr_options(seed = 0000))) %>%
  # mutate(ceres_perf= future_map(sensitivity,
  #                               .f = ~estimate_performance_xgboost(sen_df = .x,ess = ceres),
  #                               .options = furrr_options(seed = 0000))) %>%
  # mutate(demeter2_perf= future_map(sensitivity,
  #                                  .f = ~estimate_performance_xgboost(sen_df = .x,ess = demeter2),
  #                                  .options = furrr_options(seed = 0000)))
  # mutate(demeter2_full_perf= future_map(sensitivity,
  #                              .f = ~estimate_performance_xgboost(sen_df = .x,ess = demeter2_imputed),
  #                              .options = furrr_options(seed = 0000))) %>%
  #   mutate(ceres_full_perf= future_map(sensitivity,
  #                                 .f = ~estimate_performance_xgboost(sen_df = .x,ess = ceres_imputed),
  #                                 .options = furrr_options(seed = 0000))) %>%
  #   mutate(demeter2_filtered_perf= future_map(sensitivity,
  #                                    .f = ~estimate_performance_xgboost(sen_df = .x,ess = demeter2_filtered_imputed),
  #                                    .options = furrr_options(seed = 0000)))

  

# tmp1 <- ctrpv2_data$sensitivity[[1]]
# tmp <- estimate_performance_glmnet(sen_df = tmp1, ess = ces1[,1:100])
b <- Sys.time() 
b-a 
# save.image("~/drug_target/analysis/all_drug_forward_modeling2.RData")



# task

# it confuse me how to early stop without specifying test set. Is it possible to 
# early stop based on the training set? 
# seems possible according to some descriptions but have no idea.

#lets throw a xgboost tunning into the VM parallel1 and then check if we can explore
# how early stopping is working on one drug on parallel2.

