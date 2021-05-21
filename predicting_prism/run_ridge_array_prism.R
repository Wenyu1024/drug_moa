source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/get_fixed_glmnet.R')

estimate_performance_par_prism <- function(sen_df, predictor_df, fun_name){
  # .libPaths(c("/projappl/project_2003466/project_rpackages", .libPaths()))
  library(tidyverse, quietly = T)
  library(tidymodels, quietly = T)
  library(furrr, quietly = T)
  library(janitor,lib.loc ="/projappl/project_2003466/project_rpackages",quietly = T)
  df <- sen_df  %>% 
    select(depmap_id,auc) %>% 
    as_tibble() %>% 
    inner_join( 
      y= clean_names(predictor_df)%>% 
        rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
      by= c("depmap_id"="dep_map_id")) %>% 
    rename(y= auc) %>% 
    select(-depmap_id)
  
  plan(multicore)
  res <-  df %>% 
    nested_cv(outside=vfold_cv(v = 5,repeats = 3), 
              inside= vfold_cv(v = 5,repeats = 3)) %>% 
    mutate(train_data= future_map(.x = splits, .f = training,.options = furrr_options(seed = 0000))) %>% 
    mutate(test_data= future_map(.x = splits, .f = ~testing(.x),.options = furrr_options(seed = 0000))) %>% 
    mutate(final_fit = future_map2(
      .x = train_data,
      .y = inner_resamples,
      .f = fun_name,
      .options = furrr_options(seed = 0000))) %>%  
    mutate(pred_test = future_map2(
      .x = final_fit,
      .y = test_data,
      .f = ~predict(.x, new_data=.y),
      .options = furrr_options(seed = 0000))) %>% 
    mutate(spearman_cor = future_map2_dbl(
      .x= pred_test, 
      .y= test_data,
      .f= function(x,y){
        cor(x = x %>% pull(.pred), 
            y = y %>% pull(y), 
            method = "spearman")},
      .options = furrr_options(seed = 0000))) %>% 
    pull(spearman_cor)
  plan(sequential)
  return(res)
}
load("/scratch/project_2003466/prism/prism_sensitivityfitting_serverinput.RData")

library(tidyverse, quietly = T)
library(furrr, quietly = T)
print("start")

args <- as.numeric( commandArgs( TRUE ) )

job_id <- args[1]
print(job_id)

sen_df <- gdsc_df %>% 
  slice(job_id) %>% 
  pull(sensitivity) 
sen_df <- sen_df[[1]]
ces1_perf= estimate_performance_par_prism(sen_df = sen_df, predictor_df = ces1_478,fun_name = get_fixed_glmnet)
print("ces1 success")

file_name = paste0( '/scratch/project_2003466/forward_modelling/prism_spearman/drug_prism_' ,job_id,".RData" )
# save(list = c('ces1_perf'), file = file_name)

save(list = c('ces1_perf'), file = file_name)
# save(list = c('ces1_perf'), file = file_name)

print("all success")
