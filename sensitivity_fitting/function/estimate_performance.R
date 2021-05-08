# note it is very important to library inside function that need parallel, 
# with lib path specified to your own lib.loc
# especially in a multinode job when it can be tricky for machines to correctly
# 

# this function is revised to parallel on the resamples level to cope with array jobs
estimate_performance_par <- function(sen_df, predictor_df, fun_name){
  # .libPaths(c("/projappl/project_2003466/project_rpackages", .libPaths()))
  library(tidyverse, quietly = T)
  library(tidymodels, quietly = T)
  library(furrr, quietly = T)
  library(janitor,lib.loc ="/projappl/project_2003466/project_rpackages",quietly = T)
  df <- sen_df  %>% 
    as_tibble() %>% 
    inner_join( 
      y= clean_names(predictor_df)%>% 
        rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
      by= c("DepMap_ID"="dep_map_id")) %>% 
    select(-master_ccl_id, -apparent_ec50_umol, -DepMap_ID) %>% 
    rename(y= area_under_curve)
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


estimate_performance_par_gdsc <- function(sen_df, predictor_df, fun_name){
  # .libPaths(c("/projappl/project_2003466/project_rpackages", .libPaths()))
  library(tidyverse, quietly = T)
  library(tidymodels, quietly = T)
  library(furrr, quietly = T)
  library(janitor,lib.loc ="/projappl/project_2003466/project_rpackages",quietly = T)
  df <- sen_df  %>% 
    as_tibble() %>% 
    inner_join( 
      y= clean_names(predictor_df)%>% 
        rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
      by= c("DepMap_ID"="dep_map_id")) %>% 
    select(-LN_IC50 , -DepMap_ID) %>% 
    rename(y= AUC)
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
