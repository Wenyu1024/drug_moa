# This function takes a df and a model fitting function 
# to get an estimate of the model accuracy.

# The model, including the possible values of the 
# hyperparameters, is defined in the model function Rfile.

# The df contains both response variables (y) and 
# predictors (ends with "_predictors"), preparation of the
# df is done in the code to be submitted to the array jobs

# This function only provides a standardized frame for 
# estimating the accuracy (5-fold 3 replicate cv, 15 estimates)

estimate_performance_par <- function(df, fun_name){
  # .libPaths(c("/projappl/project_2003466/project_rpackages", .libPaths()))
  set.seed(0000)
  library(tidyverse, quietly = T)
  library(tidymodels, quietly = T)
  library(furrr, quietly = T)
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
