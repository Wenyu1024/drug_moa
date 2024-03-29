# This function only provides a standardized frame for 
# estimating the accuracy (5-fold 3 replicate cv, 15 estimates)

# Input is a tibble object(df) as well as model function

# The df contains both response variables (y) and 
# predictors (ends with "_predictors"), preparation of the
# df is done in the code to be submitted to the array jobs

# The model definition is coded in the model function Rfile.

.libPaths(c("/projappl/project_2003466/project_rpackages", .libPaths()))
library(tidyverse, quietly = T)
library(tidymodels, quietly = T)
library(furrr, quietly = T)

estimate_performance_par <- function(df, fun_name){
  set.seed(0000)
  res <-  df %>% 
    nested_cv(outside=vfold_cv(v = 5,repeats = 3), 
              inside= bootstraps(times = 10)) %>% 
    mutate(train_data= future_map(.x = splits, .f = training,.options = furrr_options(seed = 0000))) %>% 
    mutate(test_data= future_map(.x = splits, .f = ~testing(.x),.options = furrr_options(seed = 0000))) %>% 
    mutate(final_fit = future_map2(
      .x = train_data,
      .y = inner_resamples,
      .f = fun_name,
      .options = furrr_options(seed = 0000))) %>%  
    mutate(final_res = future_map2(
      .x = final_fit,
      .y = test_data,
      .f = function(model= .x,new_data=.y){
        predict(model, new_data= new_data) %>% 
          mutate(label= new_data$y)
      },
      .options = furrr_options(seed = 0000)
      )) %>%
    select(id, id2, final_res)
    # mutate(spearman_cor = future_map2_dbl(
    #   .x= pred_test, 
    #   .y= test_data,
    #   .f= function(x,y){
    #     cor(x = x %>% pull(.pred), 
    #         y = y %>% pull(y),
    #         method = "spearman")},
    #   .options = furrr_options(seed = 0000))) %>% 
    # mutate(final_res = future_map2(
    #   .x = final_fit,
    #   .y = test_data,
    #   .f = predict(x=.x , 
    #                 new_data= .y) %>% 
    #     select(y,.pred),
    #   .options = furrr_options(seed = 0000))) %>% 
  return(res)
}
