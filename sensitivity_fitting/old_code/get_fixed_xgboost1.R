# library(doFuture)
library(xgboost)
options(expressions = 5e5)
get_fixed_xgboost <- function(training_data,training_vfold){

  xgboost_model <- boost_tree(
    # mtry = tune(),
    trees = 1000,
    # min_n = tune(),
    tree_depth = tune(),
    learn_rate = 0.2,
    # loss_reduction = tune(),
    stop_iter = 10
  ) %>%
    set_engine("xgboost", objective = "reg:squarederror") %>% 
    set_mode("regression")
  
  xgboost_recipe <- recipe( training_data  ) %>%
    update_role(everything()) %>%
    update_role(area_under_curve, new_role = "outcome") 
    # step_normalize(all_predictors())
  
  
  
  xgboost_wflow<- workflow() %>%
    add_model(xgboost_model) %>% 
    add_recipe(xgboost_recipe)
  
  # grid specification
  xgboost_params <- xgboost_wflow %>% 
    parameters() %>% 
    # update(mtry=mtry(range = c(5L, 324L)) ) %>% 
    update(trees = tree_depth(range = c(4L, 12L))) 
  # why is there a parameter for tree when xgboost is building a unique tree
  # and the number of boosting rounds=
  
# how does mtry different from   
  
  
  xgboost_grid <- grid_regular(
    xgboost_params, 
    levels= c(5,5))
  
  # registerDoFuture()
  # plan(multiprocess)
  
  
  xgboost_search <- tune_grid(
    object = xgboost_wflow,
    resamples = training_vfold,
    grid = xgboost_grid,
    # metrics = yardstick::metric_set(rmse, rsq, mae), default is rsq
    control = tune::control_grid(verbose = F, allow_par=F)
  )
  
  xgboost_param_final <- select_by_one_std_err(xgboost_search, mtry, metric = "rsq")
  xgboost_wflow_final <- finalize_workflow(xgboost_wflow, xgboost_param_final)
  xgboost_wflow_final_fit <- fit(xgboost_wflow_final, data = training_data)
  return(xgboost_wflow_final_fit)
}
