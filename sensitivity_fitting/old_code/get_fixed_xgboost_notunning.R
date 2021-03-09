# library(doFuture)
library(xgboost)
options(expressions = 5e5)
get_fixed_xgboost <- function(training_data,training_vfold){

  xgboost_model <- boost_tree(
    # mtry = tune(),# 0.5:0.9
    sample_size = 0.9,
    trees = 1000,
    # min_n = tune(), # 1:3
    # tree_depth = tune(), # 3:6
    learn_rate = 0.1, 
    # loss_reduction = tune(), 
    stop_iter = 5, 
  
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
  
#   # grid specification
#   xgboost_params <- xgboost_wflow %>%
#     parameters() %>%
#     update(tree_depth = tree_depth(range = c(4L, 12L))) %>% 
#     update(tree_depth = tree_depth(range = c(3, 8))) %>% 
#     update(sample_size= sample_size(range = c(0.5, 0.9))) %>% 
#     update(learning_rate = learning_rate(range = c(0.05, 0.2))) %>% 
# 
# # how does mtry different from
# 
# 
#   xgboost_grid <- grid_regular(
#     xgboost_params,
#     levels= 9)
# 
# 
#   xgboost_search <- tune_grid(
#     object = xgboost_wflow,
#     resamples = training_vfold,
#     grid = xgboost_grid,
#     # metrics = yardstick::metric_set(rmse, rsq, mae), default is rsq
#     control = tune::control_grid(verbose = T, allow_par=F)
#   )
# 
#   xgboost_param_final <- select_by_one_std_err(xgboost_search,tree_depth, metric = "rsq")
#   xgboost_param_final <- select_best(xgboost_search, metric = "rsq")
#   xgboost_wflow_final <- finalize_workflow(xgboost_wflow, xgboost_param_final)
#   xgboost_wflow_final_fit <- fit(xgboost_wflow_final, data = training_data)
  xgboost_wflow_final_fit <- fit(xgboost_wflow, data = training_data)
  return(xgboost_wflow_final_fit)
}
