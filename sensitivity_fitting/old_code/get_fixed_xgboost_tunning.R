# library(doFuture)
library(xgboost)
options(expressions = 5e5)

xgboost_model <- boost_tree(

  sample_size = 0.9,
  trees = 300,
  tree_depth = tune(),
  min_n = tune(),
  learn_rate = 0.1,
  stop_iter = 3 
) %>%
  set_engine("xgboost", objective = "reg:squarederror") %>% 
  set_mode("regression")



# grid specifications
xgboost_params <- xgboost_model %>%
  parameters() %>%
  update(tree_depth = tree_depth(range = c(3, 6))) %>% 
  update(min_n = min_n(range = c(2, 4))) 
set.seed(1234)
xgboost_grid <- grid_latin_hypercube(
  xgboost_params,
  size = 3
)
# xgboost_grid <- xgboost_grid[6:10,]
get_fixed_xgboost <- function(training_data,training_vfold){

  xgboost_recipe <- recipe( training_data  ) %>%
    update_role(everything()) %>%
    update_role(area_under_curve, new_role = "outcome") 

  xgboost_wflow<- workflow() %>%
    add_model(xgboost_model) %>% 
    add_recipe(xgboost_recipe)
  
  xgboost_search <- tune_grid(
    object = xgboost_wflow,
    resamples = training_vfold,
    grid = xgboost_grid,
    control = tune::control_grid(verbose = T, allow_par=F)
  )

  xgboost_param_final <- select_by_one_std_err(xgboost_search,tree_depth, metric = "rsq")
  xgboost_param_final <- select_best(xgboost_search, metric = "rsq")
  xgboost_wflow_final <- finalize_workflow(xgboost_wflow, xgboost_param_final)
  xgboost_wflow_final_fit <- fit(xgboost_wflow_final, data = training_data)

  return(xgboost_wflow_final_fit)
}
