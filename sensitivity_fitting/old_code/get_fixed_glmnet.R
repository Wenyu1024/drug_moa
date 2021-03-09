# library(doFuture)
library(glmnet)
options(expressions = 5e5)
get_fixed_glmnet <- function(training_data,training_vfold){
  # predictors_columns <- str_subset(string = colnames(training_data), pattern = predictor_pattern)
  # formula <-  str_c(predictors_columns, collapse = " + ")
  # formula <- str_c("area_under_curve ~ ", formula,  collapse = " ")
  # formula <- as.formula(formula)
  
  glm_model <- linear_reg(penalty = tune(), mixture =tune()) %>% 
    set_engine("glmnet") %>% 
    set_mode("regression")
  
  # glm_recipe <- recipe( formula = area_under_curve ~ .,  data =  training_data  ) %>% 
  #   step_normalize(all_predictors())
  # Error: protect(): protection stack overflow
  # https://github.com/tidymodels/recipes/issues/548
  glm_recipe <- recipe( training_data  ) %>%
    update_role(everything()) %>%
    update_role(area_under_curve, new_role = "outcome") %>% 
    step_medianimpute(all_numeric()) %>% 
    step_normalize(all_predictors())
  
  
  glm_wflow <- workflow( ) %>% 
    add_model(glm_model) %>% 
    add_recipe(glm_recipe)
  
  glm_param <- parameters(penalty(range = c(-4,0), 
                                  trans = log10_trans()),
                          mixture(range = c(0,0.5)))
  
  glm_grid <- grid_regular(glm_param, levels= c(5,5))
  
  glm_search <- tune_grid(
    object = glm_wflow, 
    grid = glm_grid, 
    resamples = training_vfold,
    control = tune::control_grid(verbose = T, allow_par=T )
  )
  
  glm_param_final <- select_best(glm_search, metric = "rsq")
  glm_wflow_final <- finalize_workflow(glm_wflow, glm_param_final)
  glm_wflow_final_fit <- fit(glm_wflow_final, data = training_data)
  return(glm_wflow_final_fit)
  # return(glm_search)
}


# [ONE-TIME WARNING] Forked processing ('multicore') is disabled in future (>= 1.13.0) when running R from RStudio, because it is considered unstable. Because of this, plan("multicore") will fall back to plan("sequential"), and plan("multiprocess") will fall back to plan("multisession") - not plan("multicore") as in the past. For more details, how to control forked processing or not, and how to silence this warning in future R sessions, see ?future::supportsMulticore