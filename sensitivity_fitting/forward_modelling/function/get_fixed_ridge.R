.libPaths(c("/projappl/project_2003466/project_rpackages", .libPaths()))
library(glmnet,quietly = T)
library(tidyverse,quietly = T)
library(tidymodels,quietly = T)
options(expressions = 5e5)

glm_model <- linear_reg(penalty = tune(), mixture =0) %>% 
  set_engine("glmnet") %>% 
  set_mode("regression")

glm_param <- parameters(penalty(range = c(-4,0), 
                                trans = log10_trans()))
set.seed(0000)
glm_grid <- grid_latin_hypercube(glm_param,size = 20)

get_fixed_glmnet <- function(training_data,training_vfold,par= F){
  glm_recipe <- recipe( training_data  ) %>%
    update_role(everything()) %>%
    update_role(y, new_role = "outcome") %>% 
    step_normalize(all_predictors())
  
  glm_wflow <- workflow( ) %>% 
    add_model(glm_model) %>% 
    add_recipe(glm_recipe)
  
  glm_search <- tune_grid(
    object = glm_wflow, 
    grid = glm_grid, 
    resamples = training_vfold,
    control = tune::control_grid(verbose = F, allow_par=par)
  )
  
  glm_param_final <- select_best(glm_search, metric = "rsq")
  glm_wflow_final <- finalize_workflow(glm_wflow, glm_param_final)
  glm_wflow_final_fit <- fit(glm_wflow_final, data = training_data)
  return(glm_wflow_final_fit)
}
