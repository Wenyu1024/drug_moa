#this code combines and update the effort in run_glm.R
.libPaths(c("/projappl/project_2003466/project_rpackages", .libPaths()))
library(glmnet,quietly = T)
library(tidyverse,quietly = T)
library(tidymodels,quietly = T)
options(expressions = 5e5)

glm_model <- linear_reg(penalty = tune(), mixture =tune()) %>% 
  set_engine("glmnet") %>% 
  set_mode("regression")



glm_param <- parameters(
  penalty(range = c(-4,0), trans = log10_trans()),
                        mixture(range = c(0,1))
                        )

glm_grid <- grid_latin_hypercube(glm_param,size = 20)

# training_data= FK866 %>% select_at(10:10578) %>% slice(1:200) %>% rename(y=ic50)
get_fixed_glmnet <- function(training_data,training_vfold,par= F){
  cor_res = cor(x = training_data$y, y = training_data %>% select(-y), use="complete.obs")
  predictors_selected = colnames(cor_res)[which(abs(cor_res)  > quantile(abs(cor_res),0.9)) ]
  
  glm_recipe <- recipe( training_data  ) %>%
    # update_role(everything()) %>%
    update_role(all_of(predictors_selected), new_role = "predictor") %>%
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
