.libPaths(c("/projappl/project_2003466/project_rpackages", .libPaths()))
library(glmnet,quietly = T)
library(tidyverse,quietly = T)
library(tidymodels,quietly = T)


return_feature_importance <- function(fit){
  fit %>% 
    tidy() %>% 
    filter(term!="(Intercept)") %>% 
    arrange(term) %>% 
    pull(estimate)
}

get_feature_imp <- function(df, fun_name){
  set.seed(0000)
  df1 <-  df %>% vfold_cv(v = 5,repeats = 3)
  glm_wflow_final_fit = do.call(
    what = get_fixed_glmnet,
    args = list(training_data = df, training_vfold = df1, par= T)  
    )
  imp_mat <- return_feature_importance(glm_wflow_final_fit)
  return(imp_mat)
}
