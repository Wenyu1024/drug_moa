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

get_feature_imp <- function(sen_df, predictor_df){
  library(janitor, lib.loc = "/projappl/project_2003466/project_rpackages", quietly = T )
  df <- sen_df  %>%   as_tibble() %>% 
    inner_join( 
      y= clean_names(predictor_df)%>% 
        rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
      by= c("DepMap_ID"="dep_map_id")) %>% 
    select(-master_ccl_id, -apparent_ec50_umol, -DepMap_ID) 
  
  df1 <-  df %>% vfold_cv(v = 5,repeats = 3)
  
  glm_wflow_final_fit = get_fixed_glmnet(
    training_data = df,
    training_vfold = df1)
  
  imp_mat <- return_feature_importance(glm_wflow_final_fit)
  return(imp_mat)
}
