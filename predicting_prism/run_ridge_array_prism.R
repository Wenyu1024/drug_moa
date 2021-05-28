# 1 load data, function and packages
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/get_fixed_elasticnet.R')
# source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/get_fixed_ridge.R')
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/estimate_performance_uniform.R')
load("/scratch/project_2003466/prism/prism_sensitivityfitting_serverinput.RData")
library(tidyverse, quietly = T)
library(janitor,quietly = T)
print("start")

args <- as.numeric( commandArgs( TRUE ) )
job_id <- args[1] #841
# array job does not allow id over 999
# job_id <- 1
# id_list <- scan("/projappl/project_2003466/drug_moa/predicting_prism/test_write.txt")
# job_id <- id_list[job_id]
print(job_id)
print(prism_data$name[job_id])

# 2 generate the input format of the sensitivity data for a specific drug
sen_df0 <- prism_data %>% 
  slice(job_id) %>% 
  pull(sensitivity) 
sen_df0 <- sen_df0[[1]]

# 3 accuracy estimation
sen_df <- sen_df0  %>% 
  select(depmap_id,auc) %>%
  # select(depmap_id, ic50 ) %>% 
  as_tibble() %>% 
  rename(y= auc) %>%
  drop_na() %>% 
  group_by(depmap_id) %>% 
  summarize(y= mean(y)) %>% 
  ungroup() %>% 
  inner_join( 
    y= clean_names(ces1_478)%>% 
      rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
    by= c("depmap_id"="dep_map_id")) %>% 
  select(-depmap_id)


ces1_perf= estimate_performance_par(df=sen_df,fun_name = get_fixed_glmnet)



# 4 Fix model and derive feature importance
# return_feature_importance <- function(fit){
#   fit %>% 
#     tidy() %>% 
#     filter(term!="(Intercept)") %>% 
#     arrange(term) %>% 
#     pull(estimate)
# }
# 
# get_feature_imp <- function(df){
#   library(janitor, lib.loc = "/projappl/project_2003466/project_rpackages", quietly = T )
#   df1 <-  df %>% vfold_cv(v = 5,repeats = 3)
#   glm_wflow_final_fit = get_fixed_glmnet(
#     training_data = df,
#     training_vfold = df1)
#   
#   imp_mat <- return_feature_importance(glm_wflow_final_fit)
#   return(imp_mat)
# }
# 
# 
# ces1_feature <- get_feature_imp(df = sen_df)

# 5 write out results
file_name = paste0( '/scratch/project_2003466/forward_modelling/test_replicatibility_auc_spearman_filtering/drug_' ,job_id,".RData" )
# file_name1 = paste0( '/scratch/project_2003466/forward_modelling/prism_spearman/drug_prism_featureimp' ,job_id,".RData" )

# save(list = c('ces1_perf'), file = file_name)
# file_name= paste0( '/projappl/project_2003466/drug_moa/predicting_prism/' , job_id,".RData" )
print("all success")
save(list = c('ces1_perf'), file = file_name)
# save(list = c('ces1_feature'), file = file_name1)
print("all success")
# save.image(file = file_name)




