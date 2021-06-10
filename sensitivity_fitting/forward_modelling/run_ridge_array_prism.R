# 1 load data, function and packages
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/get_fixed_ridge.R')
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/estimate_performance_uniform.R')
load("/scratch/project_2003466/forward_modelling/prism_input.RData")
library(tidyverse, quietly = T)
# library(furrr, quietly = T)
print("start")

args <- as.numeric( commandArgs( TRUE ) )
job_id <- args[1] 
# print(job_id)
# array job does not allow id over 999
# job_id <- 1
# some of the jobs failed due to lack of memory, these failure is 
# detected and recored in the test_write.txt,
# I then use it to rerun analysis with larger memory setting
# 
# id_list <- scan("/projappl/project_2003466/drug_moa/sensitivity_fitting/forward_modelling/test_write.txt")
# job_id <- id_list[job_id]
print(job_id)

# 2 generate the input format of the sensitivity data for a specific drug
sen_drug <- data %>% 
  slice(job_id) %>% 
  pull(sensitivity) 

sen_drug <- sen_drug[[1]]

# 3 accuracy estimation
get_df_prism <- function(sen_df= sen_drug, predictor_df){
  df <- sen_df  %>% 
    as_tibble() %>% 
    inner_join( 
      y= janitor::clean_names(predictor_df) %>% 
        rename_at(vars(!contains("dep_map_id")), 
                  .fun= ~paste0(., "_predictors")), 
      by= c("depmap_id"="dep_map_id")) %>% 
    select(-depmap_id,-ec50) %>% 
    rename(y= auc) 
  return(df)
}

ces1_df <- get_df_prism(predictor_df = ces1)
ceres_df <- get_df_prism(predictor_df = ceres)
demeter2_df <- get_df_prism(predictor_df = demeter2)
exp_df <- get_df_prism(predictor_df = exp_seq)
rm(ces1, ceres,demeter2, exp_seq,data, sen_drug)

plan(multicore)
ces1_perf= estimate_performance_par(df= ces1_df,fun_name = get_fixed_glmnet)
print("ces1 success")
rm(ces1_df)
ceres_perf= estimate_performance_par(df= ceres_df ,fun_name = get_fixed_glmnet)
print("ceres success")
rm(ceres_df)
demeter2_perf= estimate_performance_par(df= demeter2_df,fun_name = get_fixed_glmnet)
print("demeter2 success")
rm(demeter2_df)
exp_perf= estimate_performance_par(df= exp_df,fun_name = get_fixed_glmnet)
print("exp success")
rm(exp_df)
plan(sequential)

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
file_name = paste0( '/scratch/project_2003466/forward_modelling/prism_spearman/drug_' ,job_id,".RData" )
# save(list = c('ces1_perf'), file = file_name)

save(list = c('ces1_perf',"ceres_perf", "demeter2_perf","exp_perf"), file = file_name)

print("all success")




