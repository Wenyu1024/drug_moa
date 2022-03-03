# 1 load data, function and packages
library(tidyverse, quietly = T)
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/forward_modelling/function/get_fixed_ridge.R')
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/forward_modelling/function/derive_feature_imp.R')
load("/scratch/project_2003466/forward_modelling/ctrp_input.RData")
library(furrr, quietly = T)
print("start")

args <- as.numeric( commandArgs( TRUE ) )
job_id <- args[1]
print(job_id)

# 2 generate the input format of the sensitivity data for a specific drug
sen_drug <- data %>% 
  slice(job_id) %>% 
  pull(sensitivity) 
sen_drug <- sen_drug[[1]]

get_df_ctrp <- function(sen_df= sen_drug, predictor_df){
  df <- sen_df  %>% 
    as_tibble() %>% 
    inner_join( 
      y= janitor::clean_names(predictor_df)%>% 
        rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
      by= c("DepMap_ID"="dep_map_id")) %>% 
    select( -apparent_ec50_umol, -DepMap_ID) %>% 
    rename(y= area_under_curve)
  return(df)
}

ces1_df <- get_df_ctrp(predictor_df = ces1)
ceres_df <- get_df_ctrp(predictor_df = ceres)
demeter2_df <- get_df_ctrp(predictor_df = demeter2)
exp_df <- get_df_ctrp(predictor_df = exp_seq)
rm(ces1, ceres,demeter2, exp_seq,data, sen_drug)
# save.image("/projappl/project_2003466/drug_moa/sensitivity_fitting/forward_modelling/function/tmp.RData")
# 3 feature imp derivation
plan(multicore)
ces1_imp= get_feature_imp(df= ces1_df,fun_name = get_fixed_glmnet)
print("ces1 success")
rm(ces1_df)
ceres_imp= get_feature_imp(df= ceres_df,fun_name = get_fixed_glmnet)
print("ceres success")
rm(ceres_df)
demeter2_imp= get_feature_imp(df= demeter2_df,fun_name = get_fixed_glmnet)
print("demeter2 success")
rm(demeter2_df)
exp_imp= get_feature_imp(df= exp_df,fun_name = get_fixed_glmnet)
print("exp success")
rm(exp_df)
plan(sequential)

# 4 write out results
save(list = c('ces1_imp',"ceres_imp", "demeter2_imp","exp_imp"), file = file_name)
file_name = paste0( '/scratch/project_2003466/forward_modelling/ctrp_spearman_feature_imp/drug_' ,job_id,".RData" )

# 5 model summary derivation
# plan(multicore)
# ces1_sum= get_model_summary(df= ces1_df,fun_name = get_fixed_glmnet)
# print("ces1 success")
# rm(ces1_df)
# ceres_sum= get_model_summary(df= ceres_df,fun_name = get_fixed_glmnet)
# print("ceres success")
# rm(ceres_df)
# demeter2_sum= get_model_summary(df= demeter2_df,fun_name = get_fixed_glmnet)
# print("demeter2 success")
# rm(demeter2_df)
# exp_sum= get_model_summary(df= exp_df,fun_name = get_fixed_glmnet)
# print("exp success")
# rm(exp_df)
# plan(sequential)
# 
# print("all success")
# file_name = paste0( '/scratch/project_2003466/forward_modelling/ctrp_model_summary/drug_' ,job_id,".RData" )
# save(list = c('ces1_sum',"ceres_sum", "demeter2_sum","exp_sum"), file = file_name)
