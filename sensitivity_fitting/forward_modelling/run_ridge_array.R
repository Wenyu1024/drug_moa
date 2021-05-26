source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/get_fixed_glmnet.R')
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/estimate_performance_uniform.R')
load("/scratch/project_2003466/glmnet_modelling_cluster/input.RData")
library(tidyverse, quietly = T)
library(furrr, quietly = T)
print("start")

args <- as.numeric( commandArgs( TRUE ) )

job_id <- args[1]
print(job_id)

sen_df <- data %>% 
  slice(job_id) %>% 
  pull(sensitivity) 
sen_df <- sen_df[[1]]

sen_df <- sen_df  %>% 
  as_tibble() %>% 
  inner_join( 
    y= clean_names(predictor_df)%>% 
      rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
    by= c("DepMap_ID"="dep_map_id")) %>% 
  select(-master_ccl_id, -apparent_ec50_umol, -DepMap_ID) %>% 
  rename(y= area_under_curve)

ces1_perf= estimate_performance_par(sen_df = sen_df, predictor_df = ces1,fun_name = get_fixed_glmnet)
print("ces1 success")
ces2_perf= estimate_performance_par(sen_df = sen_df, predictor_df = ces2,fun_name = get_fixed_glmnet)
print("ces2 success")
ceres_perf= estimate_performance_par(sen_df = sen_df, predictor_df = ceres,fun_name = get_fixed_glmnet)
print("ceres success")
demeter2_perf= estimate_performance_par(sen_df = sen_df, predictor_df = demeter2,fun_name = get_fixed_glmnet)
print("demeter2 success")
exp_perf= estimate_performance_par(sen_df = sen_df, predictor_df = exp_seq_pca,fun_name = get_fixed_glmnet)
print("exp success")

file_name = paste0( '/scratch/project_2003466/forward_modelling/ctrp_spearman/drug_' ,job_id,".RData" )
# save(list = c('ces1_perf'), file = file_name)

save(list = c('ces1_perf', 'ces2_perf',"ceres_perf", "demeter2_perf","exp_perf"), file = file_name)

print("all success")
