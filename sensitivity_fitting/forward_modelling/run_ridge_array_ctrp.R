# 1 load data, function and packages
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/forward_modelling/function/get_fixed_ridge.R')
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/forward_modelling/function/estimate_performance_par.R')
load("/scratch/project_2003466/forward_modelling/ctrp_input.RData")
library(tidyverse, quietly = T)
library(furrr, quietly = T)
library(tictoc)
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

# 3 accuracy estimation
tictoc::tic()
plan(multisession)

ces1_perf= estimate_performance_par(df= ces1_df,fun_name = get_fixed_glmnet)
print("ces1 success")
# file_name = paste0( '/scratch/project_2003466/forward_modelling/ctrp_modelres/drug_' ,job_id,".RData" )
# save('ces1_perf', file = file_name)
# plan(sequential)
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
tictoc::toc()

# # 4 write out results
file_name = paste0( '/scratch/project_2003466/forward_modelling/ctrp_modelresnew/drug_' ,job_id,".RData" )
save(list = c('ces1_perf',"ceres_perf", "demeter2_perf","exp_perf"), file = file_name)

print("all success")
