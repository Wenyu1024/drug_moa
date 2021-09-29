# 1 load data, function and packages
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/forward_modelling/function/get_fixed_ridge.R')
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/forward_modelling/function/estimate_performance_par.R')
load("/scratch/project_2003466/forward_modelling/gdsc_input.RData")
library(tidyverse, quietly = T)
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

get_df_gdsc <- function(sen_df= sen_drug, predictor_df){
  df <- sen_df  %>% 
    as_tibble() %>% 
    inner_join( 
      y= janitor::clean_names(predictor_df)%>% 
        rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
      by= c("DepMap_ID"="dep_map_id")) %>% 
    select( -LN_IC50, -DepMap_ID) %>% 
    rename(y= AUC)
  return(df)
}

ces1_df <- get_df_gdsc(predictor_df = ces1)
ceres_df <- get_df_gdsc(predictor_df = ceres)
demeter2_df <- get_df_gdsc(predictor_df = demeter2)
exp_df <- get_df_gdsc(predictor_df = exp_seq)
rm(ces1, ceres,demeter2, exp_seq,data, sen_drug)

# 3 feature imp derivation
plan(multicore)
ces1_perf= get_feature_imp(df= ces1_df,fun_name = get_fixed_glmnet)
print("ces1 success")
rm(ces1_df)
ceres_perf= get_feature_imp(df= ceres_df ,fun_name = get_fixed_glmnet)
print("ceres success")
rm(ceres_df)
demeter2_perf= get_feature_imp(df= demeter2_df,fun_name = get_fixed_glmnet)
print("demeter2 success")
rm(demeter2_df)
exp_perf= get_feature_imp(df= exp_df,fun_name = get_fixed_glmnet)
print("exp success")
rm(exp_df)
plan(sequential)

# 4 write out results
file_name = paste0( '/scratch/project_2003466/forward_modelling/gdsc_spearman_feature_imp/drug_' ,job_id,".RData" )
save(list = c('ces1_perf',"ceres_perf", "demeter2_perf","exp_perf"), file = file_name)

print("all success")

