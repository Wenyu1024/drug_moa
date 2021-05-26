# 1 load data, function and packages
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/get_fixed_glmnet.R')
source('/projappl/project_2003466/drug_moa/sensitivity_fitting/function/estimate_performance_uniform.R')
load("/scratch/project_2003466/prism/prism_sensitivityfitting_serverinput.RData")
library(tidyverse, quietly = T)
library(janitor,quietly = T)
print("start")

# args <- as.numeric( commandArgs( TRUE ) )
# job_id <- args[1] #841
# array job does not allow id over 999
job_id <- 947
# id_list <- scan("/projappl/project_2003466/drug_moa/predicting_prism/test_write.txt")
# job_id <- id_list[job_id]
print(job_id)
print(prism_data$name[job_id])

# 2 generate the input format of the sensitivity data for a specific drug
sen_df0 <- prism_data %>% 
  slice(job_id) %>% 
  pull(sensitivity) 
sen_df0 <- sen_df0[[1]]

# 3.1 IC50 based estimation
sen_df <- sen_df0  %>% 
  # select(depmap_id,auc) %>% 
  select(depmap_id, ic50 ) %>% 
  as_tibble() %>% 
  rename(y= ic50) %>%
  drop_na() %>% 
  group_by(depmap_id) %>% 
  summarize(y= mean(y)) %>% 
  ungroup() %>% 
  inner_join( 
    y= clean_names(ces1_478)%>% 
      rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
    by= c("depmap_id"="dep_map_id")) %>% 
  # rename(y= auc) %>% 
  select(-depmap_id)

cor_res = cor(x = sen_df$y, y = sen_df %>% select(-y), use="complete.obs")

select_predictors1 = which(abs(cor_res)> 0.15) +1
select_predictors2 = which(abs(cor_res)> 0.10) +1
select_predictors3 = which(abs(cor_res)> quantile(abs(cor_res),0.95))+1
print(quantile(abs(cor_res),0.95)) 
select_predictors4 = which(abs(cor_res)> quantile(abs(cor_res),0.9)) +1
print(quantile(abs(cor_res),0.9))
select_predictors5 = which(abs(cor_res)> quantile(abs(cor_res),0.75))+1
print(quantile(abs(cor_res),0.75))
select_predictors6 = which(abs(cor_res)> quantile(abs(cor_res),0.5))+1
print(quantile(abs(cor_res),0.50))


ces1_ic50_1= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors1)),fun_name = get_fixed_glmnet)
print(ces1_ic50_1)
ces1_ic50_2= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors2)),fun_name = get_fixed_glmnet)
print(ces1_ic50_2)
ces1_ic50_3= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors3)),fun_name = get_fixed_glmnet)
print(ces1_ic50_3)
ces1_ic50_4= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors4)),fun_name = get_fixed_glmnet)
print(ces1_ic50_4)
ces1_ic50_5= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors5)),fun_name = get_fixed_glmnet)
print(ces1_ic50_5)
ces1_ic50_6= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors6)),fun_name = get_fixed_glmnet)
print(ces1_ic50_6)

# 3.2 AUC based estimation
sen_df <- sen_df0  %>%
  # select(depmap_id,auc) %>%
  select(depmap_id, auc ) %>%
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
  # rename(y= auc) %>%
  select(-depmap_id)

# cor_res = cor(x = sen_df$y, y = sen_df %>% select(-y), use="complete.obs")

select_predictors1 = which(abs(cor_res)> 0.15)+1
select_predictors2 = which(abs(cor_res)> 0.10)+1
select_predictors3 = which(abs(cor_res)> quantile(abs(cor_res),0.95))+1
print(quantile(abs(cor_res),0.95))
select_predictors4 = which(abs(cor_res)> quantile(abs(cor_res),0.9))+1
print(quantile(abs(cor_res),0.9))
select_predictors5 = which(abs(cor_res)> quantile(abs(cor_res),0.75))+1
print(quantile(abs(cor_res),0.75))
select_predictors6 = which(abs(cor_res)> quantile(abs(cor_res),0.5))+1
print(quantile(abs(cor_res),0.50))

# 3 estimate out-of-bag prediction accuracy
ces1_auc_1= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors1)),fun_name = get_fixed_glmnet)
print(ces1_auc_1)
ces1_auc_2= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors2)),fun_name = get_fixed_glmnet)
print(ces1_auc_2)
ces1_auc_3= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors3)),fun_name = get_fixed_glmnet)
print(ces1_auc_3)
ces1_auc_4= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors4)),fun_name = get_fixed_glmnet)
print(ces1_auc_4)
ces1_auc_5= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors5)),fun_name = get_fixed_glmnet)
print(ces1_auc_5)
ces1_auc_6= estimate_performance_par(df=sen_df%>% select_at(c(1,select_predictors6)),fun_name = get_fixed_glmnet)
print(ces1_auc_6)


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
# file_name = paste0( '//project_2003466/forward_modelling/prism_ic50/drug_prism_elastic' ,job_id,".RData" )
# file_name1 = paste0( '/scratch/project_2003466/forward_modelling/prism_spearman/drug_prism_featureimp' ,job_id,".RData" )

# save(list = c('ces1_perf'), file = file_name)
file_name= paste0( '/projappl/project_2003466/drug_moa/predicting_prism/' , job_id,".RData" )

# save(list = c('ces1_perf'), file = file_name)
# save(list = c('ces1_feature'), file = file_name1)
print("all success")
save.image(file = file_name)
print("all success")



