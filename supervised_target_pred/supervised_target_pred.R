load("/scratch/project_2003466/glmnet_modelling/target_pred_server/targetpred_serverinput.RData")
source('/projappl/project_2003466/analysis/downstream_feature_analysis.R')
library(future,quietly = T)
library(tidyverse,quietly = T)
library(furrr,quietly = T)
plan(multicore)

acc_ces1 <- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble2,
  predictors_tibble = feature_imp_ridge_ctrp_ces1 %>% rename(drug=broad_cpd_id) ) 

acc_ces2<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble2,
  predictors_tibble = feature_imp_ridge_ctrp_ces2 %>% rename(drug=broad_cpd_id) )

acc_ceres<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble2,
  predictors_tibble = feature_imp_ridge_ctrp_ceres %>% rename(drug=broad_cpd_id) )

acc_demeter2<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble2,
  predictors_tibble = feature_imp_ridge_ctrp_demeter2 %>% rename(drug=broad_cpd_id) )

acc_sen<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble2,
  predictors_tibble = feature_sen_ctrp %>% rename(drug=broad_cpd_id) )

acc_sen_ces1<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble2,
  predictors_tibble = feature_sen_ces1_ctrp %>% rename(drug=broad_cpd_id) )

acc_seqpca<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble2,
  predictors_tibble = feature_imp_ridge_seq_pca %>% rename(drug=broad_cpd_id) )

acc_sen_seqpca <- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble2,
  predictors_tibble = feature_sen_seqpca_ctrp %>% rename(drug=broad_cpd_id) )

plan(sequential)
print("success")


save(
  list = c("acc_ces1","acc_ces2","acc_ceres","acc_demeter2","acc_sen","acc_sen_ces1","acc_seqpca","acc_sen_seqpca"),
  file = "/scratch/project_2003466/glmnet_modelling/dtc_supervised_prediction_ridge.RData"
)



#########################################################
plan(multicore)
acc_ces1 <- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_ctrp_ces1 %>% rename(drug=broad_cpd_id) )

acc_ces2<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_ctrp_ces2 %>% rename(drug=broad_cpd_id) )

acc_ceres<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_ctrp_ceres %>% rename(drug=broad_cpd_id) )

acc_demeter2<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_ctrp_demeter2 %>% rename(drug=broad_cpd_id) )

acc_sen<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble1,
  predictors_tibble = feature_sen_ctrp %>% rename(drug=broad_cpd_id) )

acc_sen_ces1<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble1,
  predictors_tibble = feature_sen_ces1_ctrp %>% rename(drug=broad_cpd_id) )

acc_seqpca<- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_seq_pca %>% rename(drug=broad_cpd_id) )

acc_sen_seqpca <- return_target_supervised_similiarity_prediction_acc1(
  target_tibble = target_tibble1,
  predictors_tibble = feature_sen_seqpca_ctrp %>% rename(drug=broad_cpd_id) )


plan(sequential)
print("success")


save(
  list = c("acc_ces1","acc_ces2","acc_ceres","acc_demeter2","acc_sen","acc_sen_ces1","acc_seqpca","acc_sen_seqpca"),
  file = "/scratch/project_2003466/glmnet_modelling/drugbank_supervised_prediction_ridge.RData"
)

##########################################################
# How about the loo_cv without tunning


plan(multicore)
acc_ces1 <- return_target_supervised_similiarity_prediction_acc(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_ctrp_ces1 %>% rename(drug=broad_cpd_id) )

acc_ces2<- return_target_supervised_similiarity_prediction_acc(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_ctrp_ces2 %>% rename(drug=broad_cpd_id) )

acc_ceres<- return_target_supervised_similiarity_prediction_acc(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_ctrp_ceres %>% rename(drug=broad_cpd_id) )

acc_demeter2<- return_target_supervised_similiarity_prediction_acc(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_ctrp_demeter2 %>% rename(drug=broad_cpd_id) )

acc_sen<- return_target_supervised_similiarity_prediction_acc(
  target_tibble = target_tibble1,
  predictors_tibble = feature_sen_ctrp %>% rename(drug=broad_cpd_id) )

acc_sen_ces1<- return_target_supervised_similiarity_prediction_acc(
  target_tibble = target_tibble1,
  predictors_tibble = feature_sen_ces1_ctrp %>% rename(drug=broad_cpd_id) )

acc_seqpca<- return_target_supervised_similiarity_prediction_acc(
  target_tibble = target_tibble1,
  predictors_tibble = feature_imp_ridge_seq_pca %>% rename(drug=broad_cpd_id) )

acc_sen_seqpca <- return_target_supervised_similiarity_prediction_acc(
  target_tibble = target_tibble1,
  predictors_tibble = feature_sen_seqpca_ctrp %>% rename(drug=broad_cpd_id) )


plan(sequential)
print("success")


save(
  list = c("acc_ces1","acc_ces2","acc_ceres","acc_demeter2","acc_sen","acc_sen_ces1","acc_seqpca","acc_sen_seqpca"),
  file = "/scratch/project_2003466/glmnet_modelling/dtc_supervised_prediction_ridge_notunning.RData"
)


########################################################################


