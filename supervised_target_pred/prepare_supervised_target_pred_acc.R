load("/scratch/project_2003466/forward_modelling/targetpred_serverinput.RData")
source('/projappl/project_2003466/drug_moa/supervised_target_pred/function/get_target_pred_accuracy_batch.R')
library(furrr,quietly = T)
library(tidyverse,quietly = T)

#217
ctrp_drug_list <- res_ctrp2 %>%
  filter(!sample_size <100) %>%
  filter(broad_cpd_id %in% drug_consensus_ctrp$drug) %>%
  select(cpd_smiles,broad_cpd_id, sample_size) %>%
  drop_na() %>%
  group_by(cpd_smiles) %>%
  slice_max(n=1, order_by=sample_size) %>%
  pull(broad_cpd_id)
# 
# #57
gdsc_drug_list <- res_gdsc2 %>%
  filter(DRUG_ID %in% drug_consensus_gdsc$drug) %>%
  filter(!sample_size <100) %>%
  select(smiles,DRUG_ID, sample_size) %>%
  drop_na() %>%
  group_by(smiles) %>%
  slice_max(n=1, order_by=sample_size) %>%
  pull(DRUG_ID)
# 
# #918  
prism_drug_list <- res_prism2  %>%
  filter(BROAD_ID %in% drug_consensus_prism$drug) %>%
  filter(!sample_size <100) %>%
  select(smiles,BROAD_ID, sample_size) %>%
  drop_na() %>%
  group_by(smiles) %>%
  slice_max(n=1, order_by=sample_size) %>%
  pull(BROAD_ID)
# 
# length(intersect(ctrp_drug_list, unique(ctrp_target_binary$drug))) #196
# length(intersect(ctrp_drug_list, unique(ctrp_target_dtc$drug))) #183
# length(intersect(gdsc_drug_list, unique(gdsc_target_binary$drug))) #49
# length(intersect(gdsc_drug_list, unique(gdsc_target_dtc$drug))) #52
# length(intersect(prism_drug_list, unique(prism_target_binary$drug))) #805
# length(intersect(prism_drug_list, unique(prism_target_dtc$drug))) #137
# 
plan(multicore)
#ctrp
acc_df_ctrp_binary_fold_cv <- get_target_pred_accuracy_batch(
  dataset= "ctrp",
  estimation_method= "fold_cv" ,
  target_source= ctrp_target_binary,
  id_list = ctrp_drug_list)

acc_df_ctrp_binary_loocv <- get_target_pred_accuracy_batch(
  dataset= "ctrp",  estimation_method= "loocv" ,
  target_source= ctrp_target_binary,
  id_list = ctrp_drug_list)

acc_df_ctrp_dtcbinary_foldcv <- get_target_pred_accuracy_batch(
  dataset= "ctrp", estimation_method= "fold_cv" ,
  target_source= ctrp_target_dtc %>% filter(binding_score >0.4) %>% select(- binding_score),
  id_list = ctrp_drug_list
)

acc_df_ctrp_dtcbinary_loocv <- get_target_pred_accuracy_batch(
  dataset= "ctrp",  estimation_method= "loocv" ,
  target_source= ctrp_target_dtc %>% filter(binding_score >0.4) %>% select(- binding_score),
  id_list = ctrp_drug_list
)

acc_df_ctrp_dtc_foldcv <- get_target_pred_accuracy_batch(
  dataset= "ctrp", estimation_method= "fold_cv" , target_source= ctrp_target_dtc%>% filter(binding_score >0),
  id_list = ctrp_drug_list)

acc_df_ctrp_dtc_loocv <- get_target_pred_accuracy_batch(
  dataset= "ctrp",  estimation_method= "loocv" ,
  target_source= ctrp_target_dtc%>% filter(binding_score >0),
  id_list = ctrp_drug_list)


#gdsc
acc_df_gdsc_binary_fold_cv <- get_target_pred_accuracy_batch(
  dataset= "gdsc",  estimation_method= "fold_cv" , target_source= gdsc_target_binary,
  id_list = gdsc_drug_list)

acc_df_gdsc_binary_loocv <- get_target_pred_accuracy_batch(
  dataset= "gdsc",  estimation_method= "loocv" ,
  target_source= gdsc_target_binary,
  id_list = gdsc_drug_list
)

acc_df_gdsc_dtcbinary_foldcv <- get_target_pred_accuracy_batch(
  dataset= "gdsc", estimation_method= "fold_cv" ,
  target_source= gdsc_target_dtc %>% filter(binding_score >0.4) %>% select(- binding_score),
  id_list = gdsc_drug_list
)

acc_df_gdsc_dtcbinary_loocv <- get_target_pred_accuracy_batch(
  dataset= "gdsc",  estimation_method= "loocv" ,
  target_source= gdsc_target_dtc %>% filter(binding_score >0.4) %>% select(- binding_score),
  id_list = gdsc_drug_list
)


acc_df_gdsc_dtc_foldcv <- get_target_pred_accuracy_batch(
  dataset= "gdsc", estimation_method= "fold_cv" ,
  target_source= gdsc_target_dtc%>% filter(binding_score >0),
  id_list = gdsc_drug_list
)

acc_df_gdsc_dtc_loocv <- get_target_pred_accuracy_batch(
  dataset= "gdsc",  estimation_method= "loocv" ,
  target_source= gdsc_target_dtc%>% filter(binding_score >0),
  id_list = gdsc_drug_list
)


#prism
acc_df_prism_binary_fold_cv <- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "fold_cv" ,
  target_source= prism_target_binary,
  id_list = prism_drug_list)

acc_df_prism_binary_loocv <- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "loocv" ,
  target_source= prism_target_binary,
  id_list = prism_drug_list)

acc_df_prism_dtcbinary_foldcv <- get_target_pred_accuracy_batch(
  dataset= "prism", estimation_method= "fold_cv" ,
  target_source= prism_target_dtc %>% filter(binding_score >0.4) %>% select(- binding_score),
  id_list = prism_drug_list
)

acc_df_prism_dtcbinary_loocv <- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "loocv" ,
  target_source= prism_target_dtc %>% filter(binding_score >0.4) %>% select(- binding_score),
  id_list = prism_drug_list
)

acc_df_prism_dtc_foldcv <- get_target_pred_accuracy_batch(
  dataset= "prism", estimation_method= "fold_cv" ,
  target_source= prism_target_dtc%>% filter(binding_score >0),
  id_list = prism_drug_list)

acc_df_prism_dtc_loocv <- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "loocv" ,
  target_source= prism_target_dtc%>% filter(binding_score >0),
  id_list = prism_drug_list)

plan(sequential)
save.image("/scratch/project_2003466/forward_modelling/targetpred_output_simplefiltering.RData")

load("/scratch/project_2003466/forward_modelling/forwardmodelling_all_new.RData")
res_prism_long <- res_prism2
res_prism2 <- res_prism_long %>%
  select(BROAD_ID, .metric, input, estimate, name, smiles,drug_category,sample_size) %>% 
  filter(.metric== "spearman coef") %>% 
  pivot_wider(names_from = input, values_from= estimate) %>% 
  mutate(acc_sen= (ceres+ces+demeter2)/3) %>% 
  drop_na() %>% 
  distinct() 
plan(multicore)
prism_drug_list <- res_prism2  %>%
  filter(BROAD_ID %in% drug_consensus_prism$drug) %>%
  filter(!sample_size <100) %>%
  select(smiles,BROAD_ID, sample_size) %>%
  drop_na() %>%
  group_by(smiles) %>%
  slice_max(n=1, order_by=sample_size) %>%
  pull(BROAD_ID)
prism_drug_list_drh <- intersect(prism_drug_list, prism_target_binary$drug)
## PRISM good bad drug exploration
prism_bad_drug_list_half <- res_prism2 %>% 
  filter(BROAD_ID %in% prism_drug_list_drh) %>% 
  mutate(acc_sum= (ceres+ces+demeter2)) %>% 
  arrange(acc_sum) %>% 
  slice(1:400) %>% 
  pull(BROAD_ID)

prism_good_drug_list_half <- res_prism2 %>% 
  filter(BROAD_ID %in% prism_drug_list_drh) %>% 
  mutate(acc_sum= (ceres+ces+demeter2)) %>% 
  arrange(acc_sum) %>% 
  slice(401:805) %>% 
  pull(BROAD_ID)


prism_bad_drug_list_300 <- res_prism2 %>% 
  filter(BROAD_ID %in% prism_drug_list_drh) %>% 
  mutate(acc_sum= (ceres+ces+demeter2)) %>% 
  arrange(acc_sum) %>% 
  slice(1:300) %>% 
  pull(BROAD_ID)

prism_good_drug_list_300 <- res_prism2 %>% 
  filter(BROAD_ID %in% prism_drug_list_drh) %>% 
  mutate(acc_sum= (ceres+ces+demeter2)) %>% 
  arrange(acc_sum) %>% 
  slice(506:805) %>% 
  pull(BROAD_ID)


prism_bad_drug_list_200 <- res_prism2 %>% 
  filter(BROAD_ID %in% prism_drug_list_drh) %>% 
  mutate(acc_sum= (ceres+ces+demeter2)) %>% 
  arrange(acc_sum) %>% 
  slice(1:200) %>% 
  pull(BROAD_ID)

prism_good_drug_list_200 <- res_prism2 %>% 
  filter(BROAD_ID %in% prism_drug_list_drh) %>% 
  mutate(acc_sum= (ceres+ces+demeter2)) %>% 
  arrange(acc_sum) %>% 
  slice(606:805) %>% 
  pull(BROAD_ID)


prism_bad_drug_list_100 <- res_prism2 %>% 
  filter(BROAD_ID %in% prism_drug_list_drh) %>% 
  mutate(acc_sum= (ceres+ces+demeter2)) %>% 
  arrange(acc_sum) %>% 
  slice(1:100) %>% 
  pull(BROAD_ID)

prism_good_drug_list_100 <- res_prism2 %>% 
  filter(BROAD_ID %in% prism_drug_list_drh) %>% 
  mutate(acc_sum= (ceres+ces+demeter2)) %>% 
  arrange(acc_sum) %>% 
  slice(706:805) %>% 
  pull(BROAD_ID)


acc_df_prism_binary_fold_cv_goodhalf <- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "fold_cv" , 
  target_source= prism_target_binary,
  id_list = prism_good_drug_list_half)

acc_df_prism_binary_fold_cv_badhalf <- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "fold_cv" , 
  target_source= prism_target_binary,
  id_list = prism_bad_drug_list_half)

acc_df_prism_binary_fold_cv_good300<- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "fold_cv" , 
  target_source= prism_target_binary,
  id_list = prism_good_drug_list_300)

acc_df_prism_binary_fold_cv_bad300 <- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "fold_cv" , 
  target_source= prism_target_binary,
  id_list = prism_bad_drug_list_300)

acc_df_prism_binary_fold_cv_good200<- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "fold_cv" , 
  target_source= prism_target_binary,
  id_list = prism_good_drug_list_200)

acc_df_prism_binary_fold_cv_bad200 <- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "fold_cv" , 
  target_source= prism_target_binary,
  id_list = prism_bad_drug_list_200)

acc_df_prism_binary_fold_cv_good100<- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "fold_cv" , 
  target_source= prism_target_binary,
  id_list = prism_good_drug_list_100)

acc_df_prism_binary_fold_cv_bad100 <- get_target_pred_accuracy_batch(
  dataset= "prism",  estimation_method= "fold_cv" , 
  target_source= prism_target_binary,
  id_list = prism_bad_drug_list_100)
plan(sequential)
save.image("/scratch/project_2003466/forward_modelling/targetpred_output_prismgoodbadexplore1.RData")
