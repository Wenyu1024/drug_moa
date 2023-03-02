# get signatures for Anna

library(tidyverse)
load("~/cluster_scratch/forward_modelling/targetpred_output_simplefiltering.RData")


ctrp_allvalid_list <- res_ctrp2 %>% 
  select(-.metric, -estimate, -variance, -input) %>% 
  distinct() %>% 
  filter(!sample_size <100) %>% 
  filter(!drug_pair) %>% 
  select(cpd_smiles,broad_cpd_id, sample_size) %>%
  drop_na() %>% 
  group_by(cpd_smiles) %>% 
  slice_max(n=1, order_by=sample_size) %>% 
  pull(broad_cpd_id)

#n= 175
gdsc_allvalid_list <- res_gdsc2 %>% 
  select(-.metric, -estimate, -variance, -input) %>% 
  distinct() %>% 
  filter(!sample_size <100) %>% 
  # select(smiles,DRUG_ID, sample_size) %>%
  # drop_na() %>% 
  group_by(DRUG_NAME) %>%
  slice_max(n=1, order_by=sample_size) %>% 
  pull(DRUG_ID)

res_prism2_long <- res_prism2
res_prism2 <- res_prism2_long %>%
  select(-variance) %>% 
  filter(.metric== "spearman coef") %>%
  pivot_wider(names_from = input, values_from= estimate) %>%
  mutate(acc_sen= (ceres+ces+demeter2)/3) %>%
  distinct() 

prism_good_fitted_drugs <- res_prism2 %>%
  filter(sample_size>100) %>%
  filter(acc_sen>0.2) %>%
  arrange(BROAD_ID)


length(intersect(ctrp_allvalid_list, unique(ctrp_target_binary$drug))) #360
length(intersect(ctrp_allvalid_list, unique(ctrp_target_dtc$drug))) #209

length(intersect(gdsc_allvalid_list, unique(gdsc_target_binary$drug))) #123
length(intersect(gdsc_allvalid_list, unique(gdsc_target_dtc$drug))) #104




ctrp_esssig <- feature_imp_ridge_ctrp_comb1 %>% filter(drug %in% ctrp_allvalid_list)
prism_esssig <- feature_imp_ridge_prism_comb1 %>% filter(drug %in% prism_good_fitted_drugs$BROAD_ID)
setwd("~/cluster_scratch/forward_modelling/feature_imp/")
dir.create("anna_collboration")
setwd("./anna_collboration")
write_csv(ctrp_esssig, "ctrp_ess_sig.csv")
write_csv(prism_esssig, "prism_ess_sig.csv")
write_csv(drug_consensus, "drug_consensus.csv")