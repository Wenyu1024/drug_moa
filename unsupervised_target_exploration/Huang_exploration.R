library(tidyverse)
source("~/cluster_wrk/drug_moa/unsupervised_target_exploration/get_drug_frac.R")
# 
# CTRP_binary_PPIold_data <- read_csv( "~/cluster_scratch/prior/CTRP_binary_PPIold_data.csv")
# GDSC_binary_PPIold_data <- read_csv( "~/cluster_scratch/prior/GDSC_binary_PPIold_data.csv")
# PRISM_binary_PPIold_data <- read_csv( "~/cluster_scratch/prior/PRISM_binary_PPIold_data.csv")
# 
# 
# network_idf <- bind_rows(
#   CTRP_binary_PPIold_data,
#   GDSC_binary_PPIold_data %>% mutate(drug= as.character(drug)),
#   PRISM_binary_PPIold_data
# ) %>%
#   filter(imptype== "ConSen-Sig") %>%
#   group_by(dataset, drug) %>%
#   mutate(rank=rank(desc(imp))) %>%
#   filter(gene %in% c("FGL1","FOXO4"))


# Use the CES score, use the CTRP dataset. Focus on FGL1, get the loocv and unsupervised target prediction for the drug.
# do not change CES score but do not limit to comparable drugs.



# check the raw data rather than long tibble since there is some filter used for the long tibble?

load("~/cluster_scratch/forward_modelling/targetpred_output_simplefiltering.RData")


feature_imp_ridge_ctrp_ces1 <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_ces1.csv")
feature_imp_ridge_gdsc_ces1 <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_ces1.csv")
feature_imp_ridge_prism_ces1 <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_ces1.csv")
# 


FGL1associateddrugs <- feature_imp_ridge_ctrp_ces1 %>%
  pivot_longer(cols= -drug, names_to = "gene", values_to="imp") %>%
  group_by(drug) %>%
  mutate(gene_rank= rank(desc(imp))) %>%
  ungroup() %>% 
  filter(gene %in% c( "NAT10", "FGL1", "FOXO4")) %>% 
  filter(!(gene_rank>106)) %>%
  pull(drug)

FGL1associateddrug_topgenedf <- feature_imp_ridge_ctrp_ces1 %>%
  filter(drug %in% FGL1associateddrugs ) %>% 
  pivot_longer(cols= -drug, names_to = "gene", values_to="imp") %>%
  group_by(drug) %>%
  mutate(gene_rank= rank(desc(imp))) %>%
  filter(gene %in% c( "NAT10", "FGL1", "FOXO4")) %>% 
  filter(!(gene_rank>106)) 



#fitting accuracy and unsupervised prediction accuracy. drug annotation. AUC_TE1, AUC_ECA

load("~/cluster_scratch/unsupervised_target_pred/unsupervised_target_pred_ConSenSig.RData")
tmp <- res_ctrp2 %>% 
  filter(broad_cpd_id %in% FGL1associateddrugs) %>% 
  select(cpd_name, broad_cpd_id,ces, exp, drug_type, sample_size,ces_variance) %>% 
  left_join(FGL1associateddrug_topgenedf %>% 
              # filter(gene=="FGL1") %>% 
              select(-imp), 
            by=c("broad_cpd_id"="drug")) %>% 
  inner_join(ctrp_data %>% select(broad_cpd_id, target, sensitivity),by="broad_cpd_id") %>% 
  mutate(TE1_auc =map_dbl(.x = sensitivity, 
                          .f=function(df){
                            auc <- df %>% filter(DepMap_ID== "ACH-000647") %>% pull(area_under_curve)
                            if(length(auc)==0) {auc <- NA}
                          return(auc)}
                          )) %>% 
  mutate(TE1_EC50 =map_dbl(.x = sensitivity, 
                           .f=function(df){
                             EC_50 <- df %>% filter(DepMap_ID== "ACH-000647") %>% pull(apparent_ec50_umol)
                             if(length(EC_50)==0) {EC_50 <- NA}
                             return(EC_50)
                             })) %>% 
  mutate(target= map_chr(target, .f = function(x){str_c(x$gene_symbol_of_protein_target, collapse = ";")})) %>% 
  select(-sensitivity) %>% 
  inner_join(ctrp_unsupervised_auc,by=c("broad_cpd_id"="drug"))
  
  
  #ACH-000647

write_csv(tmp, "FGL1.csv")



FGL1associateddrugs <- feature_imp_ridge_gdsc_ces1 %>%
  pivot_longer(cols= -drug, names_to = "gene", values_to="imp") %>%
  group_by(drug) %>%
  mutate(gene_rank= rank(desc(imp))) %>%
  ungroup() %>% 
  filter(gene %in% c( "NAT10", "FGL1", "FOXO4")) %>% 
  filter(!(gene_rank>106)) %>%
  pull(drug)

FGL1associateddrug_topgenedf <- feature_imp_ridge_gdsc_ces1 %>%
  filter(drug %in% FGL1associateddrugs ) %>% 
  pivot_longer(cols= -drug, names_to = "gene", values_to="imp") %>%
  group_by(drug) %>%
  mutate(gene_rank= rank(desc(imp))) %>%
  filter(gene %in% c( "NAT10", "FGL1", "FOXO4")) %>% 
  filter(!(gene_rank>106)) 


tmp <- res_gdsc2 %>% 
  filter(DRUG_ID %in% FGL1associateddrugs) %>% 
  select(DRUG_NAME, DRUG_ID,ces, exp,  sample_size,ces_variance) %>% 
  left_join(FGL1associateddrug_topgenedf %>% 
              # filter(gene=="FGL1") %>% 
              select(-imp), 
            by=c("DRUG_ID"="drug")) %>% 
  inner_join(gdsc_data %>% select(DRUG_ID, target, sensitivity, PATHWAY_NAME, DRUG_TYPE ),by="DRUG_ID") %>% 
  mutate(TE1_auc =map_dbl(.x = sensitivity, 
                          .f=function(df){
                            auc <- df %>% filter(DepMap_ID== "ACH-000647") %>% pull(AUC)
                            if(length(auc)==0) {auc <- NA}
                            return(auc)}
  )) %>% 
  mutate(TE1_EC50 =map_dbl(.x = sensitivity, 
                           .f=function(df){
                             LN_IC50 <- df %>% filter(DepMap_ID== "ACH-000647") %>% pull(LN_IC50)
                             if(length(LN_IC50)==0) {LN_IC50 <- NA}
                             return(LN_IC50)
                           })) %>% 
  mutate(target= map_chr(target, .f = function(x){str_c(x$gene_target, collapse = ";")})) %>% 
  select(-sensitivity) %>% 
  inner_join(gdsc_unsupervised_auc,by=c("DRUG_ID"="drug"))


write_csv(tmp, "FGL1_GDSC.csv")

# check the situation in PRISM dataset

# FGL1associateddrugs <- feature_imp_ridge_prism_ces1 %>%
#   pivot_longer(cols= -drug, names_to = "gene", values_to="imp") %>%
#   group_by(drug) %>%
#   mutate(gene_rank= rank(desc(imp))) %>%
#   ungroup() %>% 
#   filter(gene %in% c("FGL1","FOXO4")) %>% 
#   filter(!(gene_rank>10)) %>%
#   pull(drug)
# 
# FGL1associateddrug_topgenedf <- feature_imp_ridge_prism_ces1 %>%
#   filter(drug %in% FGL1associateddrugs ) %>% 
#   pivot_longer(cols= -drug, names_to = "gene", values_to="imp") %>%
#   group_by(drug) %>%
#   mutate(gene_rank= rank(desc(imp))) %>%
#   filter(!(gene_rank>10)) 
# 
# 
# 
# tmp <- res_prism2 %>% 
#   filter(BROAD_ID %in% FGL1associateddrugs) %>% 
#   select(name, BROAD_ID,ces, exp, drug_category, sample_size,ces_variance) %>% 
#   left_join(FGL1associateddrug_topgenedf %>% filter(gene %in% c("FGL1", "FOXO4")) %>% select( -imp), 
#             by=c("BROAD_ID"="drug")) %>% 
#   inner_join(prism_data %>% select(BROAD_ID, target, sensitivity),by="BROAD_ID") 
  # mutate(TE1_auc =map_dbl(.x = sensitivity, 
  #                         .f=function(df){df %>% filter(DepMap_ID== "ACH-000647") %>% pull(AUC)})) %>% 
  # mutate(TE1_EC50 =map_dbl(.x = sensitivity, 
  #                          .f=function(df){df %>% filter(DepMap_ID== "ACH-000647") %>% pull(apparent_ec50_umol)})) %>% 
  # mutate(target= map_chr(target, .f = function(x){str_c(x$gene_symbol_of_protein_target, collapse = ";")})) %>% 
  # select(-sensitivity) %>% 
  # inner_join(prism_unsupervised_auc,by=c(""="drug"))
