library(tidyverse)
load("~/cluster_scratch/forward_modelling/targetpred_output_simplefiltering.RData")
load("~/cluster_scratch/glmnet_modelling/target_pred_server/res_pathway_all.RData")

ctrp_data <- ctrp_data %>% 
  mutate(putative_target= map_chr(target, function(x){str_c(x$gene_symbol_of_protein_target,collapse = ", ")}))

gdsc_data <- gdsc_data %>% 
  mutate(putative_target= map_chr(target, function(x){str_c(x$gene_target,collapse = ", ")}))

prism_data <- prism_data %>% 
  mutate(putative_target= map_chr(target, function(x){str_c(x$target,collapse = ", ")}))

    

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

# 312
prism_allvalid_list <- res_prism2 %>%
  select(-variance) %>% 
  filter(.metric== "spearman coef") %>%
  pivot_wider(names_from = input, values_from= estimate) %>%
  mutate(acc_sen= (ceres+ces+demeter2)/3) %>%
  distinct() %>% 
  filter(sample_size>100) %>%
  filter(acc_sen>0.2) %>%
  pull(BROAD_ID)

pathways_df_ctrp <- tibble(drug=feature_imp_ridge_ctrp_comb1$drug, pathway=ctrp_kegg) %>% 
  unnest(c(pathway))%>% 
  arrange(padj) %>% 
  inner_join(
    y = ctrp_data %>% filter(drug_type== "single drug") %>% select_at(c(1,4,5,11)),
    by= c("drug"= "broad_cpd_id")) %>%
  filter(drug %in% ctrp_allvalid_list) %>% 
  filter(padj<0.05) 

pathways_df_gdsc <- 
  tibble(drug=feature_imp_ridge_gdsc_comb1$drug, pathway=gdsc_kegg) %>% 
  unnest(c(pathway))%>% 
  arrange(padj) %>% 
  inner_join(y = gdsc_data %>% select_at(c(1:4,10)),
             by= c("drug"= "DRUG_ID")) %>%
  filter(padj<0.05) %>% 
  filter(drug %in% gdsc_allvalid_list)

pathways_df_prism <- tibble(drug=feature_imp_ridge_prism_comb1$drug, pathway=prism_kegg) %>% 
  unnest(c(pathway))%>% 
  arrange(padj) %>% 
  inner_join(
    y = prism_data %>% select_at(c(1:4,9)),
    by= c("drug"= "BROAD_ID")) %>%
  filter(drug %in% prism_allvalid_list) %>% 
  filter(padj<0.05) 

GO_df_ctrp <- tibble(drug=feature_imp_ridge_ctrp_comb1$drug, pathway=ctrp_go) %>% 
  unnest(c(pathway))%>% 
  arrange(padj) %>% 
  inner_join(
    y = ctrp_data %>% filter(drug_type== "single drug") %>% select_at(c(1,4,5,11)),
    by= c("drug"= "broad_cpd_id")) %>%
  filter(drug %in% ctrp_allvalid_list) %>% 
  filter(padj<0.05) 

GO_df_gdsc <- 
  tibble(drug=feature_imp_ridge_gdsc_comb1$drug, pathway=gdsc_go) %>% 
  unnest(c(pathway))%>% 
  arrange(padj) %>% 
  inner_join(y = gdsc_data %>% select_at(c(1:4,10)),
             by= c("drug"= "DRUG_ID")) %>%
  filter(padj<0.05) %>% 
  filter(drug %in% gdsc_allvalid_list)

GO_df_prism <- tibble(drug=feature_imp_ridge_prism_comb1$drug, pathway=prism_go) %>% 
  unnest(c(pathway))%>% 
  arrange(padj) %>% 
  inner_join(
    y = prism_data %>% select_at(c(1:4,9)),
    by= c("drug"= "BROAD_ID")) %>%
  filter(drug %in% prism_allvalid_list) %>% 
  filter(padj<0.05) 



pathways_fgl1_ctrp <- pathways_df_ctrp %>% 
  mutate(FGL1_ass= map_int(.x = leadingEdge, .f = ~sum(str_detect(string = .x,pattern = "FGL1"))))

pathways_fgl1_gdsc <- pathways_df_gdsc %>% 
  mutate(FGL1_ass= map_int(.x = leadingEdge, .f = ~sum(str_detect(string = .x,pattern = "FGL1"))))

pathways_fgl1_prism <- pathways_df_prism %>% 
  mutate(FGL1_ass= map_int(.x = leadingEdge, .f = ~sum(str_detect(string = .x,pattern = "FGL1"))))

GO_fgl1_ctrp <- GO_df_ctrp %>% 
  mutate(FGL1_ass= map_int(.x = leadingEdge, .f = ~sum(str_detect(string = .x,pattern = "FGL1"))))

GO_fgl1_gdsc <- GO_df_gdsc %>% 
  mutate(FGL1_ass= map_int(.x = leadingEdge, .f = ~sum(str_detect(string = .x,pattern = "FGL1"))))

GO_fgl1_prism <- GO_df_prism %>% 
  mutate(FGL1_ass= map_int(.x = leadingEdge, .f = ~sum(str_detect(string = .x,pattern = "FGL1"))))

setwd("~/cluster_wrk/drug_moa/unsupervised_target_exploration/fig_res")
GO_fgl1_ctrp <- GO_fgl1_ctrp %>% filter(FGL1_ass==1) %>% 
  mutate(leadingEdge= map_chr(leadingEdge, ~str_c(collapse = ", ")))
GO_fgl1_gdsc <- GO_fgl1_gdsc %>% filter(FGL1_ass==1) %>% 
  mutate(leadingEdge= map_chr(leadingEdge, ~str_c(collapse = ", ")))
write_csv(GO_fgl1_ctrp,"GO_fgl1_ctrp.csv" )
write_csv(GO_fgl1_gdsc,"GO_fgl1_gdsc.csv" )
