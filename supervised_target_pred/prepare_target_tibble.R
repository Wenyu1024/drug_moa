library(tidyverse)
load("~/cluster_scratch/forward_modelling/forwardmodelling_all.RData")
ctrp_target_binary <- ctrp_data %>% 
  select(broad_cpd_id, target) %>% 
  unnest(target) %>% 
  drop_na() %>% 
  distinct() %>% 
  rename(drug=broad_cpd_id, target= gene_symbol_of_protein_target) %>% 
  filter(target!="")

gdsc_target_binary <- gdsc_data %>% 
  select(DRUG_ID, target) %>% 
  unnest(target) %>% 
  drop_na() %>% 
  rename(drug=DRUG_ID, target= gene_target) %>% 
  distinct()

prism_target_binary <- prism_data %>% 
  select(BROAD_ID, target) %>% 
  unnest(target) %>% 
  drop_na() %>% 
  distinct() %>% 
  rename(drug=BROAD_ID)

write_csv(ctrp_target_binary,
          path = "~/cluster_scratch/prior/ctrpv2_target_binary.csv")
write_csv(gdsc_target_binary,
          path = "~/cluster_scratch/prior/gdsc_target_binary.csv")
write_csv(prism_target_binary,
          path = "~/cluster_scratch/prior/prism_target_binary.csv")


ctrp_list_to_zia <- read_csv("~/cluster_scratch/prior/drug_id_list/ctrp_idlist.csv")
gdsc_list_to_zia <- read_csv("~/cluster_scratch/prior/drug_id_list/gdsc2_idlist.csv")
prism_list_to_zia <- read_csv("~/cluster_scratch/prior/drug_id_list/prism_idlist.csv")

ctrp_target_dtc <- read_csv("~/cluster_scratch/prior/dtc_return/ctrp_data.csv") %>% 
  drop_na() %>% distinct() %>% 
  inner_join(ctrp_list_to_zia %>% select(broad_cpd_id, InChIKey), by= c("standard_inchi_key"= "InChIKey") ) %>% 
  distinct() %>% 
  rename(binding_score= interaction_strength, drug= broad_cpd_id, target= gene_name ) %>%
  select(drug, target, binding_score) %>% 
  group_by(drug, target) %>% 
  summarise(binding_score= max(binding_score)) %>% 
  ungroup()

gdsc_target_dtc <- read_csv("~/cluster_scratch/prior/dtc_return/gdsc_data.csv") %>% 
  drop_na() %>% distinct() %>% 
  mutate(compound_id = as.character(compound_id) ) %>% 
  inner_join(gdsc_list_to_zia %>% select(DRUG_ID, PubCHEM), by= c("compound_id"= "PubCHEM") ) %>% 
  rename(binding_score= interaction_strength, drug= DRUG_ID, target= gene_name ) %>%
  distinct() %>% 
  select(drug, target, binding_score) %>% 
  group_by(drug, target) %>% 
  summarise(binding_score= max(binding_score)) %>% 
  ungroup()

prism_target_dtc <- read_csv("~/cluster_scratch/prior/dtc_return/prism_data.csv") %>% 
  drop_na() %>% distinct() %>% 
  inner_join(prism_list_to_zia %>% select(BROAD_ID, InChIKey), by= c("standard_inchi_key"= "InChIKey") ) %>% 
  distinct() %>% 
  rename(binding_score= interaction_strength, drug= BROAD_ID, target= gene_name ) %>%
  select(drug, target, binding_score) %>% 
  group_by(drug, target) %>% 
  summarise(binding_score= max(binding_score)) %>% 
  ungroup()


write_csv(ctrp_target_dtc,
          path = "~/cluster_scratch/prior/ctrpv2_target_dtc.csv")
write_csv(gdsc_target_dtc,
          path = "~/cluster_scratch/prior/gdsc_target_dtc.csv")
write_csv(prism_target_dtc,
          path = "~/cluster_scratch/prior/prism_target_dtc.csv")