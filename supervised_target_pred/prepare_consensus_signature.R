library(tidyverse)
drug_consensus <- read_csv( "~/cluster_scratch/L1000/consensus_signatures_challenge/consensus_signature_drugs.csv") %>% 
  rename(drug = X1)
load("~/cluster_scratch/forward_modelling/forwardmodelling_all.RData")
drug_consensus_ctrp <- drug_consensus %>% filter(drug %in% ctrp_data$broad_cpd_id)

drug_consensus_prism <- drug_consensus %>% filter(drug %in% prism_data$BROAD_ID)

#creating a cross reference table to connect gdsc id with broad id using
# drug repurposing hub

DRH <-  read_csv("/home/cloud-user/cluster_scratch/prior/drug_repurposing_hub_sample.csv") %>% 
  select(-vendor_name, -deprecated_broad_id, -smiles, -pert_iname) %>% 
  separate(col = broad_id,into = c("BRD", "id", NA, NA, NA),sep = "-", remove = T) %>% 
  unite(BRD, id, col = "BROAD_ID",  remove = T,sep= "-" ) %>% 
  distinct()%>% 
  arrange(BROAD_ID, pubchem_cid)

#In case there is two inchikeys only the first one is kept
DRH <- DRH[!duplicated(DRH$BROAD_ID),]

cross_ref_tibble <- DRH %>% 
  select(-InChIKey) %>% 
  mutate(pubchem_cid = as.character(pubchem_cid)) %>% 
  drop_na() %>% 
  distinct() %>% 
  inner_join(
    y = read_csv("~/cluster_scratch/prior/drug_id_list/gdsc2_idlist.csv") %>% 
      select(-smiles) %>%  drop_na() %>% distinct(), 
    by= c("pubchem_cid"= "PubCHEM") )

drug_consensus_gdsc <- cross_ref_tibble %>% 
  select(- pubchem_cid) %>% 
  inner_join(drug_consensus, by= c("BROAD_ID"= "drug")) %>% 
  select(-BROAD_ID) %>% 
  rename(drug= DRUG_ID)


write_csv(cross_ref_tibble,"~/cluster_scratch/prior/drug_id_list/cross_ref_tibble_gdsc")
setwd("~/cluster_scratch//L1000/consensus_signatures_challenge/")

write_csv(drug_consensus_ctrp, "feature_consensus_ctrp.csv")
write_csv(drug_consensus_gdsc,"feature_consensus_gdsc.csv")
write_csv(drug_consensus_prism,"feature_consensus_prism.csv")