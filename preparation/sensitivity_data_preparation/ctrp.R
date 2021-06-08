library(tidyverse)
library(data.table)
setwd("~/cluster_scratch/ctrp")

# download.file(
#   url = "https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip",
#   method = "wget" )

ctrpv2_cell_ref <- fread("v20.meta.per_cell_line.txt")
# depmap_cell_ref <- fread("https://ndownloader.figshare.com/files/21522000")
depmap_cell_ref <- read_csv("~/cluster_scratch/depmap_21q1/sampleinfo")
cross_cell_depmap_ctrpv2 <- full_join(ctrpv2_cell_ref,depmap_cell_ref, by= c("ccl_name" = "stripped_cell_line_name"))

ctrpv2_sensi<- fread("v20.data.curves_post_qc.txt")

# focus only on the inhibitor and drugs with target info from drugbank
# there is no need to apply to such filters.
ctrpv2_compound <- data.table::fread("v20.meta.per_compound.txt") %>% 
  # mutate(inhibitor= str_deta)
  # filter(str_detect(string = target_or_activity_of_compound, pattern = "inhibitor")) %>% 
  # filter(gene_symbol_of_protein_target != "") %>% 
  select(master_cpd_id, cpd_name, broad_cpd_id, gene_symbol_of_protein_target, target_or_activity_of_compound,cpd_smiles) %>% 
  mutate(drug_bank_target_avaliability= !(gene_symbol_of_protein_target %in%  "") ) %>%  
  separate_rows(gene_symbol_of_protein_target,sep = ";") %>% 
  group_by_at(.vars= c(1:3,5:6)) %>% 
  nest(target = gene_symbol_of_protein_target)

#clean the drug sensitivity data so that the compound and cell line information are straight forward (to remove the further use of ctrpv2_experiemnt tibble)
ctrpv2_experiment <-  fread("v20.meta.per_experiment.txt")

#add the ccl id to the ctrpv2_sensi data and remove unnecessary columns
ctrpv2_sensi <- ctrpv2_sensi %>% 
  left_join(ctrpv2_experiment %>% select(experiment_id, master_ccl_id)) %>% 
  select(master_ccl_id,master_cpd_id,area_under_curve,apparent_ec50_umol) 

rm(ctrpv2_experiment)

sensitivity <- ctrpv2_sensi %>% 
  left_join(cross_cell_depmap_ctrpv2 %>% select(master_ccl_id,DepMap_ID)) %>% 
  select(-master_ccl_id)

ctrpv2_data <- ctrpv2_compound %>% 
  select(master_cpd_id,cpd_name, broad_cpd_id,drug_bank_target_avaliability, target) %>% 
  nest_join(sensitivity) 

rm(sensitivity)

# add indicator columns
# add columns to indicate the number of obs available for reverse and forward modelling
predictor_df <- read_csv("~/cluster_scratch/ces_21q1_io/ces1_21q1_imputed.csv")

ctrpv2_data <- ctrpv2_data %>%
  mutate(drug_type= case_when(  str_detect(cpd_smiles, "\\.")~ "drug pair",
                                TRUE ~ "single drug")) %>% 
  mutate(n_forward_modelling= map_int(.x=sensitivity, 
                              .f=function(x){sum(x %>% pull(DepMap_ID) %in% predictor_df$DepMap_ID)}))
  

save(ctrpv2_data, file = "ctrpv2_data.RData")
