library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
setwd("~/cluster_scratch/ctrp")
ctrpv2_cell_ref <- fread("v20.meta.per_cell_line.txt")
# depmap_cell_ref <- fread("https://ndownloader.figshare.com/files/21522000")
depmap_cell_ref <- fread("https://ndownloader.figshare.com/files/26261569")
cross_cell_depmap_ctrpv2 <- full_join(ctrpv2_cell_ref,depmap_cell_ref, by= c("ccl_name" = "stripped_cell_line_name"))
ctrpv2_sensi<- fread("v20.data.curves_post_qc.txt")

# focus only on the inhibitor and drugs with targe info from drugbank
ctrpv2_compound <- data.table::fread("v20.meta.per_compound.txt") %>% 
  # mutate(inhibitor= str_deta)
  filter(str_detect(string = target_or_activity_of_compound, pattern = "inhibitor")) %>% 
  filter(gene_symbol_of_protein_target != "") %>% 
  select(master_cpd_id, cpd_name, broad_cpd_id, gene_symbol_of_protein_target, target_or_activity_of_compound,cpd_smiles) %>% 
  separate_rows(gene_symbol_of_protein_target,sep = ";") %>% 
  group_by_at(.vars= c(1:3,5:6)) %>% 
  nest(target= gene_symbol_of_protein_target)

#clean the drug sensitivity data so that the compound and cell line information are straight forward (to remove the further use of ctrpv2_experiemnt tibble)
ctrpv2_experiment <-  fread("v20.meta.per_experiment.txt")

#add the ccl id to the ctrpv2_sensi data and remove unnecessary columns
ctrpv2_sensi <- ctrpv2_sensi %>% 
  left_join(ctrpv2_experiment %>% select(experiment_id, master_ccl_id)) %>% 
  select(master_ccl_id,master_cpd_id,area_under_curve,apparent_ec50_umol) 

rm(ctrpv2_experiment)

sensitivity <- ctrpv2_sensi %>% left_join(cross_cell_depmap_ctrpv2 %>% select(master_ccl_id,DepMap_ID))

ctrpv2_data <- ctrpv2_compound %>% 
  select(master_cpd_id,cpd_name, broad_cpd_id,target) %>% 
  nest_join(sensitivity) 
rm(sensitivity)

save(ctrpv2_data, file = "ctrpv2_data.RData")
# comment: I should make target and sensitivity dataset independent, 
# as I am also using the DTC as a secondary source. (DTC is even going to be a prelimnary source, if the result is better)