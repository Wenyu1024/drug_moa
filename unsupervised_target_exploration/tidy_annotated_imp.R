library(tidyverse)
load("~/cluster_scratch/forward_modelling/targetpred_output_simplefiltering.RData")

#old
############################################################################################
#ctrp
CTRP_binary_PPIold_Consen_sig <- read_csv("~/cluster_scratch/prior/PPI/CTRP_binary_PPIold_Consen_sig.csv")%>% 
  pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
  mutate(dataset="CTRP", imptype= "ConSen-Sig") %>% 
  mutate_at(-1,as_factor) %>% 
  inner_join(y = feature_imp_ridge_ctrp_comb1 %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ),
             by=c("drug", "gene"))

CTRP_binary_PPIold_Conexp_sig <- read_csv("~/cluster_scratch/prior/PPI/CTRP_binary_PPIold_Conexp_sig.csv") %>% 
  pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
  mutate(dataset="CTRP", imptype= "ConExp-Sig") %>% 
  mutate_at(-1, as_factor) %>% 
  inner_join(drug_consensus_ctrp %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ))

data <- bind_rows(CTRP_binary_PPIold_Consen_sig,CTRP_binary_PPIold_Conexp_sig) %>% 
  mutate(anno_type1= 
           fct_recode(anno_type, 
                      Unknown= "NA_drug",
                      Unknown= "NA_gene",
                      Unknown= "NA_target_network",
                      Not_connected = "Not Connected",
                      "1st degree" = "1",
                      "2nd degree" = "2",
                      "3rd degree" = "3",
                      over_3= "4",
                      over_3= "5",
                      over_3= "6",
                      over_3= "7",
                      over_3= "8",
                      over_3= "9",
                      over_3= "10",
           )
  ) %>% 
  filter(drug %in% ctrp_drug_list) %>% 
  mutate(anno_type1= 
           fct_relevel(anno_type1, c("Target","1st degree","2nd degree","3rd degree","over_3","Not_connected", "Unknown" )))

write_csv(data, "~/cluster_scratch/prior/CTRP_binary_PPIold_data.csv")


#############################################################################
#GDSC
GDSC_binary_PPIold_Consen_sig <- read_csv("~/cluster_scratch/prior/PPI/GDSC_binary_PPIold_Consen_sig.csv")%>% 
  pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
  mutate(dataset="GDSC", imptype= "ConSen-Sig") %>% 
  mutate_at(-1,as_factor) %>% 
  inner_join(y = feature_imp_ridge_gdsc_comb1 %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ),
             by=c("drug", "gene"))

GDSC_binary_PPIold_Conexp_sig <- read_csv("~/cluster_scratch/prior/PPI/GDSC_binary_PPIold_Conexp_sig.csv") %>% 
  pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
  mutate(dataset="GDSC", imptype= "ConExp-Sig") %>% 
  mutate_at(-1,as_factor) %>% 
  inner_join(drug_consensus_gdsc %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ))

data <- bind_rows(GDSC_binary_PPIold_Consen_sig,GDSC_binary_PPIold_Conexp_sig) %>% 
  mutate(anno_type1= 
           fct_recode(anno_type, 
                      Unknown= "NA_drug",
                      Unknown= "NA_gene",
                      Unknown= "NA_target_network",
                      Not_connected = "Not Connected",
                      "1st degree" = "1",
                      "2nd degree" = "2",
                      "3rd degree" = "3",
                      over_3= "4",
                      over_3= "5",
                      over_3= "6",
                      over_3= "7",
                      over_3= "8",
                      over_3= "9"
                  
           )
  ) %>% 
  filter(drug %in% gdsc_drug_list) %>% 
  mutate(anno_type1= 
           fct_relevel(anno_type1, c("Target","1st degree","2nd degree","3rd degree","over_3","Not_connected", "Unknown" )))

write_csv(data, "~/cluster_scratch/prior/GDSC_binary_PPIold_data.csv")

###################################################################
#PRISM
PRISM_binary_PPIold_Consen_sig <- read_csv("~/cluster_scratch/prior/PPI/PRISM_binary_PPIold_Consen_sig.csv")%>% 
  pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
  mutate(dataset="PRISM", imptype= "ConSen-Sig") %>% 
  mutate_at(-1,as_factor) %>% 
  inner_join(y = feature_imp_ridge_prism_comb1 %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ),
             by=c("drug", "gene"))

PRISM_binary_PPIold_Conexp_sig <- read_csv("~/cluster_scratch/prior/PPI/PRISM_binary_PPIold_Conexp_sig.csv") %>% 
  pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
  mutate(dataset="PRISM", imptype= "ConExp-Sig") %>% 
  mutate_at(-1,as_factor) %>% 
  inner_join(drug_consensus_prism %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ))

data <- bind_rows(PRISM_binary_PPIold_Consen_sig,PRISM_binary_PPIold_Conexp_sig) %>% 
  mutate(anno_type1= 
           fct_recode(anno_type, 
                      Unknown= "NA_drug",
                      Unknown= "NA_gene",
                      Unknown= "NA_target_network",
                      Not_connected = "Not Connected",
                      "1st degree" = "1",
                      "2nd degree" = "2",
                      "3rd degree" = "3",
                      over_3= "4",
                      over_3= "5",
                      over_3= "6",
                      over_3= "7",
                      over_3= "8",
                      over_3= "9",
                      over_3= "10",
                      
           )
  ) %>% 
  filter(drug %in% prism_drug_list) %>% 
  mutate(anno_type1= 
           fct_relevel(anno_type1, c("Target","1st degree","2nd degree","3rd degree","over_3","Not_connected", "Unknown" )))

write_csv(data, "~/cluster_scratch/prior/PRISM_binary_PPIold_data.csv")
##################################################################################################
#new

############################################################################################
# #ctrp
# CTRP_binary_PPInew_Consen_sig <- read_csv("~/cluster_scratch/prior/PPI/CTRP_binary_PPInew_Consen_sig.csv")%>% 
#   pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
#   mutate(dataset="CTRP", imptype= "ConSen-Sig") %>% 
#   mutate_at(-1,as_factor) %>% 
#   inner_join(y = feature_imp_ridge_ctrp_comb1 %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ),
#              by=c("drug", "gene"))
# 
# CTRP_binary_PPInew_Conexp_sig <- read_csv("~/cluster_scratch/prior/PPI/CTRP_binary_PPInew_Conexp_sig.csv") %>% 
#   pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
#   mutate(dataset="CTRP", imptype= "ConExp-Sig") %>% 
#   mutate_at(-1, as_factor) %>% 
#   inner_join(drug_consensus_ctrp %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ))
# 
# 
# # function(x){
# #   w <- as.integer(x)
# #   case_when(
# #     is.na(w)~"Unknown"
# #     w>3~ "", 
# #             w<4&w>0 ~)
# # }
# 
# data <- bind_rows(CTRP_binary_PPInew_Consen_sig,CTRP_binary_PPInew_Conexp_sig) %>% 
#  mutate(anno_type1= 
#            fct_recode(anno_type, 
#                       Unknown= "NA_drug",
#                       Unknown= "NA_gene",
#                       Unknown= "NA_target_network",
#                       Not_connected = "Not Connected",
#                       "1st degree" = "1",
#                       "2nd degree" = "2",
#                       "3rd degree" = "3",
#                       over_3= "4",
#                       over_3= "5",
#                       over_3= "6",
#                       over_3= "7",
#                       over_3= "8",
#                       over_3= "9",
#                       over_3= "9",
#                       over_3= "10",
#                       over_3= "11",
#            )
#   ) %>% 
#   filter(drug %in% ctrp_drug_list) %>% 
#   mutate(anno_type1= 
#            fct_relevel(anno_type1, c("Target","1st degree","2nd degree","3rd degree","over_3","Not_connected", "Unknown" )))
# 
# write_csv(data, "~/cluster_scratch/prior/CTRP_binary_PPInew_data.csv")
# 
# #############################################################################
# #GDSC
# GDSC_binary_PPInew_Consen_sig <- read_csv("~/cluster_scratch/prior/PPI/GDSC_binary_PPInew_Consen_sig.csv")%>% 
#   pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
#   mutate(dataset="GDSC", imptype= "ConSen-Sig") %>% 
#   mutate_at(-1,as_factor) %>% 
#   inner_join(y = feature_imp_ridge_gdsc_comb1 %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ),
#              by=c("drug", "gene"))
# 
# GDSC_binary_PPInew_Conexp_sig <- read_csv("~/cluster_scratch/prior/PPI/GDSC_binary_PPInew_Conexp_sig.csv") %>% 
#   pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
#   mutate(dataset="GDSC", imptype= "ConExp-Sig") %>% 
#   mutate_at(-1,as_factor) %>% 
#   inner_join(drug_consensus_gdsc %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ))
# 
# data <- bind_rows(GDSC_binary_PPInew_Consen_sig,GDSC_binary_PPInew_Conexp_sig) %>% 
#   mutate(anno_type1= 
#            fct_recode(anno_type, 
#                       Unknown= "NA_drug",
#                       Unknown= "NA_gene",
#                       Unknown= "NA_target_network",
#                       Not_connected = "Not Connected",
#                       "1st degree" = "1",
#                       "2nd degree" = "2",
#                       "3rd degree" = "3",
#                       over_3= "4",
#                       over_3= "5",
#                       over_3= "6",
#                       over_3= "7",
#                       over_3= "8",
#                       over_3= "9",
#                       over_3= "10",
#                       over_3= "11"
#                       
#            )
#   ) %>% 
#   filter(drug %in% gdsc_drug_list) %>% 
#   mutate(anno_type1= 
#            fct_relevel(anno_type1, c("Target","1st degree","2nd degree","3rd degree","over_3","Not_connected", "Unknown" )))
# 
# write_csv(data, "~/cluster_scratch/prior/GDSC_binary_PPInew_data.csv")
# 
# ###################################################################
# #PRISM 
# PRISM_binary_PPInew_Consen_sig <- read_csv("~/cluster_scratch/prior/PPI/PRISM_binary_PPInew_Consen_sig.csv")%>% 
#   pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
#   mutate(dataset="PRISM", imptype= "ConSen-Sig") %>% 
#   mutate_at(-1,as_factor) %>% 
#   inner_join(y = feature_imp_ridge_prism_comb1 %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ),
#              by=c("drug", "gene"))
# 
# PRISM_binary_PPInew_Conexp_sig <- read_csv("~/cluster_scratch/prior/PPI/PRISM_binary_PPInew_Conexp_sig.csv") %>% 
#   pivot_longer(cols = -drug, names_to= "gene",values_to="anno_type" ) %>% 
#   mutate(dataset="PRISM", imptype= "ConExp-Sig") %>% 
#   mutate_at(-1,as_factor) %>% 
#   inner_join(drug_consensus_prism %>% pivot_longer(cols = -drug, names_to= "gene",values_to="imp" ))
# 
# data <- bind_rows(PRISM_binary_PPInew_Consen_sig,PRISM_binary_PPInew_Conexp_sig) %>% 
#   mutate(anno_type1= 
#            fct_recode(anno_type, 
#                       Unknown= "NA_drug",
#                       Unknown= "NA_gene",
#                       Unknown= "NA_target_network",
#                       Not_connected = "Not Connected",
#                       "1st degree" = "1",
#                       "2nd degree" = "2",
#                       "3rd degree" = "3",
#                       over_3= "4",
#                       over_3= "5",
#                       over_3= "6",
#                       over_3= "7",
#                       over_3= "8",
#                       over_3= "9",
#                       over_3= "10",
#                       over_3= "11",
#                       over_3= "12"
#                       
#            )
#   ) %>% 
#   filter(drug %in% prism_drug_list) %>% 
#   mutate(anno_type1= 
#            fct_relevel(anno_type1, c("Target","1st degree","2nd degree","3rd degree","over_3","Not_connected", "Unknown" )))
# 
# write_csv(data, "~/cluster_scratch/prior/PRISM_binary_PPInew_data.csv")
# 
