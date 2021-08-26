# plot_target_feature_scatter <- function(geneleve_drug_representation_tibble, dtc_tibble, idx_included){ 
#   geneleve_drug_representation_tibble %>%
#     pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
#     inner_join(gdsc_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
#     filter(importance >0) %>%
#     # filter(binding_score >0) %>%
#     ggplot(aes(x= importance,y= binding_score)) +
#     geom_point(size= 0.01) +
#     theme_classic()
#   }
# 

####################################
# ctrp
#1 ces1 based feature importance

tmp <- feature_imp_ridge_ctrp_ces1 %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(ctrp_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value >0) %>%
  filter(binding_score >0) 

cor.test(tmp$binding_score, tmp$feature_value) # cor 0.21, p 3.306e-11  
# The higher the feature importance value based on CES1, the higher the binding score
tmp %>% 
  ggplot(aes(x= feature_value,y= binding_score)) +
  geom_point(size= 0.01) +
  theme_classic()

# 
tmp <- feature_imp_ridge_ctrp_ces1 %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(ctrp_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value < 0) %>%
  filter(binding_score >0) 

cor.test(tmp$binding_score, tmp$feature_value) # cor -0.03, p 0.07



#2 drug perturbated gene expression  (L1000)
tmp <- drug_consensus_ctrp %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(ctrp_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value <0) %>% #genes whose expression is downregulated
  filter(binding_score >0) 
cor.test(tmp$binding_score, tmp$feature_value) #-0.11, p = 0.008

#the more the gene is down-regulated, the higher the binding score

# tmp <- drug_consensus_ctrp %>%
#   pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
#   inner_join(ctrp_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
#   filter(feature_value >0) %>%
#   filter(binding_score >0) 
# cor.test(tmp$binding_score, tmp$feature_value) #-0.10, p = 0.008

###############################
# 3 gene expression based feature imp

tmp <- feature_imp_ridge_ctrp_exp %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(ctrp_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value <0) %>%  # feature<0 choose gene whose expression contribute positively to the sensitivity
  filter(binding_score >0) 

cor.test(tmp$binding_score, tmp$feature_value) # -0.12, p = 3.306e-11
# 


#gdsc
######################################################################
#1 ces1 based feature importance

tmp <- feature_imp_ridge_gdsc_ces1 %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(gdsc_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value >0) %>%
  filter(binding_score >0) 

cor.test(tmp$binding_score, tmp$feature_value) # cor 0.32, p 2.2e-16 
# The higher the feature importance value based on CES1, the higher the binding score
tmp %>% 
  ggplot(aes(x= feature_value,y= binding_score)) +
  geom_point(size= 0.01) +
  theme_classic()

# 
tmp <- feature_imp_ridge_gdsc_ces1 %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(gdsc_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value < 0) %>%
  filter(binding_score >0) 

cor.test(tmp$binding_score, tmp$feature_value) # cor -0.07, p 0.002



#2 drug perturbated gene expression  (L10000)
tmp <- drug_consensus_gdsc %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(gdsc_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value <0) %>% #genes whose expression is downregulated
  filter(binding_score >0) 
cor.test(tmp$binding_score, tmp$feature_value) #-0.11, p = 3.306e-11

#the more the gene is down-regulated, the higher the binding score

# tmp <- drug_consensus_gdsc %>%
#   pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
#   inner_join(gdsc_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
#   filter(feature_value >0) %>%
#   filter(binding_score >0) 
# cor.test(tmp$binding_score, tmp$feature_value) #-0.10, p = 0.008

# 3 gene expression based feature imp

tmp <- feature_imp_ridge_gdsc_exp %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(gdsc_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value <0) %>%  # feature<0 choose gene whose expression contribute positively to the sensitivity
  filter(binding_score >0) 

cor.test(tmp$binding_score, tmp$feature_value) # -0.15, p = 3.815e-09
# 


# PRISM
###########################
#1 ces1 based feature importance

tmp <- feature_imp_ridge_prism_ces1 %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(prism_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value >0) %>%
  filter(binding_score >0) 

cor.test(tmp$binding_score, tmp$feature_value) # cor 0.09
# The higher the feature importance value based on CES1, the higher the binding score
tmp %>% 
  ggplot(aes(x= feature_value,y= binding_score)) +
  geom_point(size= 0.01) +
  theme_classic()

# 
tmp <- feature_imp_ridge_prism_ces1 %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(prism_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value < 0) %>%
  filter(binding_score >0) 

cor.test(tmp$binding_score, tmp$feature_value) # cor 0.01, p 0.51


#2 drug perturbed gene expression  (L1000)
tmp <- drug_consensus_prism %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(prism_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value <0) %>% #genes whose expression is downregulated
  filter(binding_score >0) 
cor.test(tmp$binding_score, tmp$feature_value) #-0.12, p = 0.004

#the more the gene is down-regulated, the higher the binding score

# tmp <- drug_consensus_prism %>%
#   pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
#   inner_join(prism_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
#   filter(feature_value >0) %>%
#   filter(binding_score >0) 
# cor.test(tmp$binding_score, tmp$feature_value) #-0.10, p = 0.008

# 3 gene expression based feature imp

tmp <- feature_imp_ridge_prism_exp %>%
  pivot_longer(cols= -drug, names_to = "gene",values_to= "feature_value") %>%
  inner_join(prism_target_dtc ,by =c("drug"= "drug", "gene"="target")) %>%
  filter(feature_value <0) %>%  # feature<0 choose gene whose expression contribute positively to the sensitivity
  filter(binding_score >0) 

cor.test(tmp$binding_score, tmp$feature_value) # -0.03, p = 0.047




