library(tidyverse)
library(tidymodels)

## 1 read in model based features
load("~/cluster_scratch/forward_modelling/forwardmodelling_all_new.RData")

#confirm whether the gene is matched
# tmp <- tibble(gene)
# tmp <- map(res_feature$ces1_perf,function(x){bind_cols(tmp, tibble(x))} )
# gene <- sort(colnames(ces1)[-1])
# transform_mat_to_tibble <- function(mat){as_tibble(mat, rownames= "drug")}


# ### 1.1 CTRP
# 
# ctrp_imp_ces1 <- ctrp_imp_ceres <-ctrp_imp_demeter2 <-ctrp_imp_exp <- matrix(nrow = 545,ncol = length(gene),dimnames = list(ctrp_data$broad_cpd_id, gene))
# 
# for (job_id in 1:545){
#   file_name <- paste0( '~/cluster_scratch/forward_modelling/ctrp_spearman_feature_imp/drug_' ,job_id,".RData" )
#   load(file_name)
#   ctrp_imp_ces1[job_id,] <- ces1_imp
#   ctrp_imp_ceres[job_id,] <- ceres_imp
#   ctrp_imp_demeter2[job_id,] <- demeter2_imp
#   ctrp_imp_exp[job_id,] <- exp_imp
# }
# 
# write_csv(transform_mat_to_tibble(ctrp_imp_ces1),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_ces1.csv")
# write_csv(transform_mat_to_tibble(ctrp_imp_ceres),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_ceres.csv")
# write_csv(transform_mat_to_tibble(ctrp_imp_demeter2),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_demeter2.csv")
# 
# 
# ### 1.2 GDSC
# 
# gdsc_imp_ces1 <- gdsc_imp_ceres <-gdsc_imp_demeter2 <-gdsc_imp_exp <- matrix(nrow = 198,ncol = length(gene),dimnames = list(gdsc_data$DRUG_ID, gene))
# 
# for (job_id in 1:198){
#   file_name <- paste0( '~/cluster_scratch/forward_modelling/gdsc_spearman_feature_imp/drug_' ,job_id,".RData" )
#   load(file_name)
#   gdsc_imp_ces1[job_id,] <- ces1_perf
#   gdsc_imp_ceres[job_id,] <- ceres_perf
#   gdsc_imp_demeter2[job_id,] <- demeter2_perf
#   gdsc_imp_exp[job_id,] <- exp_perf
# }
# 
# write_csv(transform_mat_to_tibble(gdsc_imp_ces1),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_ces1.csv")
# write_csv(transform_mat_to_tibble(gdsc_imp_ceres),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_ceres.csv")
# write_csv(transform_mat_to_tibble(gdsc_imp_demeter2),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_demeter2.csv")
# 
# 
# ### 1.3 PRISM
# 
# prism_imp_ces1 <- prism_imp_ceres <-prism_imp_demeter2 <-prism_imp_exp <- matrix(nrow = 1448,ncol = length(gene),dimnames = list(prism_data$BROAD_ID, gene))
# 
# for (job_id in 1:1448){
#   file_name <- paste0( '~/cluster_scratch/forward_modelling/prism_spearman_feature_imp/drug_' ,job_id,".RData" )
#   load(file_name)
#   prism_imp_ces1[job_id,] <- ces1_perf
#   prism_imp_ceres[job_id,] <- ceres_perf
#   prism_imp_demeter2[job_id,] <- demeter2_perf
#   prism_imp_exp[job_id,] <- exp_perf
# }
# 
# write_csv(transform_mat_to_tibble(prism_imp_ces1),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_ces1.csv")
# write_csv(transform_mat_to_tibble(prism_imp_ceres),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_ceres.csv")
# write_csv(transform_mat_to_tibble(prism_imp_demeter2),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_demeter2.csv")



feature_imp_ridge_ctrp_ces1 <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_ces1.csv")

feature_imp_ridge_ctrp_ceres <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_ceres.csv")

feature_imp_ridge_ctrp_demeter2 <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_demeter2.csv")

feature_imp_ridge_gdsc_ces1 <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_ces1.csv")

feature_imp_ridge_gdsc_ceres <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_ceres.csv")

feature_imp_ridge_gdsc_demeter2 <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_demeter2.csv")

feature_imp_ridge_prism_ces1 <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_ces1.csv")

feature_imp_ridge_prism_ceres <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_ceres.csv")

feature_imp_ridge_prism_demeter2 <- read_csv("~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_demeter2.csv")

## 2 combine the sensitivity signatures 
feature_imp_ridge_ctrp_comb1 <- 
  bind_cols(
    feature_imp_ridge_ctrp_ces1 %>% select(drug),
    as_tibble(
      (
        as.matrix(feature_imp_ridge_ctrp_ces1 %>% select(-drug)) +
          as.matrix(feature_imp_ridge_ctrp_ceres %>% select(-drug)) +
          as.matrix(feature_imp_ridge_ctrp_demeter2 %>% select(-drug))
      )/3
    )
  )

feature_imp_ridge_gdsc_comb1 <- 
  bind_cols(
    feature_imp_ridge_gdsc_ces1 %>% select(drug),
    as_tibble(
      (
        as.matrix(feature_imp_ridge_gdsc_ces1 %>% select(-drug)) +
          as.matrix(feature_imp_ridge_gdsc_ceres %>% select(-drug)) +
          as.matrix(feature_imp_ridge_gdsc_demeter2 %>% select(-drug)) 
      )/3
    )
  )

feature_imp_ridge_prism_comb1 <- 
  bind_cols(
    feature_imp_ridge_prism_ces1 %>% select(drug),
    as_tibble(
      (as.matrix(feature_imp_ridge_prism_ces1 %>% select(-drug)) +
         as.matrix(feature_imp_ridge_prism_ceres %>% select(-drug)) +
         as.matrix(feature_imp_ridge_prism_demeter2 %>% select(-drug))  
      )/3
    )
  )


## 3 read in target labels 
ctrp_target_binary <- read_csv("~/cluster_scratch/prior/ctrpv2_target_binary.csv")
gdsc_target_binary <- read_csv("~/cluster_scratch/prior/gdsc_target_binary.csv")
prism_target_binary <- read_csv("~/cluster_scratch/prior/prism_target_binary.csv")
ctrp_target_dtc <- read_csv("~/cluster_scratch/prior/ctrpv2_target_dtc.csv")
gdsc_target_dtc <- read_csv("~/cluster_scratch/prior/gdsc_target_dtc.csv")
prism_target_dtc <- read_csv("~/cluster_scratch/prior/ctrpv2_target_dtc.csv")

## 4 read in additional feature preparation
#fingerprint
setwd("~/cluster_scratch/fingerprints/new/fingerprint/result")
feature_maccs_ctrp <- read_csv("maccs_finger_ctrp.csv") %>% select(-smiles, -InChIKey,-pubchem_cid) %>% rename(drug =broad_cpd_id)
feature_extended_ctrp <- read_csv("extended_finger_ctrp.csv") %>% select(-smiles, -InChIKey,-pubchem_cid) %>% rename(drug =broad_cpd_id)

feature_maccs_gdsc <- read_csv("maccs_finger_gdsc.csv")%>% select(-smiles,-PubCHEM) %>% rename(drug =DRUG_ID) 
feature_maccs_gdsc <-feature_maccs_gdsc%>% distinct()
feature_extended_gdsc <- read_csv("extended_finger_gdsc.csv")%>% select(-smiles, -PubCHEM) %>% rename(drug =DRUG_ID)
feature_extended_gdsc <-feature_extended_gdsc%>% distinct()

feature_maccs_prism <- read_csv("maccs_finger_prism.csv")%>% select(-smiles, -InChIKey,-pubchem_cid) %>% rename(drug =BROAD_ID)
feature_extended_prism <- read_csv("extended_finger_prism.csv") %>% select(-smiles, -InChIKey,-pubchem_cid) %>% rename(drug =BROAD_ID)

## 5. read in consensus signature
setwd("~/cluster_scratch//L1000/consensus_signatures_challenge/")
drug_consensus_ctrp <- read_csv("feature_consensus_ctrp.csv")
drug_consensus_prism <- read_csv("feature_consensus_prism.csv")
drug_consensus_gdsc <- read_csv("feature_consensus_gdsc.csv")

##6 save workspace
rm(feature_imp_ridge_ctrp_ces1,feature_imp_ridge_ctrp_ceres,feature_imp_ridge_ctrp_demeter2,feature_imp_ridge_gdsc_ces1,feature_imp_ridge_gdsc_ceres,feature_imp_ridge_gdsc_demeter2,feature_imp_ridge_prism_ces1,feature_imp_ridge_prism_ceres,feature_imp_ridge_prism_demeter2)  
save.image("~/cluster_scratch/forward_modelling/targetpred_serverinput.RData")

##7 Output a Rdata for Rosanne to test additional modelling on the combined dataset
# combined_feature_tibble <- feature_imp_ridge_prism_comb1 %>% 
#   inner_join(feature_extended_prism) %>% 
#   inner_join(drug_consensus_prism %>% rename_at(-1,tolower)) %>% 
#   filter(drug %in% prism_drug_list)
# 
# combined_feature_mat <- combined_feature_tibble %>% 
#   arrange(drug) %>% 
#   select(-drug) %>% 
#   mutate_all(.funs = scale)
# 
# cor_mat_prism_spearman <- 
#   cor(t(combined_feature_mat),method = "spearman") %>% replace(is.na(.), 0)
# cor_mat_prism_spearman[cor_mat_prism_spearman< 0 ] <- 0 
# 
# row.names(cor_mat_prism_spearman) <- colnames(cor_mat_prism_spearman) <- combined_feature_tibble$drug
# 
# obj_list <- c("feature_imp_ridge_prism_comb1","feature_extended_prism", "drug_consensus_prism", "cor_mat_prism_spearman", "prism_target_binary", "prism_target_dtc", "prism_drug_list" )
# 




# Some previous analysis
# Write A function to calculate jaccard similarity based on fingerprint
# extended_cor <- feature_extended_finger_print %>% 
#   pivot_longer(cols = -broad_cpd_id,names_to= "name", values_to= "value") %>% 
#   pivot_wider(id_cols = name, names_from= broad_cpd_id,  values_from= value) %>%
#   arrange(name) %>% 
#   select(-name) %>% 
#   as.matrix() %>%
#   prabclus::jaccard()
# cor(method = "spearman")


# hist(extended_cor) 
# extended_cor1 <- (extended_cor - 0.5)* (-2)
# extended_cor1 <- -extended_cor + mean(extended_cor)
# hist(extended_cor1) 

# maccs_cor <- feature_maccs_finger_print %>% 
#   pivot_longer(cols = -broad_cpd_id,names_to= "name", values_to= "value") %>% 
#   arrange(broad_cpd_id) %>% 
#   pivot_wider(id_cols = name, names_from= broad_cpd_id,  values_from= value) %>% 
#   select(-name) %>% 
#   as.matrix() %>% 
#   prabclus::jaccard()
#   # cor(method = "spearman")

# hist(maccs_cor)
# maccs_cor1 <- -maccs_cor + mean(maccs_cor)
# hist(maccs_cor1) 
# hist(sen_cor)

# ces1_cor <- cor(t(
#   x = feature_imp_ridge_ctrp_ces1 %>% 
#     rename(drug=broad_cpd_id) %>% 
#     arrange(drug) %>% 
#     select(-drug) %>% 
#     mutate_all(.funs = scale)), method = "spearman")
# ces1_cor[ces1_cor<0] <- 0
# colnames(ces1_cor) <- row.names(ces1_cor) <- sort(feature_imp_ridge_ctrp_ces1$broad_cpd_id)
# 
# sen_cor <- cor(t(
#   x = feature_sen_ctrp %>% 
#     rename(drug=broad_cpd_id) %>% 
#     arrange(drug) %>% 
#     select(-drug) %>% 
#     mutate_all(.funs = scale)), method = "spearman")
# sen_cor[sen_cor<0] <- 0
# colnames(sen_cor) <- row.names(sen_cor) <- sort(feature_imp_ridge_ctrp_ces1$broad_cpd_id)


# ces1_extended_cor <- (ces1_cor + extended_cor)/2
# ces1_extended_sen_cor <- (ces1_cor + extended_cor+ sen_cor)/3
# ces1_maccs_cor <- (ces1_cor + maccs_cor)/2
# ces1_maccs_sen_cor <- (ces1_cor + maccs_cor+ sen_cor)/3

# save.image("~/cluster_scratch/glmnet_modelling/target_pred_server/targetpred_serverinput.RData")



