# Earlier calculation were based on annotated_long_df and were limited to drugs with both 
# ConSen-Sig and ConExp-Sig.

# Here for exploring the application potential of ConSen-Sig persa and independantly for each drug,
# we recalculated the auc for all the drugs.
library(furrr)
library(tidyverse)
library(tictoc)
# set.seed(0)
load("~/cluster_scratch/forward_modelling/targetpred_output_simplefiltering.RData")
get_auc_for_imp_tibble <- function(imp_tibble=feature_imp_ridge_prism_comb1,
                                   target_tibble=prism_target_binary){
    
  gene_list <- colnames(imp_tibble)[-1]
  drug_list <- unique(target_tibble$drug)
  get_auc_for_drug_i <- function(i){
    drug_i = imp_tibble$drug[i]
    pred_vec <- imp_tibble[i,-1]
    pred_vec <- ifelse(test = pred_vec > 0, yes = pred_vec, no= 0 )
    pred_vec <- unlist(pred_vec)
  
    if (!(drug_i %in% drug_list)){return(NA)} else{
      true_target_list <-  (target_tibble %>% filter(drug== drug_i) %>% pull(target))
      if (length(intersect(gene_list, true_target_list))==0){return(NA)} else{
        truth_vec <- factor(x = (gene_list %in% true_target_list), levels = c(TRUE,FALSE))
        auc= yardstick::roc_auc_vec(truth = truth_vec,estimate = pred_vec)
        return(auc)
      }
    }
    }
  auc_vec <- future_map_dbl(.x = 1:nrow(imp_tibble),
                     .f = get_auc_for_drug_i
                     , .options = future_options(seed = 0)
                     )
  return(auc_vec)
}

# debug(get_auc_for_imp_tibble)
tic()
plan(multisession, workers=7)

ctrp_unsupervised_auc <-
  get_auc_for_imp_tibble(imp_tibble = feature_imp_ridge_ctrp_comb1, 
                         target_tibble = ctrp_target_binary)
ctrp_unsupervised_auc <- tibble(drug= feature_imp_ridge_ctrp_comb1$drug, 
                                auc= ctrp_unsupervised_auc)

gdsc_unsupervised_auc <-
  get_auc_for_imp_tibble(imp_tibble = feature_imp_ridge_gdsc_comb1,
                         target_tibble = gdsc_target_binary)

gdsc_unsupervised_auc <- tibble(drug= feature_imp_ridge_gdsc_comb1$drug,
                                auc= gdsc_unsupervised_auc)


prism_unsupervised_auc <-
  get_auc_for_imp_tibble(imp_tibble = feature_imp_ridge_prism_comb1,
                       target_tibble = prism_target_binary)

prism_unsupervised_auc <- tibble(drug= feature_imp_ridge_prism_comb1$drug,
                                auc= prism_unsupervised_auc)

# plan(sequential)
toc()

save(list = c("ctrp_unsupervised_auc", "gdsc_unsupervised_auc", "prism_unsupervised_auc"),
     file = "~/cluster_scratch/unsupervised_target_pred/unsupervised_target_pred_ConSenSig.RData")
# save(list = c("ctrp_unsupervised_auc"),
#      file = "~/cluster_scratch/unsupervised_target_pred/unsupervised_target_pred_ConSenSig.RData")


# Why the value is different? floating point difference???
# load("~/cluster_scratch/unsupervised_target_pred/unsupervised_target_pred_method_comparison.RData")
# auc_A <- ctrp_unsupervised_auc
# load("~/cluster_scratch/unsupervised_target_pred/unsupervised_target_pred_ConSenSig.RData")
# auc_B <- ctrp_unsupervised_auc
# 
# tmp <- auc_A %>% filter(imptype=="ConSen-Sig") %>% inner_join(auc_B)
# 
# hist((tmp$AUC- tmp$auc))
# # https://stackoverflow.com/questions/5683533/parallelism-subtly-different-floating-point-results
# #	BRD-A35588707 	0.6117881 / 0.6152434
# # BRD-K06593056 ConSen-Sig 0.2590504 0.2593674



