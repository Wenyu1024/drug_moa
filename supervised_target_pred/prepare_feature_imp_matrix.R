library(tidyverse)
load("~/cluster_scratch/forward_modelling/forwardmodelling_all.RData")
transform_mat_to_tibble <- function(mat){as_tibble(mat, rownames= "drug")}
#confirm whether the gene is matched
# tmp <- tibble(gene)
# tmp <- map(res_feature$ces1_perf,function(x){bind_cols(tmp, tibble(x))} )
gene <- sort(colnames(ces1)[-1])



### 1.1 CTRP

ctrp_imp_ces1 <- ctrp_imp_ceres <-ctrp_imp_demeter2 <-ctrp_imp_exp <- matrix(nrow = 545,ncol = length(gene),dimnames = list(ctrp_data$broad_cpd_id, gene))

for (job_id in 1:545){
  file_name <- paste0( '~/cluster_scratch/forward_modelling/ctrp_spearman_feature_imp/drug_' ,job_id,".RData" )
  load(file_name)
  ctrp_imp_ces1[job_id,] <- ces1_imp
  ctrp_imp_ceres[job_id,] <- ceres_imp
  ctrp_imp_demeter2[job_id,] <- demeter2_imp
  ctrp_imp_exp[job_id,] <- exp_imp
}

write_csv(transform_mat_to_tibble(ctrp_imp_ces1),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_ces1.csv")
write_csv(transform_mat_to_tibble(ctrp_imp_ceres),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_ceres.csv")
write_csv(transform_mat_to_tibble(ctrp_imp_demeter2),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_demeter2.csv")
write_csv(transform_mat_to_tibble(ctrp_imp_exp),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_ctrp_exp.csv")


### 1.2 GDSC

gdsc_imp_ces1 <- gdsc_imp_ceres <-gdsc_imp_demeter2 <-gdsc_imp_exp <- matrix(nrow = 198,ncol = length(gene),dimnames = list(gdsc_data$DRUG_ID, gene))

for (job_id in 1:198){
  file_name <- paste0( '~/cluster_scratch/forward_modelling/gdsc_spearman_feature_imp/drug_' ,job_id,".RData" )
  load(file_name)
  gdsc_imp_ces1[job_id,] <- ces1_imp
  gdsc_imp_ceres[job_id,] <- ceres_imp
  gdsc_imp_demeter2[job_id,] <- demeter2_imp
  gdsc_imp_exp[job_id,] <- exp_imp
}

write_csv(transform_mat_to_tibble(gdsc_imp_ces1),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_ces1.csv")
write_csv(transform_mat_to_tibble(gdsc_imp_ceres),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_ceres.csv")
write_csv(transform_mat_to_tibble(gdsc_imp_demeter2),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_demeter2.csv")
write_csv(transform_mat_to_tibble(gdsc_imp_exp),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_gdsc_exp.csv")


### 1.3 PRISM

prism_imp_ces1 <- prism_imp_ceres <-prism_imp_demeter2 <-prism_imp_exp <- matrix(nrow = 1448,ncol = length(gene),dimnames = list(prism_data$BROAD_ID, gene))

for (job_id in 1:1448){
  file_name <- paste0( '~/cluster_scratch/forward_modelling/prism_spearman_feature_imp/drug_' ,job_id,".RData" )
  load(file_name)
  prism_imp_ces1[job_id,] <- ces1_imp
  prism_imp_ceres[job_id,] <- ceres_imp
  prism_imp_demeter2[job_id,] <- demeter2_imp
  prism_imp_exp[job_id,] <- exp_imp
}

write_csv(transform_mat_to_tibble(prism_imp_ces1),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_ces1.csv")
write_csv(transform_mat_to_tibble(prism_imp_ceres),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_ceres.csv")
write_csv(transform_mat_to_tibble(prism_imp_demeter2),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_demeter2.csv")
write_csv(transform_mat_to_tibble(prism_imp_exp),"~/cluster_scratch/forward_modelling/feature_imp/feature_imp_ridge_prism_exp.csv")


