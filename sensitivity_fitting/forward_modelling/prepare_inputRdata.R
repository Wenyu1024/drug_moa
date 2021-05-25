# load("~/cluster_scratch/ctrpv2_data.RData")
# 
# ces1 <- readr::read_csv("~/cluster_scratch/impute_ess_20q4/ces1_478.csv")
# data <- ctrpv2_data %>%
#   ungroup() %>%  
#   slice(1:60)
# save.image("~/cluster_scratch/glmnet_modelling_cluster/input_test.RData")
# 

load("~/cluster_scratch/ctrpv2_data.RData")
ces1<- read_csv("~/cluster_scratch/impute_ess_20q4/ces1_478.csv")
ces2 <- read_csv("~/cluster_scratch/impute_ess_20q4/ces2_478.csv")
ceres  <- read_csv("~/cluster_scratch/impute_ess_20q4/ceres_478.csv")
demeter2 <- read_csv("~/cluster_scratch/impute_ess_20q4/demeter2_478.csv")
exp_seq_pca <- read_csv("~/cluster_scratch/impute_and_derive_exp_based_feature_importance/exp_seq_pca.csv")

data <- ctrpv2_data %>%
  ungroup() 
rm(ctrpv2_data)
save.image("~/cluster_scratch/glmnet_modelling_cluster/input.RData")
