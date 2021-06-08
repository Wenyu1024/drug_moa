# load("~/cluster_scratch/ctrpv2_data.RData")
# 
# ces1 <- readr::read_csv("~/cluster_scratch/impute_ess_20q4/ces1_478.csv")
# data <- ctrpv2_data %>%
#   ungroup() %>%  
#   slice(1:60)
# save.image("~/cluster_scratch/glmnet_modelling_cluster/input_test.RData")
# 

# load("~/cluster_scratch/ctrp/ctrpv2_data.RData")
setwd("~/cluster_scratch/ces_21q1_io")
ces1<- read_csv("./ces1_21q1_imputed.csv")
ceres  <- read_csv("./ceres_21q1_imputed.csv")
demeter2 <- read_csv("./demeter2_21q1_imputed.csv")
# exp_seq_pca <- read_csv("~/cluster_scratch/impute_and_derive_exp_based_feature_importance/exp_seq_pca.csv")
exp_seq <- read_csv("./expseq_21q1_imputed.csv")

# data <- ctrpv2_data %>%
#   ungroup() 
# rm(ctrpv2_data)
# save.image("~/cluster_scratch/forward_modelling/ctrp_input.RData")


load("~/cluster_scratch/gdsc/gdsc2.RData")
data <- gdsc_df %>%
  ungroup() 
rm(gdsc_df)
save.image("~/cluster_scratch/forward_modelling/gdsc_input.RData")

