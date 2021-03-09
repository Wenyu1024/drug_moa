source('/projappl/project_2003466/function/get_fixed_glmnet.R')
source('/projappl/project_2003466/function/get_feature_imp.R')
load("/scratch/project_2003466/glmnet_modelling_cluster/input.RData")
library(tidyverse, quietly = T)
library(furrr, quietly = T)
library(snow, quietly = T)
print("start")

cl <- getMPIcluster()
plan(cluster, workers = cl)

res_feature <- data %>% 
  mutate(ces1_perf= future_map(.x = sensitivity, .progress = FALSE,
                               .f = ~get_feature_imp(sen_df = .x,predictor_df = ces1),
                               .options = furrr_options(seed = 0000))) %>% 
  mutate(ces2_perf= future_map(.x = sensitivity, .progress = FALSE,
                               .f = ~get_feature_imp(sen_df = .x,predictor_df = ces2),
                               .options = furrr_options(seed = 0000))) %>%
  mutate(ceres_perf= future_map(.x = sensitivity, .progress = FALSE,
                                .f = ~get_feature_imp(sen_df = .x,predictor_df = ceres),
                                .options = furrr_options(seed = 0000))) %>%
  mutate(demeter2_perf= future_map(.x = sensitivity, .progress = FALSE,
                                   .f = ~get_feature_imp(sen_df = .x,predictor_df = demeter2),
                                   .options = furrr_options(seed = 0000)))  %>%
  mutate(exp_perf= future_map(.x = sensitivity, .progress = FALSE,
                                   .f = ~get_feature_imp(sen_df = .x,predictor_df = exp_seq_pca),
                                   .options = furrr_options(seed = 0000)))  %>%
  select(broad_cpd_id, ces1_perf, ces2_perf, ceres_perf, demeter2_perf,exp_perf)

stopCluster(cl)
print("success")
save(res_feature, file = "/scratch/project_2003466/glmnet_modelling_cluster/ridge_feature_multinode2.RData")

