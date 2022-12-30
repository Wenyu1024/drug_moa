library(furrr)
library(future)
library(tidyverse)
setwd("/scratch/project_2003466/forward_modelling/noncancer_drugtarget_validation/")
load("null_dist_input.RData")

cor_res_rdm_perm_score3_prop0.25 <- get_cor_res_rdm(10, 3,0.025, cor_type= "cor_unsup" )

tictoc::tic()
plan(multisession,workers= 6)
par_grid_unsup <- expand_grid( 
  perm_num=1000,
  score_type=c(1,2,3),
  threshold_prop=c(1,0.05, 0.03, 0.025, 0.02, 0.01,0.005)) %>% 
  mutate(res= future_pmap(
    .l = list(perm_num, score_type, threshold_prop), 
    .f = get_cor_res_rdm,
    .options = furrr_options(seed=0000)
  )) 
plan(sequential)
tictoc::toc()
save(list = c("par_grid_unsup"),file = "null_dist_unsup.RData")
