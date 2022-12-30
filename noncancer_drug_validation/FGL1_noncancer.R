library(tidyverse)
# prism_unsupervised <- read_csv("./fig_res/prism312_unsupervised.csv")
# load("~/cluster_scratch/forward_modelling/targetpred_output_simplefiltering.RData")
prism_unsupervised_noncancer_FGL1 <- prism_unsupervised %>% 
  pivot_longer(cols=-target, names_to="drug",values_to="sig_value") %>% 
  group_by(drug) %>% 
  slice_max(order_by= sig_value,n = 100) %>% 
  ungroup() %>% 
  filter(target=="FGL1") 

prism_unsupervised_noncancer_FOXO4 <- prism_unsupervised %>% 
  pivot_longer(cols=-target, names_to="drug",values_to="sig_value") %>% 
  group_by(drug) %>% 
  slice_min(order_by= sig_value,n = 100) %>%
  ungroup() %>% 
  filter(target=="FOXO4")
