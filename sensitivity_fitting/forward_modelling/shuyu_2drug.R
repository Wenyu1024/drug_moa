library(tidyverse)
load("~/cluster_scratch/forward_modelling/ctrp_input.RData")
tmp_ctrp <- data %>% filter(cpd_name %in% c("carboplatin","olaparib"))
load("~/cluster_scratch/forward_modelling/gdsc_input.RData")
tmp_gdsc <- data %>% filter(DRUG_NAME == "Olaparib")
load("~/cluster_scratch/forward_modelling/prism_input.RData")
tmp_prism <- data %>% filter(name %in% c("carboplatin","olaparib"))
rm(data)

setwd("~/cluster_scratch/depmap_21q1/")
cell_info <- read_csv("sampleinfo")


save(list=c("tmp_ctrp","tmp_gdsc","tmp_prism", "cell_info"), file = "shuyu_2drug.RData")
