setwd("/scratch/project_2003466/ces_21q1_io/")
library(tidyverse)
data <- read_csv("ces_input_21q1.csv")

lm <- lm(ceres~ demeter2+mut+exp_seq+cn+exp_array+DepMap_ID ,data = data)
pred <- predict.lm(lm)

save.image("tmp.RData")
load("tmp.RData")
# tmp = proc.time() - ptm
# print(tmp)
data <- data %>% 
  select_at(c(1:4,6)) %>% 
  mutate(ces1= pred) 

# Only include genes whose observation over 3/4  
gene_list <- data %>% 
  group_by(gene) %>% 
  summarize(n_obs= n()) %>% 
  ungroup() %>% 
  filter(n_obs > length(unique(data$DepMap_ID))*0.75 ) %>% 
  pull(gene)

data <- data %>% filter(gene %in% gene_list)
ces1_437 <- data %>% select(DepMap_ID, gene, ces1) %>% pivot_wider(id_cols = DepMap_ID, names_from= gene, values_from= ces1)
ceres_437 <- data %>% select(DepMap_ID, gene, ceres) %>% pivot_wider(id_cols = DepMap_ID, names_from= gene, values_from= ceres)
demeter2_437 <- data %>% select(DepMap_ID, gene, demeter2) %>% pivot_wider(id_cols = DepMap_ID, names_from= gene, values_from= demeter2)
expseq_437 <- data %>% select(DepMap_ID, gene, exp_seq ) %>% pivot_wider(id_cols = DepMap_ID, names_from= gene, values_from= exp_seq)

impute_ess_df <- 
  function(df){
    library(missForest)
    mat <- as.matrix(df[,-1])
    doParallel::registerDoParallel(cores = 39) # set based on number of CPU cores
    doRNG::registerDoRNG(seed = 1234)
    imputed_X <- missForest(xmis = mat, maxiter = 10, verbose = T, parallelize = 'variables')$ximp
    df1 <- bind_cols(df[,1], as.data.frame(imputed_X))
    return(df1)
  }

save.image("tmp1.RData")
set.seed(1234)
ces1_437 <- impute_ess_df(ces1_437)
write_csv(ces1_437, "ces1_21q1_imputed.csv")
ceres_437 <- impute_ess_df(ceres_437)
write_csv(ceres_437, "ceres_21q1_imputed.csv")
demeter2_437 <- impute_ess_df(demeter2_437)
write_csv(demeter2_437, "demeter2_21q1_imputed.csv")
expseq_437 <-  impute_ess_df(expseq_437)
write_csv(expseq_437, "expseq_21q1_imputed.csv")
save.image("tmp2.RData")
