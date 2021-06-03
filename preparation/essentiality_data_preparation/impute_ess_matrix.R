# the aim of this code is to impute ess matrix for four types of ess score
# in addition, I want to also try demeter and ceres original data so I will impute these
# data as well.



library(tidyverse)
setwd("/home/cloud-user/cluster_scratch/ces_io/")

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
set.seed(1234)

tmp1 <- impute_ess_df(ceres)
tmp2 <- impute_ess_df(demeter2)
tmp3 <- impute_ess_df(ces1)


write_csv(tmp1, "ceres_21q1_imputed.csv")
write_csv(tmp2, "demeter2_21q1_imputed.csv")
write_csv(tmp3, "ces1_21q1_imputed.csv")


