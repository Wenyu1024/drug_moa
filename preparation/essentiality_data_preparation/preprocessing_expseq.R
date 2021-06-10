# The Exp matrix is sparse, leading to glmnet model failure.
# Lets use a PCA first to do a dimension reduction.

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

exp_seq_imputed <- read_csv("~/cluster_scratch/ces_21q1_io/expseq_21q1_imputed.csv")

tmp <- exp_seq_imputed[,2:10624]
tmp <- prcomp(tmp)
tmp$sdev
transform_dev_to_percent_acc <- function(x){
  x_acc <- NULL
  for (i in 1:length(x)){
    x_acc[i] <- sum(x[1:i])/sum(x)
    # print(i)
    # print(x_acc[i])
  }
  return(x_acc)
}
tmp1 <- transform_dev_to_percent_acc(tmp$sdev)
tmp1 <- tmp$x
sum(transform_dev_to_percent_acc(tmp$sdev)<0.95) # 174
exp_seq_pca <- exp_seq_imputed %>% select(DepMap_ID) %>% bind_cols( as.data.frame(tmp1[,1:1203]))
write_csv(exp_seq_pca, path = "~/cluster_scratch/impute_and_derive_exp_based_feature_importance/exp_seq_pca.csv")
