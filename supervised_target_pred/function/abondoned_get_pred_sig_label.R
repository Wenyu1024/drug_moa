# library(tidyverse,quietly = T)
# get_predictions <- function(sig_data,train_drugs, pred_drugs, train_label) {
#   sig_data <- sig_data %>% 
#     filter(drug %in% union(train_drugs, pred_drugs))
#   train_cor_mat <- cor(t(sig_data %>% select(-drug)))
#   colnames(train_cor_mat) <- sig_data$drug
#   row.names(train_cor_mat) <- sig_data$drug
#   diag(train_cor_mat) <- 0
#   target_mat <- get_target_mat(
#     target_tibble = train_label %>%
#       filter(drug %in% train_drugs) %>%
#       arrange(drug)
#   )
#   train_idx <- match(x = train_drugs, table = sig_data$drug, nomatch = NA)
#   test_idx <- match(x = pred_drugs, table = sig_data$drug, nomatch = NA)
#   
#   pred <- no_tunning_weighted_averaging(
#     target_mat = target_mat,
#     cor_mat = train_cor_mat,
#     cor_test = train_cor_mat[,train_idx],
#     test_idx = test_idx,
#     pred_new = T)
#   min_max_normalization <- function(x){(x-min(x))/ (max(x)-min(x)) }
#   pred_normalized= apply(pred,1 ,min_max_normalization)
#   return(pred)
# }
