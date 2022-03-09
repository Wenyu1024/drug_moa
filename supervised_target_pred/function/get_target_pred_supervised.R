source("~/cluster_wrk/drug_moa/supervised_target_pred/function/get_target_matrix_from_long_target_tibble.R")
source("~/cluster_wrk/drug_moa/supervised_target_pred/function/no_tunning_weighted_averaging.R")

get_predictions <- function(sig_data,train_drugs, pred_drugs, train_label, similarity= "spearman") {
  train_drugs <- sort(train_drugs)
  pred_drugs <- sort(pred_drugs)
  sig_data <- sig_data %>% 
    filter(drug %in% union(train_drugs, pred_drugs)) %>% 
    arrange(drug)
  # feature_name <- function(x){deparse(substitute(sig_data))}
  # istahimoto = str_detect(string = "feature_name", pattern = "extended")| str_detect(string = "feature_name", pattern = "maccs")
  if (similarity== "tahimoto"){
     train_cor_mat <- sig_data %>% 
      arrange(drug) %>% 
      select(-drug) 
     train_cor_mat = 1- as.matrix(vegan::vegdist(x = ( train_cor_mat),method = "jaccard",upper = F))
    
  } else {
     train_cor_mat <- sig_data %>% 
      arrange(drug) %>% 
      select(-drug) %>% 
      mutate_all(.funs = scale)
     train_cor_mat = cor(t( train_cor_mat),method = similarity) %>% replace(is.na(.), 0)
     train_cor_mat[train_cor_mat< 0 ] <- 0 
  }
  colnames(train_cor_mat) <- sig_data$drug
  row.names(train_cor_mat) <- sig_data$drug
  diag(train_cor_mat) <- 0
  target_mat <- get_target_mat(
    target_tibble = train_label %>%
      filter(drug %in% train_drugs) %>%
      arrange(drug)
  )
  train_idx <- match(x = train_drugs, table = sig_data$drug, nomatch = NA)
  test_idx <- match(x = pred_drugs, table = sig_data$drug, nomatch = NA)
  
  pred <- no_tunning_weighted_averaging(
    target_mat = target_mat,
    cor_test = train_cor_mat[,train_idx],
    pred_new = T
    )
  colnames(pred) <- pred_drugs
  return(pred)
}


