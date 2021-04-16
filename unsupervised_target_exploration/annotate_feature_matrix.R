# The function takes three input and 
# output a long table of feature importance matrix with target network annotation:
# the function may take several minute to run, depend on how many rows of the input are there.
library(tidyverse)
library(tidymodels)
annotate_feature_matrix <- function(network, feature_imp_mat, target){
  network_full <- network %>% 
    rename(gene1= gene2, gene2= gene1) %>% 
    bind_rows(network)
  target_gene <- target$target_gene
  drug_list <- sort(unique(feature_imp_mat$drug))
  gene_list <-  colnames(feature_imp_mat)[-1]
  target_list <- NULL
  i=1
  res_list <- as.list(1:length(drug_list))
  # lets use a list column to assign the result and withdraw it later
  for (drug_i in drug_list){
    target_list <- NULL
    # tic()
    target_i <-  target %>% filter(drug %in%  drug_i) %>% pull(target_gene)
    first_n <- setdiff(network_full %>% filter( gene1 %in% target_i) %>% pull(gene2) , target_i)
    first_n <- c(first_n,"NA")
    second_n <- setdiff(network_full %>% filter( gene1 %in% first_n) %>% pull(gene2) , union(target_i,first_n))
    second_n <- c(second_n,"NA")
    third_n <- setdiff(
      x = network_full %>% filter(gene1 %in% second_n) %>% pull(gene2),
      y = union(x = union(target_i, first_n),y = second_n)
    )
    third_n <- c(third_n,"NA")
    for (gene in  gene_list ) { if 
      (gene %in%  target_i)   {target_cat <- "0"} else if 
      (gene %in%  first_n)  {target_cat <- "1"} else if
      (gene %in%  second_n) {target_cat <- "2"} else if 
      (gene %in%  third_n)  {target_cat <- "3"} else {target_cat <- "4"}
      target_list <- c(target_list,target_cat)}  
    res_list[i] <- list(target_list)
    print(i)
    # check how to change this print out to a progress bar and think how does input size
    # impacting this
    # https://github.com/r-lib/progress
    i=i+1
  }
  
  target_category <- expand.grid( gene_list,drug_list)
  target_category <- cbind(target_category,unlist(res_list))
  colnames(target_category) <- c( "gene","drug", "target_cat")
  target_category <- as.data.frame(target_category) %>% as_tibble()
  
  res <- feature_imp_mat%>% 
    pivot_longer(cols= -drug, names_to="gene",values_to="imp") %>% 
    inner_join(target_category,by=c("drug", "gene"))
  
  return(res)
}
