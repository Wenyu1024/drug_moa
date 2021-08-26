# The function takes three input and 
# output a long table of feature importance matrix with target network annotation:
# the function may take several minute to run, depend on how many rows of the input are there.
library(tidyverse)
library(igraph)
library(tictoc)
library(furrr)

load("/scratch/project_2003466/forward_modelling/targetpred_output_simplefiltering.RData")
setwd("/scratch/project_2003466/prior/PPI")

# load("~/cluster_scratch/forward_modelling/targetpred_output_simplefiltering.RData")
# setwd("~/cluster_scratch/prior/PPI/")



# kegg_tibble <- read_csv("kegg_tibble.csv")
PPI_new <- read_csv("PPI_alllinks_simplified.csv") 
PPI_old <- read_csv( "PPI_old_simplified.csv" ) 



# annotate_feature_matrix <- function(network, drug_gene_matrix, target_source){
#   network_full <- network %>%
#     rename(gene1= gene2, gene2= gene1) %>%
#     bind_rows(network)
#   target_gene <- target_source$target
#   drug_list <- sort(unique(drug_gene_matrix$drug))
#   gene_list <-  colnames(drug_gene_matrix)[-1]
#   target_list <- NULL
#   i=1
#   res_list <- as.list(1:length(drug_list))
#   # lets use a list column to assign the result and withdraw it later
#   for (drug_i in drug_list){
#     target_list <- NULL
#     # tic()
#     target_i <-  target_source %>% filter(drug %in%  drug_i) %>% pull(target)
#     first_n <- setdiff(network_full %>% filter( gene1 %in% target_i) %>% pull(gene2) , target_i)
#     first_n <- c(first_n,"NA")
#     second_n <- setdiff(network_full %>% filter( gene1 %in% first_n) %>% pull(gene2) , union(target_i,first_n))
#     second_n <- c(second_n,"NA")
#     third_n <- setdiff(
#       x = network_full %>% filter(gene1 %in% second_n) %>% pull(gene2),
#       y = union(x = union(target_i, first_n),y = second_n)
#     )
#     third_n <- c(third_n,"NA")
#     for (gene in  gene_list ) { if
#       (gene %in%  target_i)   {target_cat <- "0"} else if
#       (gene %in%  first_n)  {target_cat <- "1"} else if
#       (gene %in%  second_n) {target_cat <- "2"} else if
#       (gene %in%  third_n)  {target_cat <- "3"} else {target_cat <- "4"}
#       target_list <- c(target_list,target_cat)}
#     res_list[i] <- list(target_list)
#     print(i)
#     # check how to change this print out to a progress bar and think how does input size
#     # impacting this
#     # https://github.com/r-lib/progress
#     i=i+1
#   }
  
  
# lets write a new function that can annotate a drug gene pair according to the given network and target source
  
# get_network_anno_for_drug_gene_pair_full <- function(drug_i, gene_i, network, target_source){
#   # network is a table with two columns named gene1 and gene2
#   # target_source is a table with two columns named drug and gene
# 
#   network_gene_list <- unique(network$gene1, network$gene2)
#   target_list <- target_source %>% filter(drug==drug_i) %>% pull(target)
# 
#   if (length(target_list)==0) {return("NA_target")  }    #drug target not avalible
# 
#   if (gene %in% target_list) {return("target")}
# 
#   # note a drug can have multiple NAs!
#   if (( length(intersect(network_gene_list, target_gene_list)) ==0 ) | !(gene_i %in% network_gene_list) )
#     {return("NA_network") } else {
#       target_list <- intersect(target_list, network_gene_list)
#       network <- graph_from_data_frame(network,directed = F)
#       node_list <- names(V(network))
#       tmp <- distances(graph, v = match(gene_i, node_list), to = match(target_list, node_list))
#       tmp1 <- min(tmp)-1
#       return(tmp1)
#     }
#   }
# 

# can this analysis be generalized ?
# perhaps not only to drug-target-network
#but also other similiar structure where neighbour and network degree information is important



  
# annotate_feature_matrix <- function(network, drug_gene_matrix, target_source) {
#   graph_network <- graph_from_data_frame(network,directed = F)
#   graph_node_list <- names(V(graph_network))
#   n= nrow(drug_gene_matrix)
#   m= ncol(drug_gene_matrix)-1
#   drug_list <- drug_gene_matrix$drug
#   gene_list <- colnames(drug_gene_matrix)[-1]
#   target_list <- unique(target_source$target)
#   drug_with_target <- unique(target_source$drug)
# 
#   NA_drugs_idx <- (1:n)[!(drug_list %in%  drug_with_target)]# index drug without target
#   NA_genes_idx <- (1:m)[!(gene_list %in% c(target_list, graph_node_list))] # index genes that are neither in targetlist or in networknodelist
#   remaining_drugs_idx <- setdiff(1:n, NA_drugs_idx)
#   remaining_genes_idx <- setdiff(1:m, NA_genes_idx)
# 
# 
#   #use character matrix to store result.
#   res <- matrix(nrow = n, ncol = m)
#   colnames(res) <- gene_list
#   row.names(res) <- drug_list
# 
# 
#   if (length(NA_drugs_idx)!=0){  res[NA_drugs_idx,] <- "NA_drug"}
#   if (length(NA_genes_idx)!=0){  res[,NA_genes_idx] <- "NA_gene"}
#   print(length(remaining_genes_idx))
#   print(length(remaining_drugs_idx))
#   print(remaining_drugs_idx)
# 
#   for (drug_i in remaining_drugs_idx) {
#     print(drug_i)
#     drug_name= drug_list[drug_i]
#     target_specific <- target_source %>% filter(drug== drug_name) %>% pull(target)
#     for (gene_i in remaining_genes_idx) {
#       gene_name <- gene_list[gene_i]
# 
# 
#       if (gene_name %in% target_specific) {
#         res[drug_i,gene_i] <- "Target"
#       }
#       else if (!(gene_name %in% graph_node_list)) {
#         res[drug_i,gene_i] <- "NA_gene"}
#       else if(length(intersect(target_specific, graph_node_list)) == 0) {
#           res[drug_i,gene_i] <- "NA_target_network"}
#       else {
#         tmp <- distances(graph_network,
#                          v = match(gene_name, graph_node_list),
#                          to = na.omit(match(target_specific, graph_node_list ))
#         )
#         tmp1 <- min(tmp)-1
#         if (is.infinite(tmp1)){tmp1 <- "Not Connected"}
#         res[drug_i,gene_i] <- as.character(tmp1)
#       }
#     }
#   }
#   res <-   res %>%
#   as_tibble(rownames = NA) %>%
#   rownames_to_column(var = "drug")
#   return(res)
# }


# annotate_feature_matrix_par <- function(network, drug_gene_matrix, target_source) {
#   graph_network <- graph_from_data_frame(network,directed = F)
#   graph_node_list <- names(V(graph_network))
#   n= nrow(drug_gene_matrix)
#   m= ncol(drug_gene_matrix)-1
#   drug_list <- drug_gene_matrix$drug
#   gene_list <- colnames(drug_gene_matrix)[-1]
#   target_list <- unique(target_source$target)
#   drug_with_target <- unique(target_source$drug)
#   
#   NA_drugs_idx <- (1:n)[!(drug_list %in%  drug_with_target)]# index drug without target
#   NA_genes_idx <- (1:m)[!(gene_list %in% c(target_list, graph_node_list))] # index genes that are neither in targetlist or in networknodelist
#   remaining_drugs_idx <- setdiff(1:n, NA_drugs_idx)
#   remaining_genes_idx <- setdiff(1:m, NA_genes_idx)
#   
#   
#   #use character matrix to store result.
#   res <- matrix(nrow = n, ncol = m)
#   colnames(res) <- gene_list
#   row.names(res) <- drug_list
#   
#   
#   if (length(NA_drugs_idx)!=0){  res[NA_drugs_idx,] <- "NA_drug"}
#   if (length(NA_genes_idx)!=0){  res[,NA_genes_idx] <- "NA_gene"}
#   print(length(remaining_genes_idx))
#   # print(length(remaining_drugs_idx))
#   # print(remaining_drugs_idx)
#   
#   w <-   foreach (remaining_drugs_idx) %dopar% {
#     drug_name= drug_list[remaining_drugs_idx]
#     target_specific <- target_source %>% filter(drug== drug_name) %>% pull(target)
#     res0 <- NULL
#     for (gene_i in remaining_genes_idx) {
#       # print(gene_i)
#       gene_name <- gene_list[gene_i]
#       
#       if (gene_name %in% target_specific) { 
#         res0 <- c(res0,"Target")} 
#       else if (!(gene_name %in% graph_node_list)) {
#         res0 <- c(res0,"NA_gene")}  
#       else if(length(intersect(target_specific, graph_node_list)) == 0) {
#         res0 <- c(res0,"NA_target_network")} 
#       else {
#         tmp <- distances(graph_network, 
#                          v = match(gene_name, graph_node_list), 
#                          to = na.omit(match(target_specific, graph_node_list ))
#         )
#         tmp1 <- min(tmp)-1
#         if (is.infinite(tmp1)){tmp1 <- "Not Connected"}
#         tmp1 <- as.character(tmp1)
#         res0 <- c(res0,tmp1)
#       }
#       # print(res0)
#     }
#     res0
#   }
#   # print(res)
# 
#   for(i in  1:length(remaining_drugs_idx)){
#     idx <- remaining_drugs_idx[i]
#     res[idx, remaining_genes_idx] <- w[[i]]
#     print(idx)
#     print(w[[i]])
#   }
#   
#   res <-   res %>%
#     as_tibble(rownames = NA) %>%
#     rownames_to_column(var = "drug")
#   return(res)
# }


annotate_feature_matrix_par <- function(network, drug_gene_matrix, target_source) {
  graph_network <- graph_from_data_frame(network,directed = F)
  graph_node_list <- names(V(graph_network))
  n= nrow(drug_gene_matrix)
  m= ncol(drug_gene_matrix)-1
  drug_list <- drug_gene_matrix$drug
  gene_list <- colnames(drug_gene_matrix)[-1]
  target_list <- unique(target_source$target)
  drug_with_target <- unique(target_source$drug)
  
  NA_drugs_idx <- (1:n)[!(drug_list %in%  drug_with_target)]# index drug without target
  NA_genes_idx <- (1:m)[!(gene_list %in% c(target_list, graph_node_list))] # index genes that are neither in targetlist or in networknodelist
  remaining_drugs_idx <- setdiff(1:n, NA_drugs_idx)
  remaining_genes_idx <- setdiff(1:m, NA_genes_idx)
  
  
  #use character matrix to store result.
  res <- matrix(nrow = n, ncol = m)
  colnames(res) <- gene_list
  row.names(res) <- drug_list
  
  
  if (length(NA_drugs_idx)!=0){  res[NA_drugs_idx,] <- "NA_drug"}
  if (length(NA_genes_idx)!=0){  res[,NA_genes_idx] <- "NA_gene"}
  print(length(remaining_genes_idx))
  # print(length(remaining_drugs_idx))
  # print(remaining_drugs_idx)
  
  w <-   future_map(
    .x= remaining_drugs_idx, 
    .f= function(remaining_drugs_idx){
              
      drug_name= drug_list[remaining_drugs_idx]
      target_specific <- target_source %>% filter(drug== drug_name) %>% pull(target)
      res0 <- NULL
      for (gene_i in remaining_genes_idx) {
        # print(gene_i)
        gene_name <- gene_list[gene_i]
        
        if (gene_name %in% target_specific) { 
          res0 <- c(res0,"Target")} 
        else if (!(gene_name %in% graph_node_list)) {
          res0 <- c(res0,"NA_gene")}  
        else if(length(intersect(target_specific, graph_node_list)) == 0) {
          res0 <- c(res0,"NA_target_network")} 
        else {
          tmp <- distances(graph_network, 
                           v = match(gene_name, graph_node_list), 
                           to = na.omit(match(target_specific, graph_node_list ))
          )
          tmp1 <- min(tmp)
          if (is.infinite(tmp1)){tmp1 <- "Not Connected"}
          tmp1 <- as.character(tmp1)
          res0 <- c(res0,tmp1)
        }
        # print(res0)
      }
      return(res0)
      } 
    )
  
  for(i in  1:length(remaining_drugs_idx)){
    idx <- remaining_drugs_idx[i]
    res[idx, remaining_genes_idx] <- w[[i]]
    # print(idx)
    # print(w[[i]])
  }
  
  res <-   res %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "drug")
  
  return(res)
}


# tmp <- annotate_feature_matrix(network = kegg_tibble,
#                         drug_gene_matrix = feature_imp_ridge_ctrp_comb1,
#                         target_source = ctrp_target_binary)

    



# use the name attributes    
# At the implementation level, a vertex sequence is simply a vector containing numeric vertex ids, but it has a special class attribute which makes it possible to perform graph specific operations on it, like selecting a subset of the vertices based on graph structure, or vertex attributes.

    
      
  
#   
#   target_category <- expand.grid( gene_list,drug_list)
#   target_category <- cbind(target_category,unlist(res_list))
#   colnames(target_category) <- c( "gene","drug", "target_cat")
#   target_category <- as.data.frame(target_category) %>% as_tibble()
#   
#   res <- drug_gene_matrix%>% 
#     pivot_longer(cols= -drug, names_to="gene",values_to="imp") %>% 
#     inner_join(target_category,by=c("drug", "gene"))
#   
#   return(res)
# }

#for getting igraph representation of GO, check:http://igraph.wikidot.com/r-recipes#toc0
# GO is just gene collection without interaction pairs.


# 1 for Sensitivity based signatures:

# ## kegg
# CTRP_binary_kegg_Consen_sig <-
#   annotate_feature_matrix(
#     network = kegg_tibble,
#     drug_gene_matrix = feature_imp_ridge_ctrp_comb1,
#     target_source=  ctrp_target_binary
#   ) 
# 
# write_csv(CTRP_binary_kegg_Consen_sig,"CTRP_binary_kegg_Consen_sig.csv")
# 
# 
# GDSC_binary_kegg_Consen_sig <-
#   annotate_feature_matrix(
#     network = kegg_tibble,
#     drug_gene_matrix = feature_imp_ridge_gdsc_comb1,
#     target_source=  gdsc_target_binary
#   ) 
# 
# write_csv(GDSC_binary_kegg_Consen_sig,"GDSC_binary_kegg_Consen_sig.csv")
# 
# PRISM_binary_kegg_Consen_sig <-
#   annotate_feature_matrix(
#     network = kegg_tibble,
#     drug_gene_matrix = feature_imp_ridge_prism_comb1,
#     target_source=  prism_target_binary
#   ) 
# 
# write_csv(PRISM_binary_kegg_Consen_sig,"PRISM_binary_kegg_Consen_sig.csv")

##PPI old
# tic()
# CTRP_binary_PPIold_Consen_sig <-
#   annotate_feature_matrix(
#     network = PPI_old,
#     drug_gene_matrix = feature_imp_ridge_ctrp_comb1[1:8,],
#     target_source=  ctrp_target_binary
#   ) 
# toc()



# tic()
# plan(sequential) # 
# CTRP_binary_PPIold_Consen_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = feature_imp_ridge_ctrp_comb1[1:30,1:1000],
#     target_source=  ctrp_target_binary
#   )
# toc()
# 
# tic()
# plan(multicore, workers =4) #
# CTRP_binary_PPIold_Consen_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = feature_imp_ridge_ctrp_comb1[1:30,1:1000],
#     target_source=  ctrp_target_binary
#   ) 
# toc()
# 
# tic()
# plan(multicore, workers =8) #fastest 
# CTRP_binary_PPIold_Consen_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = feature_imp_ridge_ctrp_comb1[1:30,1:1000],
#     target_source=  ctrp_target_binary
#   ) 
# toc()
# 
# tic()
# plan(multisession, workers =4) 
# CTRP_binary_PPIold_Consen_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = feature_imp_ridge_ctrp_comb1[1:30,1:1000],
#     target_source=  ctrp_target_binary
#   ) 
# toc()
# 
# 
# tic()
# plan(multisession, workers =8) 
# CTRP_binary_PPIold_Consen_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = feature_imp_ridge_ctrp_comb1[1:30,1:1000],
#     target_source=  ctrp_target_binary
#   ) 
# toc()

# plan(multicore, workers =40)
# tic()
# CTRP_binary_PPIold_Consen_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = feature_imp_ridge_ctrp_comb1,
#     target_source=  ctrp_target_binary
#   )
# toc()
# write_csv(CTRP_binary_PPIold_Consen_sig,"CTRP_binary_PPIold_Consen_sig.csv")
# 
# 
# tic()
# GDSC_binary_PPIold_Consen_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = feature_imp_ridge_gdsc_comb1,
#     target_source=  gdsc_target_binary
#   )
# toc()
# write_csv(GDSC_binary_PPIold_Consen_sig,"GDSC_binary_PPIold_Consen_sig.csv")
# 
# tic()
# PRISM_binary_PPIold_Consen_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = feature_imp_ridge_prism_comb1,
#     target_source=  prism_target_binary
#   )
# toc()
# write_csv(PRISM_binary_PPIold_Consen_sig,"PRISM_binary_PPIold_Consen_sig.csv")
# 
# plan(sequential)

# 
# ##PPI new
plan(multicore, workers=40)
tic()
CTRP_binary_PPInew_Consen_sig <-
  annotate_feature_matrix_par(
    network = PPI_new,
    drug_gene_matrix = feature_imp_ridge_ctrp_comb1,
    target_source=  ctrp_target_binary
  )

write_csv(CTRP_binary_PPInew_Consen_sig,"CTRP_binary_PPInew_Consen_sig.csv")


GDSC_binary_PPInew_Consen_sig <-
  annotate_feature_matrix_par(
    network = PPI_new,
    drug_gene_matrix = feature_imp_ridge_gdsc_comb1,
    target_source=  gdsc_target_binary
  )

write_csv(GDSC_binary_PPInew_Consen_sig,"GDSC_binary_PPInew_Consen_sig.csv")


PRISM_binary_PPInew_Consen_sig <-
  annotate_feature_matrix_par(
    network = PPI_new,
    drug_gene_matrix = feature_imp_ridge_prism_comb1,
    target_source=  prism_target_binary
  )

write_csv(PRISM_binary_PPInew_Consen_sig,"PRISM_binary_PPInew_Consen_sig.csv")


plan(sequential, workers=40)
toc()

# 2for ConSen sig

# ## KEGG
# CTRP_binary_KEGG_Conexp_sig <-
#   annotate_feature_matrix_par(
#     network = kegg_tibble,
#     drug_gene_matrix = drug_consensus_ctrp,
#     target_source=  ctrp_target_binary
#   )
# write_csv(CTRP_binary_KEGG_Conexp_sig,"CTRP_binary_KEGG_Conexp_sig.csv")
# 
# GDSC_binary_KEGG_Conexp_sig <-
#   annotate_feature_matrix_par(
#     network = kegg_tibble,
#     drug_gene_matrix = drug_consensus_gdsc,
#     target_source=  gdsc_target_binary
#   )
# write_csv(GDSC_binary_KEGG_Conexp_sig,"GDSC_binary_KEGG_Conexp_sig.csv")
# 
# PRISM_binary_KEGG_Conexp_sig <-
#   annotate_feature_matrix_par(
#     network = kegg_tibble,
#     drug_gene_matrix = drug_consensus_prism,
#     target_source=  prism_target_binary
#   )
# write_csv(PRISM_binary_KEGG_Conexp_sig,"PRISM_binary_KEGG_Conexp_sig.csv")

# tic()
# ## PPI_old
# CTRP_binary_PPIold_Conexp_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = drug_consensus_ctrp,
#     target_source=  ctrp_target_binary
#   )
# write_csv(CTRP_binary_PPIold_Conexp_sig,"CTRP_binary_PPIold_Conexp_sig.csv")
# 
# GDSC_binary_PPIold_Conexp_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = drug_consensus_gdsc,
#     target_source=  gdsc_target_binary
#   )
# write_csv(GDSC_binary_PPIold_Conexp_sig,"GDSC_binary_PPIold_Conexp_sig.csv")
# 
# PRISM_binary_PPIold_Conexp_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_old,
#     drug_gene_matrix = drug_consensus_prism,
#     target_source=  prism_target_binary
#   )
# write_csv(PRISM_binary_PPIold_Conexp_sig,"PRISM_binary_PPIold_Conexp_sig.csv")
# 
# 
# ## PPI_new
# CTRP_binary_PPInew_Conexp_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_new,
#     drug_gene_matrix = drug_consensus_ctrp,
#     target_source=  ctrp_target_binary
#   )
# write_csv(CTRP_binary_PPInew_Conexp_sig,"CTRP_binary_PPInew_Conexp_sig.csv")
# 
# GDSC_binary_PPInew_Conexp_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_new,
#     drug_gene_matrix = drug_consensus_gdsc,
#     target_source=  gdsc_target_binary
#   )
# write_csv(GDSC_binary_PPInew_Conexp_sig,"GDSC_binary_PPInew_Conexp_sig.csv")
# 
# PRISM_binary_PPInew_Conexp_sig <-
#   annotate_feature_matrix_par(
#     network = PPI_new,
#     drug_gene_matrix = drug_consensus_prism,
#     target_source=  prism_target_binary
#   )
# write_csv(PRISM_binary_PPInew_Conexp_sig,"PRISM_binary_PPInew_Conexp_sig.csv")
# toc()