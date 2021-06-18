# The function takes three input and 
# output a long table of feature importance matrix with target network annotation:
# the function may take several minute to run, depend on how many rows of the input are there.
library(tidyverse)
library(tidymodels)
load("~/cluster_scratch/forward_modelling/targetpred_output.RData")
annotate_feature_matrix <- function(network, feature_imp_mat, target_source){
  network_full <- network %>% 
    rename(gene1= gene2, gene2= gene1) %>% 
    bind_rows(network)
  target_gene <- target_source$target
  drug_list <- sort(unique(feature_imp_mat$drug))
  gene_list <-  colnames(feature_imp_mat)[-1]
  target_list <- NULL
  i=1
  res_list <- as.list(1:length(drug_list))
  # lets use a list column to assign the result and withdraw it later
  for (drug_i in drug_list){
    target_list <- NULL
    # tic()
    target_i <-  target_source %>% filter(drug %in%  drug_i) %>% pull(target)
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


# Call the function with KEGG dataset and 

kegg_tibble <- read_csv("~/cluster_scratch/prior/kegg_tibble.csv")
setwd("~/cluster_scratch/prior/")

ces1_glm_feature_annotated_ctrpbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_ctrp_ces1,
    target_source=  ctrp_target_binary
  )
write_csv(ces1_glm_feature_annotated_ctrpbinary_kegg,"ces1_glm_feature_annotated_ctrpbinary_kegg.csv")

ceres_glm_feature_annotated_ctrpbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_ctrp_ceres,
    target_source=  ctrp_target_binary
  )
write_csv(ceres_glm_feature_annotated_ctrpbinary_kegg,"ceres_glm_feature_annotated_ctrpbinary_kegg.csv")

demeter_glm_feature_annotated_ctrpbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_ctrp_demeter2 ,
    target_source=  ctrp_target_binary
  )
write_csv(demeter_glm_feature_annotated_ctrpbinary_kegg,"demeter_glm_feature_annotated_ctrpbinary_kegg.csv")

exp_glm_feature_annotated_ctrpbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_ctrp_exp ,
    target_source=  ctrp_target_binary
  )
write_csv(exp_glm_feature_annotated_ctrpbinary_kegg,"exp_glm_feature_annotated_ctrpbinary_kegg.csv")

# use dtc binarized tibble insteand
ces1_glm_feature_annotated_ctrpdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_ctrp_ces1,
    target_source=  ctrp_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(ces1_glm_feature_annotated_ctrpdtc_kegg,"ces1_glm_feature_annotated_ctrpdtc_kegg.csv")

ceres_glm_feature_annotated_ctrpdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_ctrp_ceres,
    target_source=  ctrp_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(ceres_glm_feature_annotated_ctrpdtc_kegg,"ceres_glm_feature_annotated_ctrpdtc_kegg.csv")

demeter_glm_feature_annotated_ctrpdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_ctrp_demeter2 ,
    target_source=  ctrp_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(demeter_glm_feature_annotated_ctrpdtc_kegg,"demeter_glm_feature_annotated_ctrpdtc_kegg.csv")

exp_glm_feature_annotated_ctrpdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_ctrp_exp,
    target_source=  ctrp_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(exp_glm_feature_annotated_ctrpdtc_kegg,"exp_glm_feature_annotated_ctrpdtc_kegg.csv")

















# gdsc


ces1_glm_feature_annotated_gdscbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_gdsc_ces1,
    target_source=  gdsc_target_binary
  )
write_csv(ces1_glm_feature_annotated_gdscbinary_kegg,"ces1_glm_feature_annotated_gdscbinary_kegg.csv")

ceres_glm_feature_annotated_gdscbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_gdsc_ceres,
    target_source=  gdsc_target_binary
  )
write_csv(ceres_glm_feature_annotated_gdscbinary_kegg,"ceres_glm_feature_annotated_gdscbinary_kegg.csv")

demeter_glm_feature_annotated_gdscbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_gdsc_demeter2 ,
    target_source=  gdsc_target_binary
  )
write_csv(demeter_glm_feature_annotated_gdscbinary_kegg,"demeter_glm_feature_annotated_gdscbinary_kegg.csv")

exp_glm_feature_annotated_gdscbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_gdsc_exp ,
    target_source=  gdsc_target_binary
  )
write_csv(exp_glm_feature_annotated_gdscbinary_kegg,"exp_glm_feature_annotated_gdscbinary_kegg.csv")

# use dtc binarized tibble insteand
ces1_glm_feature_annotated_gdscdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_gdsc_ces1,
    target_source=  gdsc_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(ces1_glm_feature_annotated_gdscdtc_kegg,"ces1_glm_feature_annotated_gdscdtc_kegg.csv")

ceres_glm_feature_annotated_gdscdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_gdsc_ceres,
    target_source=  gdsc_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(ceres_glm_feature_annotated_gdscdtc_kegg,"ceres_glm_feature_annotated_gdscdtc_kegg.csv")

demeter_glm_feature_annotated_gdscdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_gdsc_demeter2 ,
    target_source=  gdsc_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(demeter_glm_feature_annotated_gdscdtc_kegg,"demeter_glm_feature_annotated_gdscdtc_kegg.csv")

exp_glm_feature_annotated_gdscdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_gdsc_exp,
    target_source=  gdsc_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(exp_glm_feature_annotated_gdscdtc_kegg,"exp_glm_feature_annotated_gdscdtc_kegg.csv")

# prism


ces1_glm_feature_annotated_prismbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_prism_ces1,
    target_source=  prism_target_binary
  )
write_csv(ces1_glm_feature_annotated_prismbinary_kegg,"ces1_glm_feature_annotated_prismbinary_kegg.csv")

ceres_glm_feature_annotated_prismbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_prism_ceres,
    target_source=  prism_target_binary
  )
write_csv(ceres_glm_feature_annotated_prismbinary_kegg,"ceres_glm_feature_annotated_prismbinary_kegg.csv")

demeter_glm_feature_annotated_prismbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_prism_demeter2 ,
    target_source=  prism_target_binary
  )
write_csv(demeter_glm_feature_annotated_prismbinary_kegg,"demeter_glm_feature_annotated_prismbinary_kegg.csv")

exp_glm_feature_annotated_prismbinary_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_prism_exp ,
    target_source=  prism_target_binary
  )
write_csv(exp_glm_feature_annotated_prismbinary_kegg,"exp_glm_feature_annotated_prismbinary_kegg.csv")

# use dtc binarized tibble insteand
ces1_glm_feature_annotated_prismdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_prism_ces1,
    target_source=  prism_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(ces1_glm_feature_annotated_prismdtc_kegg,"ces1_glm_feature_annotated_prismdtc_kegg.csv")

ceres_glm_feature_annotated_prismdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_prism_ceres,
    target_source=  prism_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(ceres_glm_feature_annotated_prismdtc_kegg,"ceres_glm_feature_annotated_prismdtc_kegg.csv")

demeter_glm_feature_annotated_prismdtc_kegg <-
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_prism_demeter2 ,
    target_source=  prism_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(demeter_glm_feature_annotated_prismdtc_kegg,"demeter_glm_feature_annotated_prismdtc_kegg.csv")

exp_glm_feature_annotated_prismdtc_kegg <- 
  annotate_feature_matrix(
    network = kegg_tibble,
    feature_imp_mat = feature_imp_ridge_prism_exp,
    target_source=  prism_target_dtc %>% filter(binding_score>0.4) %>% select(-binding_score)
  )
write_csv(exp_glm_feature_annotated_prismdtc_kegg,"exp_glm_feature_annotated_prismdtc_kegg.csv")
