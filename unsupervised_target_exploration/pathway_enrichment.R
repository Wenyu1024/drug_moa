# show the drug perturbation at Pathway level:
# gene-level features provided a way to establish a pathway perturbation score,
# 
# Where we can map a bulk-level drug perturbation effect to perturbation of specific pathways.

# pathways from
# http://www.gsea-msigdb.org/gsea/downloads.jsp  
# this analysis is not so heavy. It can be finished with just one node in several minutes

library(tidyverse)
library(fgsea)
library(furrr)

# feature_imp_ridge_ctrp_ces1 <- read_csv("~/cluster_scratch/glmnet_modelling/result_largevm/feature_imp_ridge_ctrp_ces1.csv")
load("/scratch/project_2003466/forward_modelling/targetpred_output_simplefiltering.RData")
setwd("/scratch/project_2003466/prior/")
pathway_go <- gmtPathways("c5.go.v7.4.symbols.gmt")
pathway_kegg <- gmtPathways("c2.cp.kegg.v7.4.symbols.gmt")

return_pathway_for_genefeaturematrix_par <- function(data,pathways){
  # res_pathway <- vector(length = nrow(data), mode = "list")
  res_pathway <- future_map(
    .options = furrr_options(seed = 1),
    .x = 1:(nrow(data)),
    .f = function(x){
      rank <- unlist(data[x,-1])
      fgseaRes <- fgsea(pathways, rank, minSize=15, maxSize=500)
      return(fgseaRes)})
  return(res_pathway)
  }

# GO CTRP
plan(multisession, workers=40)

ctrp_kegg <- return_pathway_for_genefeaturematrix_par(
  data= feature_imp_ridge_ctrp_comb1[,-1],pathways = pathway_kegg)
print("ctrp_kegg finished")
ctrp_go <- return_pathway_for_genefeaturematrix_par(
  data= feature_imp_ridge_ctrp_comb1[,-1],pathways = pathway_go)
print("ctrp_go finished")

gdsc_kegg <- return_pathway_for_genefeaturematrix_par(
  data= feature_imp_ridge_gdsc_comb1[,-1],pathways = pathway_kegg)
print("gdsc_kegg finished")
gdsc_go <- return_pathway_for_genefeaturematrix_par(
  data= feature_imp_ridge_gdsc_comb1[,-1],pathways = pathway_go)
print("gdsc_go finished")

prism_kegg <- return_pathway_for_genefeaturematrix_par(
  data= feature_imp_ridge_prism_comb1[,-1],pathways = pathway_kegg)
print("prism_kegg finished")
prism_go <- return_pathway_for_genefeaturematrix_par(
  data= feature_imp_ridge_prism_comb1[,-1],pathways = pathway_go)
print("prism_go finished")
plan(sequential)

# save.image("/scratch/project_2003466/glmnet_modelling/target_pred_server/res_pathway_all.RData")
setwd("/scratch/project_2003466/glmnet_modelling/target_pred_server")
save(list = c("ctrp_kegg","ctrp_go", "gdsc_kegg", "gdsc_go", "prism_kegg", "prism_go" ),
     file =  "res_pathways_all.RData")

# return_pathway_for_genefeaturematrix_par <- function(data,pathways){
#   # res_pathway <- vector(length = nrow(data), mode = "list")
#   res_pathway <- future_map(
#     .options = furrr_options(seed = 1),
#     .x = 1:(nrow(data)),
#     .f = function(x){  
#       rank <- abs(unlist(data[x,]))  # use the absolute value for the rank.
#       fgseaRes <- fgsea(pathways, rank, minSize=15, maxSize=500)
#       return(fgseaRes)})
#   return(res_pathway)
# }
# 
# # GO CTRP
# plan(multisession, workers=40)
# 
# ctrp_kegg <- return_pathway_for_genefeaturematrix_par(
#   data= feature_imp_ridge_ctrp_comb1[,-1],pathways = pathway_kegg)
# print("ctrp_kegg finished")
# ctrp_go <- return_pathway_for_genefeaturematrix_par(
#   data= feature_imp_ridge_ctrp_comb1[,-1],pathways = pathway_go)
# print("ctrp_go finished")
# 
# gdsc_kegg <- return_pathway_for_genefeaturematrix_par(
#   data= feature_imp_ridge_gdsc_comb1[,-1],pathways = pathway_kegg)
# print("gdsc_kegg finished")
# gdsc_go <- return_pathway_for_genefeaturematrix_par(
#   data= feature_imp_ridge_gdsc_comb1[,-1],pathways = pathway_go)
# print("gdsc_go finished")
# 
# prism_kegg <- return_pathway_for_genefeaturematrix_par(
#   data= feature_imp_ridge_prism_comb1[,-1],pathways = pathway_kegg)
# print("prism_kegg finished")
# prism_go <- return_pathway_for_genefeaturematrix_par(
#   data= feature_imp_ridge_prism_comb1[,-1],pathways = pathway_go)
# print("prism_go finished")
# plan(sequential)
# 
# 
# setwd("/scratch/project_2003466/glmnet_modelling/target_pred_server")
# save(list = c("ctrp_kegg","ctrp_go", "gdsc_kegg", "gdsc_go", "prism_kegg", "prism_go" ), 
#      file =  "res_pathways_all_abs.RData")
# 
