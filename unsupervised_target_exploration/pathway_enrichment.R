# show the drug perturbation at Pathway level:
# CES1 based features provided a way to establish a pathway perturbation score,
# 
# Where we can map a bulk-level drug perturbation effect to perturbation of specific pathways.

# pathways from
# http://www.gsea-msigdb.org/gsea/downloads.jsp  

library(fgsea)
feature_imp_ridge_ctrp_ces1 <- read_csv("~/cluster_scratch/glmnet_modelling/result_largevm/feature_imp_ridge_ctrp_ces1.csv")

pathways <- gmtPathways("/home/cloud-user/cluster_scratch/prior/c5.go.v7.4.symbols.gmt")
res_pathway <- vector(length = 365, mode = "list")
for (i in 1:365){
  rank <- unlist(feature_imp_ridge_ctrp_ces1[i,-1])
  fgseaRes <- fgsea(pathways, rank, minSize=15, maxSize=500)
  res_pathway[[i]] <- fgseaRes
  print(i)
}
save(res_pathway,file =  "~/cluster_scratch/glmnet_modelling/target_pred_server/res_pathway_go.RData")


pathways <- gmtPathways("/home/cloud-user/cluster_scratch/prior/c2.cp.kegg.v7.4.symbols.gmt")
res_pathway <- vector(length = 365, mode = "list")
for (i in 1:365){
  rank <- unlist(feature_imp_ridge_ctrp_ces1[i,-1])
  fgseaRes <- fgsea(pathways, rank, minSize=15, maxSize=500)
  res_pathway[[i]] <- fgseaRes
  print(i)
}
save(res_pathway,file =  "~/cluster_scratch/glmnet_modelling/target_pred_server/res_pathway_kegg.RData")