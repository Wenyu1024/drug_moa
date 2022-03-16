load("~/cluster_scratch/forward_modelling/targetpred_output_simplefiltering.RData")
# Navitoclax 
match("BRD-K82746043", colnames(ctrpess_drugbank_predictions))
supervised_drugbankpredictions_Navitoclax_ctrp <- ctrpess_drugbank_predictions[,423]

match("BRD-K82746043", colnames(ctrpess_dtc_predictions))
supervised_dtcpredictions_Navitoclax_ctrp <- ctrpess_dtc_predictions[,423]

unsupervised_predictions_Navitoclax_ctrp <- feature_imp_ridge_ctrp_comb1 %>% filter(drug== "BRD-K82746043") %>% 
  pivot_longer(-drug,names_to = "gene", values_to="value") %>% 
  mutate(value_rank= rank(-value,ties.method = "min"))


match("1011", colnames(gdscess_drugbank_predictions))
supervised_drugbankpredictions_Navitoclax_gdsc <- gdscess_drugbank_predictions[,7]
match("1011", colnames(gdscess_dtc_predictions))
supervised_dtcpredictions_Navitoclax_gdsc <- gdscess_dtc_predictions[,7]

unsupervised_predictions_Navitoclax_gdsc = feature_imp_ridge_gdsc_comb1 %>% filter(drug== 1011) %>% 
  pivot_longer(-drug,names_to = "gene", values_to="value") %>% 
  mutate(value_rank= rank(-value,ties.method = "min"))

save(list = c("supervised_drugbankpredictions_Navitoclax_ctrp", "supervised_dtcpredictions_Navitoclax_ctrp",
              "unsupervised_predictions_Navitoclax_ctrp","unsupervised_predictions_Navitoclax_gdsc",
              "supervised_dtcpredictions_Navitoclax_gdsc","supervised_drugbankpredictions_Navitoclax_gdsc"),
     file= "~/cluster_wrk/drug_moa/supervised_target_pred/Navitoclax_Alina.RData")


#Bortezomib BRD-K88510285
match("BRD-K88510285", colnames(ctrpess_drugbank_predictions))
supervised_drugbankpredictions_Bortezomib_ctrp <- ctrpess_drugbank_predictions[,446]
match("BRD-K88510285", colnames(ctrpess_dtc_predictions))
supervised_dtcpredictions_Bortezomib_ctrp <- ctrpess_dtc_predictions[,446]

unsupervised_predictions_Bortezomib_ctrp <- feature_imp_ridge_ctrp_comb1 %>% filter(drug== "BRD-K88510285") %>% 
  pivot_longer(-drug,names_to = "gene", values_to="value") %>% 
  mutate(value_rank= rank(-value,ties.method = "min"))


match("1191", colnames(gdscess_drugbank_predictions))
supervised_drugbankpredictions_Bortezomib_gdsc <- gdscess_drugbank_predictions[,38]
match("1191", colnames(gdscess_dtc_predictions))
supervised_dtcpredictions_Bortezomib_gdsc <- gdscess_dtc_predictions[,38]

unsupervised_predictions_Bortezomib_gdsc = feature_imp_ridge_gdsc_comb1 %>% filter(drug== 1191) %>% 
  pivot_longer(-drug,names_to = "gene", values_to="value") %>% 
  mutate(value_rank= rank(-value,ties.method = "min"))

save(list = c("supervised_drugbankpredictions_Bortezomib_ctrp", "supervised_dtcpredictions_Bortezomib_ctrp",
              "unsupervised_predictions_Bortezomib_ctrp","unsupervised_predictions_Bortezomib_gdsc",
              "supervised_dtcpredictions_Bortezomib_gdsc","supervised_drugbankpredictions_Bortezomib_gdsc"),
     file= "~/cluster_wrk/drug_moa/supervised_target_pred/Bortezomib_Alina.RData")
