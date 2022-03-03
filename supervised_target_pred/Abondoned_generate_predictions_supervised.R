# fix a training set and then use this trainset to predict target
# for all the avalible drugs

# for acc evaluation we use hold out set. now we are applying fitted models,
# now this is not needed as we are not evaluating accuracy but instead applying the models.

# the training drugs are then all the drugs with good sen-sig and drug targe
# the drugs to be predicted are all the drugs with good sen-sig.

#question, how to choose the right drug?

train_consen_sig <- bind_rows(
  feature_imp_ridge_prism_comb1 %>%
    filter(drug %in% prism_good_fitted_drugs_annotated) %>%
    arrange(drug),
  
  feature_imp_ridge_prism_comb1 %>%
    filter(drug %in% prism_good_fitted_nononcologydrugs) %>%
    arrange(drug)
)


train_cor_mat <- cor(t(train_consen_sig %>% select(-drug)))
colnames(train_cor_mat) <- train_consen_sig$drug
row.names(train_cor_mat) <- train_consen_sig$drug
target_mat_268 <- get_target_mat(
  target_tibble = prism_target_binary %>%
    filter(drug %in% prism_good_fitted_drugs_annotated) %>%
    arrange(drug))


pred_268_supervised <- map_dfc(
  .x = 1:268,
  .f = ~no_tunning_weighted_averaging(
    target_mat = target_mat_268,
    cor_mat = train_cor_mat[1:268,1:268],
    test_idx = .x,
    pred_all = T)
) 


colnames(pred_268_supervised) <- row.names(target_mat_268)
cor_df_predicted_target <- 
  cor(pred_268_supervised, use = "pairwise.complete.obs", method = "spearman") %>%
  as_tibble() %>% 
  mutate(drug1=row.names(target_mat_268)) %>%
  pivot_longer(cols = -drug1,names_to= "drug2") %>% 
  arrange(drug1, drug2)

cor_df_unpred_target <- 
  train_cor_mat[1:268,1:268] %>% 
  as_tibble() %>% 
  mutate(drug1=row.names(target_mat_268)) %>%
  pivot_longer(cols = -drug1,names_to= "drug2") %>% 
  arrange(drug1, drug2)


cor_df_sen <- prism_data %>% filter(BROAD_ID %in% row.names(target_mat_268)) %>% 
  select(BROAD_ID, sensitivity) %>% 
  unnest(sensitivity) %>% 
  select(-ec50) %>% 
  drop_na() %>% 
  pivot_wider(id_cols = depmap_id,names_from= BROAD_ID, values_from=auc) %>% 
  select(-depmap_id) %>% 
  cor(use = "pairwise.complete.obs",method = "spearman") 

colnames(cor_df_sen) <- row.names(cor_df_sen)
cor_df_sen <- as_tibble(cor_df_sen)  %>% 
  # rowid_to_column()
  mutate(drug1= row.names(cor_df_sen)) %>% 
  pivot_longer(cols = -drug1,names_to= "drug2") %>% arrange(drug1,drug2)
