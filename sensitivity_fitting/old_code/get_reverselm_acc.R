# reveres modeling with 465 ess and 365 ctrpv2 sensi

# lets forget about previous code and try to revise based on the current framework created for glmnet.






options(expressions = 5e5)
get_lm <- function(training_data){
  lm_model <-
    linear_reg() %>%
    set_mode("regression") %>%
    set_engine("lm")
  
  lm_recipe <- recipe( training_data  ) %>%
    update_role(everything()) %>%
    update_role(area_under_curve, new_role = "outcome") %>% 
    step_medianimpute(all_numeric()) %>% 
    step_normalize(all_predictors())
  
  
  lm_wflow <- workflow( ) %>% 
    add_model(lm_model) %>% 
    add_recipe(lm_recipe)
  
  lm_wflow_final_fit <- fit(lm_wflow, data = training_data)
  return(lm_wflow_final_fit)
}

load("~/drug_target/data/clean_ctrpv2_data.RData")

ceres_465_imputed <- read_csv("~/cluster_wrk/impute_ess/imputed/ceres_465_imputed.csv")
demeter2_465_imputed <- read_csv("~/cluster_wrk/impute_ess/imputed/demeter2_465_imputed.csv")
ces1_465_imputed <- read_csv("~/cluster_wrk/impute_ess/imputed/ces1_465_imputed.csv")
ces2_465_imputed <- read_csv("~/cluster_wrk/impute_ess/imputed/ces2_465_imputed.csv")
ceres_full_imputed <- read_csv("~/cluster_wrk/impute_ess/imputed/ceres_full_imputed.csv")
demeter2_full_imputed <- read_csv("~/cluster_wrk/impute_ess/imputed/demeter2_full_imputed.csv")
demeter2_filtered_imputed <- read_csv("~/cluster_wrk/impute_ess/imputed/demeter2_filtered_imputed.csv")


# ess <- demeter2_465_imputed
# target <- ctrpv2_data$target[[1]]$gene_symbol_of_protein_target

estimate_performance_reverselm <- function(sen_df, ess, target){
  target <- target$gene_symbol_of_protein_target
  idx <- sum(target %in% colnames(ess)) < length(target)
  if (idx){res= NA} else{
    df <- sen_df  %>%   as_tibble() %>% 
      inner_join( 
        y= janitor::clean_names(
          ess %>% 
            select(one_of(c("DepMap_ID", target)  )) %>% 
            rename_at(vars(!contains("DepMap_ID")), .fun= ~paste0(., "_predictors"))
        ), 
        by= c("DepMap_ID"="dep_map_id")) %>% 
      select(-master_ccl_id, -apparent_ec50_umol, -DepMap_ID) 
    
    # set.seed(0000)
    df1 <-  df %>% 
      vfold_cv(v = 5) %>% 
      mutate(train_data= map(.x = splits, .f = training)) %>% 
      mutate(test_data= map(.x = splits, .f = ~testing(.x))) %>% 
      mutate(lm_wflow_final_fit = 
               map(.x = train_data,
                   .f = ~get_lm(training_data = .x))) %>%  
      mutate(pred_test = map2(.x = lm_wflow_final_fit,
                              .y = test_data,
                              .f = ~predict(.x, new_data=.y))) %>% 
      mutate(spearman_cor =  map2_dbl(.x= pred_test, .y= test_data,
                                      .f= function(x,y){
                                        cor(x = x %>% pull(.pred), y = y %>% pull(area_under_curve), method = "spearman")}))
    res= df1 %>% pull(spearman_cor)
    
  }
  
  
  return(res)
}


# save.image("~/drug_target/analysis/reverse_modeling_data.RData")
library(furrr)
set.seed(0000)
future::plan(multiprocess)
data <- ctrpv2_data %>%
  ungroup() 
  # slice(1:10) 

a <- Sys.time()
res2 <- data %>% 
  mutate(ces1_perf= future_map2(sensitivity, target,
                          .f = ~estimate_performance_reverselm(sen_df = .x,target =.y ,ess = ces1_465_imputed))) %>% 
  mutate(ces2_perf= future_map2(sensitivity, target,
                             .f = ~estimate_performance_reverselm(sen_df = .x,target =.y ,ess = ces2_465_imputed))) %>% 
  mutate(ceres_perf= future_map2(sensitivity, target,
                             .f = ~estimate_performance_reverselm(sen_df = .x,target =.y ,ess = ceres_465_imputed))) %>% 
  mutate(demeter2_perf= future_map2(sensitivity, target,
                                   .f = ~estimate_performance_reverselm(sen_df = .x,target =.y ,ess = demeter2_465_imputed))) %>% 
  mutate(ceres_full_perf= future_map2(sensitivity, target,
                           .f = ~estimate_performance_reverselm(sen_df = .x,target =.y ,ess = ceres_full_imputed))) %>% 
  mutate(demeter2_full_perf= future_map2(sensitivity, target,
                           .f = ~estimate_performance_reverselm(sen_df = .x,target =.y ,ess = demeter2_full_imputed))) %>% 
  mutate(demeter2_filtered_perf= future_map2(sensitivity, target,
                           .f = ~estimate_performance_reverselm(sen_df = .x,target =.y ,ess = demeter2_filtered_imputed)))

b <- Sys.time()
b-a



