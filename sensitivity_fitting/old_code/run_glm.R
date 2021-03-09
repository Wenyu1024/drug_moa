#this code combines and update the effort in run_glm.R

library(glmnet)
library(furrr)
library(tidyverse)
library(tidymodels)
options(expressions = 5e5)

glm_model <- linear_reg(penalty = tune(), mixture =tune()) %>% 
  set_engine("glmnet") %>% 
  set_mode("regression")

glm_param <- parameters(penalty(range = c(-4,0), 
                                trans = log10_trans()),
                        mixture(range = c(0,1)))

glm_grid <- grid_latin_hypercube(glm_param,size = 20)



get_fixed_glmnet <- function(training_data,training_vfold){

  glm_recipe <- recipe( training_data  ) %>%
    update_role(everything()) %>%
    update_role(area_under_curve, new_role = "outcome") %>% 
    step_normalize(all_predictors())
  
  glm_wflow <- workflow( ) %>% 
    add_model(glm_model) %>% 
    add_recipe(glm_recipe)
  
  glm_search <- tune_grid(
    object = glm_wflow, 
    grid = glm_grid, 
    resamples = training_vfold,
    control = tune::control_grid(verbose = T, allow_par=F)
  )
  
  glm_param_final <- select_best(glm_search, metric = "rsq")
  glm_wflow_final <- finalize_workflow(glm_wflow, glm_param_final)
  glm_wflow_final_fit <- fit(glm_wflow_final, data = training_data)
  return(glm_wflow_final_fit)
}

return_feature_importance <- function(fit){
  fit %>% 
    tidy() %>% 
    filter(term!="(Intercept)") %>% 
    arrange(term) %>% 
    pull(estimate)
}

estimate_performance_glmnet <- function(sen_df, ess){
  df <- sen_df  %>%   as_tibble() %>% 
    inner_join( 
      y= janitor::clean_names(ess)%>% 
        rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
      by= c("DepMap_ID"="dep_map_id")) %>% 
    select(-master_ccl_id, -apparent_ec50_umol, -DepMap_ID) 
  
  df1 <-  df %>% 
    nested_cv(outside=vfold_cv(v = 5,repeats = 1), inside= vfold_cv(v = 5,repeats = 1)) %>% 
    mutate(train_data= map(.x = splits, .f = training)) %>% 
    mutate(test_data= map(.x = splits, .f = ~testing(.x))) %>% 
    mutate(glm_wflow_final_fit = 
             map2(.x = train_data,.y = inner_resamples,
                  .f = ~get_fixed_glmnet(training_data = .x, training_vfold = .y))) %>%  
    mutate(pred_test = map2(.x = glm_wflow_final_fit,
                            .y = test_data,
                            .f = ~predict(.x, new_data=.y))) %>% 
    mutate(spearman_cor =  map2_dbl(.x= pred_test, .y= test_data,
                                    .f= function(x,y){
                                      cor(x = x %>% pull(.pred), y = y %>% pull(area_under_curve), method = "spearman")}))%>% 
    mutate(imp_vec= map(.x = glm_wflow_final_fit, 
                                     .f = ~return_feature_importance(.x)))
  res <- df1 %>% select(spearman_cor, imp_vec)
  return(res)
}



load("/home/cloud-user/cluster_wrk/data/ctrpv2_data.RData")
ces1<- read_csv("/home/cloud-user/cluster_scratch/impute_ess_20q4/ces1_478.csv")
ces2 <- read_csv("/home/cloud-user/cluster_scratch/impute_ess_20q4/ces2_478.csv")
ceres <- read_csv("/home/cloud-user/cluster_scratch/impute_ess_20q4/ceres_478.csv")
demeter2 <- read_csv("/home/cloud-user/cluster_scratch/impute_ess_20q4/demeter2_478.csv")
data <- ctrpv2_data %>%
  ungroup() %>% 
  slice(266:365)

# colnames(ces1)[1] <- "DepMap_ID"
# ces1_small <- ces1[,1:10]
# sen_df <- ctrpv2_data$sensitivity[[1]]
# training_data <- df1$train_data[[1]]
# training_vfold <- df1$inner_resamples[[1]]
# tmp <- estimate_performance_glmnet(sen_df = sen_df, ess = ces1_small)
# res2 <- data %>%
#   mutate(ces1_perf= future_map(sensitivity,
#                                .f = ~estimate_performance_glmnet(sen_df = .x,ess = ces1_small)))



future::plan(multicore,workers=35)
# future::plan(sequential)
a <- Sys.time()
res2 <- data %>% 
  mutate(ces1_perf= future_map(sensitivity, 
                               .f = ~estimate_performance_glmnet(sen_df = .x,ess = ces1),
                               .options = furrr_options(seed = 0000))) %>%  
  mutate(ces2_perf= future_map(sensitivity,
                               .f = ~estimate_performance_glmnet(sen_df = .x,ess = ces2),
                               .options = furrr_options(seed = 0000))) %>%
  mutate(ceres_perf= future_map(sensitivity,
                                .f = ~estimate_performance_glmnet(sen_df = .x,ess = ceres),
                                .options = furrr_options(seed = 0000))) %>%
  mutate(demeter2_perf= future_map(sensitivity,
                                   .f = ~estimate_performance_glmnet(sen_df = .x,ess = demeter2),
                                   .options = furrr_options(seed = 0000)))
future::plan(sequential)
b <- Sys.time() 
b-a 
save(res2,file="/home/cloud-user/cluster_scratch/glmnet_modelling/result_largevm/res2_all4_266_365.RData")