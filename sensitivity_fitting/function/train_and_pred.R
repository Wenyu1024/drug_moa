library(tidyverse, quietly = T,lib.loc ="/projappl/project_2003466/project_rpackages")
library(tidymodels, quietly = T,lib.loc ="/projappl/project_2003466/project_rpackages")

lm_model <-
  linear_reg() %>%
  set_mode("regression") %>%
  set_engine("lm")

train_and_pred <- function(df){
  ceres_columns <- str_subset(string = colnames(df$training_data[[1]]), pattern = "_ceres")
  demeter_columns <- str_subset(string = colnames(df$training_data[[1]]), pattern = "_demeter")
  ces_columns <- str_subset(string = colnames(df$training_data[[1]]), pattern = "_ces")
  # if (is_empty(ceres_columns) | is_empty(demeter_columns) )  next
  formula_ceres <- as.formula(str_c("area_under_curve ~ ", str_c(ceres_columns, collapse = " + "), collapse = " "))
  formula_demeter2 <- as.formula(str_c("area_under_curve ~ ", str_c(demeter_columns, collapse = " + "), collapse = " "))
  formula_ces <- as.formula(str_c("area_under_curve ~ ", str_c(ces_columns, collapse = " + "), collapse = " "))
  # formula_all <- as.formula(area_under_curve ~ .)
  
  tmp <-   df %>% 
    mutate(fit_lm_ceres  = map(training_data, ~fit(lm_model, formula_ceres, data = .x))) %>% 
    mutate(fit_lm_demeter2  = map(training_data, ~fit(lm_model, formula_demeter2, data = .x))) %>% 
    mutate(fit_lm_ces = map(training_data, ~fit(lm_model,  formula_ces, data = .x))) %>% 
    mutate(model_pred_ceres = map2(fit_lm_ceres, testing_data, ~predict(.x, new_data = .y))) %>% 
    mutate(model_pred_demeter2 = map2(fit_lm_demeter2, testing_data, ~predict(.x, new_data = .y))) %>% 
    mutate(model_pred_ces = map2(fit_lm_ces, testing_data, ~predict(.x, new_data = .y))) %>% 
    mutate(res= pmap(list(testing_data, model_pred_ceres, model_pred_demeter2,model_pred_ces),
                     function(first, second, third,fourth){
                       tibble(real = first$area_under_curve, 
                              pred_ceres=second$.pred,
                              pred_demeter2=third$.pred,
                              pred_ces=fourth$.pred,                             
                       )}
    )) 
  
  
  get_perf <- function(x){
    x1= enquo(x)
    tmp %>%
      select(id, res) %>%
      unnest(res) %>%
      group_by(id) %>%
      metrics(truth = real, estimate = !!x1) %>%
      ungroup() %>%
      select(-3)
  }
  
  perf_ceres <- get_perf("pred_ceres") %>% rename(ceres = .estimate)
  perf_demeter2 <- get_perf("pred_demeter2") %>% rename(demeter = .estimate)
  perf_ces<- get_perf("pred_ces") %>% rename(ces = .estimate)
  
  res <-  perf_ceres %>%
    inner_join(perf_demeter2,by = c("id", ".metric")) %>%
    inner_join(perf_ces,by = c("id", ".metric"))
  return(res)
  
  # tmp <- tmp %>% mutate(ceres= perf_cer)
  # return(tmp)
}

