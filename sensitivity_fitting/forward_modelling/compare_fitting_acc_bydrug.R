library(furrr)
library(tidyverse)
library(tidymodels)
load("/scratch/project_2003466/forward_modelling/performance_testres/input.RData")

tictoc::tic()
plan(multicore,workers= 40)

# performance_difference_test_res1 <- combined_performance2 %>% 
#   mutate(test = furrr::future_map(
#     .x = data, 
#     .f = function(df){
#       df <- df %>% 
#         filter(input %in% c("ces", "ceres")) %>% 
#         pivot_wider(id_cols = one_of(c("id", "id2")),
#                     names_from= input, values_from=.estimate) %>% 
#         mutate(difference= ces- ceres) %>% 
#         t_test(response= difference,mu= 0)
#       return(df)
#       # t.test(x = df$ces, y = df$ceres,paired = TRUE,
#       #               alternative = "two.sided")
#     })) %>% 
#   # mutate(test= map(test,.f = tidy)) # %>%
#   select(-data) %>% 
#   unnest(test) 
# 
# performance_difference_test_res2 <- combined_performance2 %>% 
#   mutate(test = furrr::future_map(
#     .x = data, 
#     .f = function(df){
#       df <- df %>% 
#         filter(input %in% c("ces", "demeter2")) %>% 
#         pivot_wider(id_cols = one_of(c("id", "id2")),
#                     names_from= input, values_from=.estimate) %>% 
#         mutate(difference= ces- demeter2) %>% 
#         t_test(response= difference,mu= 0)
#       return(df)
#     })) %>% 
#   select(-data) %>% 
#   unnest(test) 

performance_difference_test_res1 <- combined_performance2 %>% 
  mutate(test = furrr::future_map(
    .x = data, 
    .f = function(df){
      df <- df %>% 
        filter(input %in% c("ces", "ceres")) %>% 
        pivot_wider(id_cols = one_of(c("id", "id2")),
                    names_from= input, values_from=.estimate) %>% 
        mutate(difference= ceres- ces) %>% 
        t_test(response= difference,mu= 0)
      return(df)
    })) %>% 
  select(-data) %>% 
  unnest(test) 

performance_difference_test_res2 <- combined_performance2 %>% 
  mutate(test = furrr::future_map(
    .x = data, 
    .f = function(df){
      df <- df %>% 
        filter(input %in% c("ceres", "demeter2")) %>% 
        pivot_wider(id_cols = one_of(c("id", "id2")),
                    names_from= input, values_from=.estimate) %>% 
        mutate(difference= ceres- demeter2) %>% 
        t_test(response= difference,mu= 0)
      return(df)
    })) %>% 
  select(-data) %>% 
  unnest(test) 

plan(sequential) 
tictoc::toc()

save(list = c("performance_difference_test_res1","performance_difference_test_res2"), 
     file = "/scratch/project_2003466/forward_modelling/performance_testres/output_ceresbased.RData")
