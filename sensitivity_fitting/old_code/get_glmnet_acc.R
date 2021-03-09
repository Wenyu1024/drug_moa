library(furrr)
library(tidyverse)
library(tidymodels)
source("get_fixed_glmnet.R")
# sen_df refers to the cellular sensitivity data of a drug


estimate_performance_glmnet <- function(sen_df, ess, inside_parallel=F){
    df <- sen_df  %>%   as_tibble() %>% 
      inner_join( 
        y= janitor::clean_names(ess)%>% 
          rename_at(vars(!contains("dep_map_id")), .fun= ~paste0(., "_predictors")), 
        by= c("DepMap_ID"="dep_map_id")) %>% 
      select(-master_ccl_id, -apparent_ec50_umol, -DepMap_ID) 
    
    # if (ess= "ceres + demeter2") {df1 <- }
    # if (ess= "ces1 + ces2")  {df1 <- }
    

    # set.seed(0000)
    df1 <-  df %>% 
      nested_cv(outside=vfold_cv(v = 5), inside= vfold_cv(v = 5)) %>% 
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
        cor(x = x %>% pull(.pred), y = y %>% pull(area_under_curve), method = "spearman")}))
      res= df1 %>% pull(spearman_cor)
      return(df1)
     
}





# tmp <- estimate_performance_glmnet(sen_df = tmp_bottom$sensitivity, ess = ces1[,1:100])


# load("~/drug_target/data/ctrpv2_data.RData")
# load("~/drug_target/data/new_ess_modeling.RData")
source('~/drug_target/analysis/get_fixed_glmnet.R')

# parallel with purrr http://zevross.com/blog/2019/02/12/dramatically-speed-up-your-r-purrr-functions-with-the-furrr-package/
set.seed(0000)
future::plan(multiprocess)
data <- ctrpv2_data %>%
  ungroup() %>%
  slice(1:2) 

a <- Sys.time()
res2 <- data %>% 
  mutate(ces1_perf= future_map(sensitivity, 
                        .f = ~estimate_performance_glmnet(sen_df = .x,ess = ces1))) %>% 
  mutate(ces2_perf= future_map(sensitivity,
                        .f = ~estimate_performance_glmnet(sen_df = .x,ess = ces2))) %>%
  mutate(ceres_perf= future_map(sensitivity,
                        .f = ~estimate_performance_glmnet(sen_df = .x,ess = ceres))) %>%
  mutate(demeter2_perf= future_map(sensitivity,
                        .f = ~estimate_performance_glmnet(sen_df = .x,ess = demeter2)))


# tmp1 <- ctrpv2_data$sensitivity[[1]]
# tmp <- estimate_performance_glmnet(sen_df = tmp1, ess = ces1[,1:100])
b <- Sys.time() 
b-a #Time difference of 5.273595 mins
save.image("~/drug_target/analysis/all_drug_forward_modeling2.RData")


#Fold3: internal: A correlation computation is required, but `estimate` is constant 
# It happens when your model produces a single predicted value. R2 needs the variance (which is then zero) and produces an NA value.
# note not all fold return this result so it doesn\t matter, since the problematic fold don\t influence choosing the best model


# load("~/drug_target/analysis/all_drug_forward_modeling1.RData")

tmp <- res %>% select(broad_cpd_id,ces1_perf:demeter2_perf) %>% 
  unnest() %>% 
  pivot_longer(cols= -"broad_cpd_id",names_to= "ess_type",values_to= "spear_cor") 

tmp1 <- tmp  %>% 
  group_by(broad_cpd_id, ess_type) %>% 
  summarize(ACC= mean(spear_cor)) %>% 
  ungroup() %>% filter(ess_type== "ces1_perf") 




tmp %>% ggplot2::ggplot(mapping = aes(x= broad_cpd_id, y = spear_cor, color= ess_type)) +
  geom_boxplot(position= position_dodge(width = 1))+
  geom_point(position=position_dodge(width = 1)) +
  coord_flip() +
  theme_classic()


tmp %>% ggplot2::ggplot(mapping = aes(x= ess_type, y = spear_cor)) +
  geom_boxplot(position= position_dodge(width = 1))+
  theme_classic()
# idea, uncertainty varies across drugs. maybe give an uncertainty estimate for the identified drug target?

# Time difference of 1.790067 hours for 2 drugs on 4 core pc
# Time difference of 1.784935 hours for 10 drugs on 40 cores 5*5 nestcv                                            
# Time difference of 16.17045 hours for 165 drugs on 40 cores 5'5 nestcv
# Time difference of 9.5164 hours for 150 drugs on 40 cores 5*5 nestcv.
# warning to be confirmed :




