setwd("/scratch/project_2003466/ces_io/")
load("ces_prepare_21q1_intermediate.RData")
library(tidyverse)
df_long_scaled_noarray <- 
  wide_to_long(ceres, "ceres",cells1) %>% 
  inner_join(wide_to_long(demeter2, "demeter2",cells1),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(exp_seq, "exp_seq",cells1),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(mut, "mut",cells1) ,by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(cn, "cn",cells1) ,by=c("DepMap_ID", "gene")) %>% 
  drop_na() %>% 
  group_by(DepMap_ID) %>% 
  mutate_at(.vars = 3:7,.funs = scale) %>% 
  ungroup()


df_long_scaled <- 
  wide_to_long(ceres, "ceres",cells2) %>% 
  inner_join(wide_to_long(demeter2, "demeter2",cells2),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(exp_array, "exp_array",cells2),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(exp_seq, "exp_seq",cells2),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(mut, "mut",cells2) ,by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(cn, "cn",cells2) ,by=c("DepMap_ID", "gene"))%>% 
  drop_na() %>% 
  group_by(DepMap_ID) %>% 
  mutate_at(.vars = 3:8,.funs = scale) %>% 
  ungroup()


write_csv(df_long_scaled, "ces_input_21q1.csv")
write_csv(df_long_scaled_noarray,"ces_input_noarray_21q1.csv")
write_csv(ceres, "ceres_21q1.csv")
write_csv(demeter2, "demeter2_21q1.csv")
save.image("ces_prepare_21q1.RData")

