setwd("/scratch/project_2003466/ces_21q1_io/")
library(tidyverse)
data <- read_csv("ces_input_21q1.csv")

# install.packages("speedglm")
#cells <- unique(data$DepMap_ID)
#data1 <- data %>% dplyr::filter(DepMap_ID %in% cells[1:20])

# ptm <- proc.time()
lm <- lm(ceres~ demeter2+mut+exp_seq+cn+exp_array+DepMap_ID ,data = data)
pred <- predict.lm(lm)
# tmp = proc.time() - ptm
# print(tmp)
data <- data %>% 
  select_at(c(1:4,6)) %>% 
  mutate(ces1= pred) 
  
ces1_437 <- data %>% select(DepMap_ID, gene, ces1) %>% pivot_wider(id_cols = DepMap_ID, names_from= gene, values_from= ces1)
ceres_437 <- data %>% select(DepMap_ID, gene, ceres) %>% pivot_wider(id_cols = DepMap_ID, names_from= gene, values_from= ceres)
demeter2_437 <- data %>% select(DepMap_ID, gene, demeter2) %>% pivot_wider(id_cols = DepMap_ID, names_from= gene, values_from= demeter2)
expseq_437 <- data %>% select(DepMap_ID, gene, ceres) %>% pivot_wider(id_cols = DepMap_ID, names_from= gene, values_from= expseq_437)

write_csv(ces1_437,path = "ces_added_longdf_21q1.csv")
write_csv(ces1_437,path = "ces_added_longdf_21q1.csv")
write_csv(ces1_437,path = "ces_added_longdf_21q1.csv")
write_csv(ces1_437,path = "ces_added_longdf_21q1.csv")