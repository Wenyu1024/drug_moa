data <- data.table::fread("/scratch/project_2003466/ces_io/ces_input_20q4.csv")

#library(dplyr)
# install.packages("speedglm")
#cells <- unique(data$DepMap_ID)
#data1 <- data %>% dplyr::filter(DepMap_ID %in% cells[1:20])

ptm <- proc.time()
lm <- lm(demeter2~ ceres+mut+exp_seq+cn+exp_array+DepMap_ID ,data = data)
pred <- predict.lm(lm)
tmp = proc.time() - ptm
print(tmp)
save(list= c( "pred"), file="ces2.RData")
write.csv(pred, "ces2_pred.csv", row.names=F)

