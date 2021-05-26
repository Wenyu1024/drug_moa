library(caret)
rm(list=ls(all=TRUE))
# all = read.csv("secondary-screen-dose-response-curve-parameters.csv")

load("FK866.RData")
FK866 = FK866[which(duplicated(FK866$ccle_name)==F),] # remove duplicate

data = data.frame(FK866[, c(11:10578)])
y = FK866$ic50
y[is.na(y)] = mean(y, na.rm = T)


res = mat.or.vec(20,2)
i=1
for(i in 1:20){
  set.seed(i)
  print(i)
  trainIndex <- createDataPartition(y, p = 0.8, list = F) # 80% training 20% testing
  train.data = data[trainIndex,]
  test.data = data[-trainIndex,]
  cor_res = cor(y[trainIndex], train.data, use="complete.obs")
  select = which(cor_res > 0.15 | cor_res< -0.15) # focus only on the correlated genes
  
  
  trctrl <- trainControl(method = "cv", number = 5, verboseIter = F)
  fit1 = train(x = train.data[, select], y = y[trainIndex], method = "glmnet", trControl = trctrl, verbose = FALSE, trace = FALSE) 
  
  
  # fit1
  pred = predict(fit1, test.data[,select])
  
  res[i,1] = cor(pred, y[-trainIndex])
  # mean(abs(pred-y[-trainIndex]))
  # plot(pred, y[-trainIndex])
  
  # fit2 = train(x = train.data[, select], y = y[trainIndex], method = "avNNet", trControl = trctrl, verbose = FALSE, trace = FALSE) 
  fit2 = train(x = train.data[, select], y = y[trainIndex], method = "glmnet", trControl = trctrl, verbose = FALSE, trace = FALSE, tuneGrid = expand.grid(alpha = 0, lambda = seq(0,1,0.1))) 
  fit2
  pred2 = predict(fit2, test.data[,select])
  
  res[i,2] = cor(pred2, y[-trainIndex])
  # mean(abs(pred2 - y[-trainIndex]))
  # plot(pred2, y[-trainIndex])
  
  
  # fit3 = train(x = train.data[, select], y = y[trainIndex], method = "avNNet", trControl = trctrl, verbose = FALSE, trace = FALSE) 
  # fit3
  # pred3 = predict(fit3, test.data[,select])
  # 
  # res[i,3] = cor(pred3, y[-trainIndex])
  # mean(abs(pred3 - y[-trainIndex]))
  # plot(pred3, y[-trainIndex])
}

colnames(res) = c("glmnet","ridge")
res
