#!/usr/bin/env Rscript
library(data.table)
require(gbm)

set.seed(2509)
# Read the data
test <- read.table("TA99Set",header=TRUE)

model <- readRDS("model_FINAL.rds")

pred_test <- data.frame(predict(model, n.trees=1000, newdata = test[,c(6, 3, 4, 7, 2, 11, 10, 5, 12, 14, 9, 8, 13)]))
test$PredictedActivity <- round(pred_test[,1],digits=3)
fwrite(test, file = "TA99Predictions.txt", sep="\t")
