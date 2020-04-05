#!/usr/bin/env Rscript
library(data.table)
require(gbm)

set.seed(2509)
# Read the data
test <- read.table("TB254Set",header=TRUE)

model <- readRDS("model_FINAL.rds")

pred_test <- data.frame(predict(model, n.trees=1000, newdata = test[,c(6, 3, 4, 7, 2, 11, 10, 5, 12, 14, 9, 8, 13)]))
test$PredictedActivity <- pred_test[,1]
fwrite(test, file = "TBB254Predictions.txt", sep="\t")
