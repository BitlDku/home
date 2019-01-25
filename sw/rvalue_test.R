#####################################################################
# Test R value R(f)
#####################################################################

noFold=5
noFS=40

setwd("D:/Users/Lee/Desktop/RFS_rproject")
source("RvalueLib.r")

rawData = read.csv("madelon.csv", header=FALSE)

noRow = nrow(rawData)     # number of rows
noCol = ncol(rawData)     # number of columns
noFeat = noCol-1;

# cross validation
library(cvTools)

# creat fold list
folds<-cvFolds(noRow, K=noFold)

acc=0
for (f in 1:noFold) {
  # creat data set
  train <- rawData[folds$subsets[folds$which != f], ]
  test <- rawData[folds$subsets[folds$which == f], ]
  
  # Write modeling and evaluation code here.
  trainCls = train[,1]          # class info
  trainData = train[,2:noCol]

  rval=c() # declare result vector
  for (i in 1:noFeat) {
    rval[i] = Rvalue.f(cbind(trainData[,i]), cbind(trainCls),k=7, th=3) 
  }
  ranker=sort(rval, method="sh", index.return=TRUE)$ix[c(1:noFS)]

  testCls = test[, 1]
  noTest = length(testCls)
  #testData

  #select feature
  newTrainData <- train[, ranker]
  newTestData <- test[, ranker]

  predict=knn(newTrainData, newTestData, trainCls, k=7)
  #print(predict)
  hit=0
  for(i in 1:noTest) {
    if(predict[i]==testCls[i]) {
      hit=hit+1
    }
  }
  tmpAcc=0
  tmpAcc=hit/noTest
  print(tmpAcc)
  
  acc=acc+tmpAcc
}

acc=acc/noFold
cat("Total Accuracy: ", acc)