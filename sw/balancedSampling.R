######################################################
# Implement balanced sampling method
# 2018.02.19 BITL
# Parameters
#   ds   : target dataset 
#   cl   : class label vector
#   ts.n : ratio of test instances over whole instances 
#   k    : k of k-nearest neighbor
# return
#   indexes of training/test instances  
######################################################

library(FNN)               # for using get.knn
library(caret)
balanced.sampling <- function (ds, cl, ts.n=0.2, k=5) {
  classLabel <- unique(cl)
  ts.index <- c()
  
  # find k nearest neighbors for all instances
  knns <- get.knn(ds, k)      
  knns.idx <- data.frame(knns[["nn.index"]])
  # count different neighbors
  diff.neighbors <- c()
  for (j in 1:nrow(ds)) {
    this.class <- cl[j]
    neigbbors.class <- cl[c(as.integer(knns.idx[j,]))]
    diff.neighbors[j] <- length(which(neigbbors.class != this.class))
  }
  
  for (i in 1:length(classLabel)) {
    thisClass.idx = which(cl==classLabel[i])
    no = as.integer(length(thisClass.idx)*ts.n)
    if (no==0) no = 1
    this.diff.neighbors <- diff.neighbors[thisClass.idx]
    merged <- cbind(thisClass.idx, this.diff.neighbors)
    
    sampled.By.rvalue <- createDataPartition(y=merged[,2], p=ts.n)$Resample1
    ts.index <- c(ts.index, merged[sampled.By.rvalue,1])
  }
  
  tr.index <- setdiff(1:nrow(ds), ts.index)
  result <- list(tr=tr.index, ts=ts.index)
  return (result)
}

