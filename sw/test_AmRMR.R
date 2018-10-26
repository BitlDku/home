
#_____________________________________________________________________________________#
setwd(choose.dir()) # set AmRMR source file, R-value source file, isolet dataset and this file in chosen directory.
source("AmRMR.R")
library("class")
normalize <- function(x, class = FALSE, cl = ncol(x)){ 
  ds=x
  if(class == TRUE){
    ds = ds[, -cl]
  }
  nor = as.data.frame(lapply(ds, function(x){return((x- min(x))/(max(x) - min(x)))}))
  return(nor)
}

samplize <- function(x, is_valid = FALSE, p = c(4, 1)){
  ds        <- x
  nds <- nrow(ds)
  partition <- c("train", "test")
  idx       <- c()
  rtn       <- c()
  test <- data.frame(); train <- data.frame(); valid <- data.frame();
  if(is_valid == FALSE){
    idx     <- sample(x = partition, size = nds, replace = TRUE, prob = p)
    test    <- ds[idx == "test", ]
    train   <- ds[idx == "train", ]
    rtn     <- list(train=train, test=test)
    return(rtn)
  }else{
    partition = c(partition, "valid")
    idx     <- sample(x = partition, size = nds, replace = TRUE, prob = p)
    test    <- ds[idx == "test", ]
    train   <- ds[idx == "train", ]
    valid   <- ds[idx == "valid", ]
    rtn     <- list(train=train, test=test, valid = valid)
    return(rtn)
  }
}
#_____________________________________________________________________________________#


# step 1. get target data and set seed value.
target = read.csv("isolet.csv")
set.seed(521)

# Let me analysis test dataset using KNN algorithm.
# step 2. KNN analysis without AmRMR Feature Selection.
temp = samplize(x = target, is_valid = FALSE, p = c(3, 1))

train.ds = normalize(temp$train[, -1])
train.cl = as.factor(temp$train[, 1])
test.ds  = normalize(temp$test[, -1])
test.cl  = as.factor(temp$test[, 1])

nor.knn = knn(train = train.ds, test = test.ds, cl = train.cl, k = 5)
nor.acc = mean(nor.knn == test.cl)

# step 3. KNN analysis with AmRMR.
# AmRMR will choose 25 optimal features 
target.20 = runMRMR.boost(ds = target[, -1], cl = target[, 1], cnt = 25, method = "MIQ")
View(target[,c(1, target.20)])
temp = samplize(x = target[,c(1, target.20)], is_valid = FALSE, p = c(3, 1))
train.ds.20 = normalize(temp$train[, -1])
train.cl.20 = as.factor(temp$train[, 1])
test.ds.20  = normalize(temp$test[, -1])
test.cl.20  = as.factor(temp$test[, 1])

AmRMR.knn = knn(train = train.ds.20, test = test.ds.20, cl = train.cl.20, k = 5)
AmRMR.acc = mean(AmRMR.knn == test.cl.20)

#step 4. show result
cat("KNN Accuracy without AmRMR : ", nor.acc, "\n", 
    "KNN Accuracy with AmRMR    : ", AmRMR.acc, "\n", sep = "")


