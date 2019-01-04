setwd("c:\\works")

#rm(list = ls())
install.packages("FNN") # for balanced sampling
install.packages("caret") # for balanced sampling
install.packages("class") # for knn

source("balancedSampling.R")

file = "dataset.csv"
whole.ds = read.csv(paste(file,sep=""))
ds = whole.ds[,-1]
cl = whole.ds[,1]
   
# balanced sampling
idx = balanced.sampling(ds,cl, ts.n=0.25, k = 6)


TRAIN.DS = ds[idx$tr,]
TRAIN.CL = cl[idx$tr]
TEST.DS = ds[idx$ts,]
TEST.CL = cl[idx$ts]
#knn ---------------------------------
rs = knn(TRAIN.DS, TEST.DS, TRAIN.CL, k=5)
knn.acc <- mean(rs == TEST.CL)

cat("Knn :",knn.acc, "\n")








