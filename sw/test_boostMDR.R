###################################################################
# Test boosting MDR 
# 
###################################################################
# Before run below code, change proper work directory 
#setwd("D:/mytest") 
 
## declare dataset --------------------------------------
library(mbmdr)
data(simSNP)   # it is in MDR package
cl = simSNP[,1]
ds = simSNP[,3:12]

## Test boostMDR.mormal() --------------------------------------
## easily understand the boostMDR algorithm
source("boostMDR_normal.R")
dim = 3            # dimension of features
sRate = 0.02       # search range in combinations of SNPs

mdrResult1 = boostMDR.normal(ds, cl, dim, sRate) 
print(mdrResult1$models)

## Test boostMDR.OptiS() --------------------------------------
## optimized and single core version
source("boostMDR_optimized_singlecore.R")
dim = 3            # dimension of features
sRate = 0.02       # search range in combinations of SNPs
noCore = 1         # no of core

mdrResult2 = boostMDR.OptiS(ds, cl, dim, sRate) 
print(mdrResult2$models)


## Test boostMDR.OptiM() --------------------------------------
## optimized and multicore version
source("boostMDR_optimized_multicore.R")
dim = 3            # dimension of features
sRate = 0.02        # search range in combinations of SNPs
noCore = 1         # no of core

mdrResult3 = boostMDR.OptiM(ds, cl, dim, sRate, noCore) 
print(mdrResult3$models)

## Feature selection (filtering) --------------------------------------
## boostMDR contains simple feature selection function based on mutual information
## parpameter 'select' assigns the number of features that we select
newds = mi.filter(ds, cl, select=100)