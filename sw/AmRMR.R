#######################################################
# MRMR BOOST implement 
# Developed by Bio-IT Lab. in Dankook University, Korea
# Verson 1.0 (2017.01.28) 
#########################################################
require("infotheo")      # mutual infroamtion, entropy
source("RvalueLib2.R")
##============================================================================
## claculate block entropy for given single feature
## input : ds=dataset, cl=class label, cnt=no of chosen feature,
##         method="MID" or "MIQ", prt="Y" or "N"
## ds shold have descrete values            
##============================================================================
runMRMR.boost = function (ds, cl, cnt, method="MID", prt="N") {
    fea = c()
    noCol = ncol(ds)
    KMAX = min(1000,noCol)
    cnt = min(KMAX,cnt)
    
    #calcualte MI 
    rds = c(1:noCol)
    names(rds) = c(1:noCol)

    for (i in 1:noCol) {
      rds[i] = Rvalue.f(ds[,i], cl)
    }
    rds = sort(rds, decreasing = TRUE)
    rds = rds[1:KMAX]
    idx = as.numeric(names(rds))
    
    fea[1] = idx[1]
#    rds = rds[2:KMAX]
    idx = idx[-1]
    
    if (prt == "Y") {
        cat("Feature 1 :", fea[1], "\n") 
    }
    
    for (i in 2:cnt) {
        length.fea = length(fea)
        length.idx = length(idx)
  
        mrmr = c()
        for (j in 1:length.idx) {
            t_mi = Rvalue.f(cbind(ds[,c(fea)],ds[,idx[j]]), cl,k=5)   # relavance            
            cor_array = c()
#           for (k in 1:length.fea) {
              thisIdx = idx[j]

# PLAN A-------------------              
#              tmp = abs(cor(ds[,fea[k]], ds[,thisIdx]))
#              if (tmp > 0.5) {
#                cor_array[k] = tmp  
#              } else {
#                cor_array[k] = 0
#              }
#            }
#            c_mi =  max(cor_array)                               # redundancy
# -------------------------
# PLAN B-------------------              
#              cor_array = abs(cor(ds[,fea], ds[, thisIdx], method = "pearson"))
#              cor_value = cor_array[cor_array>=0.5]
#              c_mi = mean(cor_value)
#              if(is.nan(c_mi)){ ## NaN process
#                c_mi = mean(cor_array)
#              }
# -------------------------
# PLAN C-------------------              
              cor_array = abs(cor(ds[,fea], ds[, thisIdx], method = "pearson"))
              c_mi = max(cor_array)
              if(c_mi < 0.5){ ## NaN process
                c_mi = mean(cor_array)
              }
# -------------------------
# PLAN D-------------------              
#              cor_array = abs(cor(ds[,fea], ds[, thisIdx], method = "pearson"))
#              c_mi = mean(cor_array)
# -------------------------
              
            if (method == "MID") {    
               mrmr[j] = t_mi - c_mi
            } else if (method == "MIQ") {
               if ((t_mi / c_mi)>t_mi) {
                  mrmr[j] = t_mi 
               } else {
                  mrmr[j] = t_mi / c_mi
               }
               
            }
        }
        maxIdx = min(which(mrmr == max(mrmr)))
        if (prt == "Y") {
            cat("Feature", i, ":", idx[maxIdx], "\n") 
        }
        fea[i] = idx[maxIdx]
        rds = rds[-maxIdx]           # remove chosen feature     
        idx = idx[-maxIdx]         
    } 
    
    return (fea)
}
