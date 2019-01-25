###################################################################
# Boosting MDR
# multicore-optimized version
# 2014.5.22
###################################################################

library(MDR)                   # for mdr
library(combinat)              # for combn 
library(parallel)              # for parallel processing on multi core cpu

PAIRxy = NULL

#######################################################################
mi.filter = function (ds, cl, select) {
    nocol = ncol(ds)
    if (select >= nocol) {
       print ("'select' should be less than number columns of 'ds'")
       return (ds)
    }

    mi.value = c()
    for (i in 1:nocol) {
       mi.value[i] = infotheo::mutinformation(ds[,i], cl, method="emp") 
    }
    tmp = order(mi.value, decreasing = TRUE)[1:select]
    
    return(ds[,c(tmp)])
}


#######################################################################
PCT = function (ds.subset) {
   nocol = length(ds.subset)

   # sum of SNP pair value
   combi.snp = t(combinat::combn(c(1:nocol), 2))
   combi.snp.length = nrow(combi.snp) 
   C = 0                                    # sum of SNP pair 
   for(x in 1:combi.snp.length) {
      C = C + PAIRxy[ds.subset[combi.snp[x,1]],ds.subset[combi.snp[x,2]]]
   }
   m.value =  C

   return (m.value)
}

#############################################################

boostMDR.OptiM = function(ds, cl, dim=2,sRate=0.02, noCore=1) {
    ## Step. 1 evaluate pair of SNPs -------------------------------------
    nocol = ncol(ds)      

    PAIRxy <<- matrix(nrow=nocol, ncol=nocol)
    cPair = NULL
    for (i in 1:(nocol-1)) {
      for (j in (i+1):nocol) {
         cPair = rbind(cPair, c(i,j))
      }
    }
    cPair.length = nrow(cPair)
    
    parallel.PAIRxy = function(m) {
       i = cPair[m,1]; j = cPair[m,2]          
       total = cbind(cl, ds[,i],ds[,j])
       colnames(total) = c("CL", "SNP1", "SNP2")
       
       tmpcase = total[total[,1]==0,]
       tmpcontrol = total[total[,1]==1,]
    
       diff.value = xtabs(~SNP1+SNP2, tmpcase)-xtabs(~SNP1+SNP2, tmpcontrol)
       diff.value = sum(abs(diff.value))/nrow(total)
    }
    ss = simplify2array(mclapply(1:cPair.length, parallel.PAIRxy, mc.cores=noCore))       # MCAPPLY
    for (m in 1:cPair.length) {
      i = cPair[m,1]; j = cPair[m,2]          
      PAIRxy[i,j] <<- ss[m] ; PAIRxy[j,i] <<- ss[m]
    }

    ## Step. 2 evaluate combinations of SNPs -----------------------------------
    
    # 2.1 make combinations of SNPs
    combi = combinat::combn(c(1:nocol), dim) 
    my.select = as.integer(ncol(combi)*sRate)        # no of selected best feature
    if (my.select < min(50, ncol(combi))) {
        my.select = min(50, ncol(combi))
    }
    length.combi = ncol(combi)
    
    # recommand good pair of features -------------------
    parallel.PCT = function(m) {
        fs = c(combi[,m])  # combination of feature
        cur.PCT = PCT(c(fs))
    }              
    result= simplify2array(mclapply(1:length.combi, parallel.PCT, mc.cores=noCore))   # MCAPPLY
    tmp = order(result, decreasing = TRUE)[1:my.select]   
    
    # 2.2 perform mdr --------------------------------------
    
    whole = cbind(cl,ds)
    scombi = t(combi[,c(tmp)])
    rs = mdr(whole , scombi , x=5 , ratio = 1)
    
    return(rs)
}  
