###################################################################
# Boosting MDR
# understandable version
# 2014.5.22
###################################################################

library(MDR)                 # for mdr 
library(combinat)            # cor combn 

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
   for(i in 1:combi.snp.length) {
      S1 = ds.subset[combi.snp[i,1]]
      S2 = ds.subset[combi.snp[i,2]]
   
      total = cbind(cl, ds[,S1],ds[,S2])
      colnames(total) = c("CL", "SNP1", "SNP2")
        
      tmpcase = total[total[,1]==0,]
      tmpcontrol = total[total[,1]==1,]
        
      diff.value = xtabs(~SNP1+SNP2, tmpcase)-xtabs(~SNP1+SNP2, tmpcontrol)
      diff.value = sum(abs(diff.value))/nrow(total)   

      C = C + diff.value
   }
   m.value =  C

   return (m.value)
}


#######################################################################
boostMDR.normal = function(ds, cl, dim=2,sRate=0.02) {
    nocol = ncol(ds)
    
    # make combinations of SNPs
    combi = combinat::combn(c(1:nocol), dim) 
    my.select = as.integer(ncol(combi)*sRate)        # no of selected best feature
    if (my.select < min(50, ncol(combi))) {
        my.select = min(50, ncol(combi))
    }
    length.combi = ncol(combi)
    
    result = c()
    # recommand good pair of features -------------------
    for ( m in 1:length.combi) {
        fs = c(combi[,m])  # combination of feature
        result[m] = PCT(c(fs))
    }              
    tmp = order(result, decreasing = TRUE)[1:my.select]   
    #combi[,c(tmp)]
    
    ## perform mdr --------------------------------------
    
    whole = cbind(cl,ds)
    scombi = t(combi[,c(tmp)])
    rs = mdr(whole , scombi , x=5 , ratio = 1)
    
    return(rs)
}  

