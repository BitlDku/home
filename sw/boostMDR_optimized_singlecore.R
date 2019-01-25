###################################################################
# Boosting MDR
# single core - optimized version
# 2014.5.22
###################################################################

library(MDR)                   # for mdr
library(combinat)              # for combn 

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


#######################################################################
boostMDR.OptiS = function(ds, cl, dim=2,sRate=0.02) {
    ## Step. 1 evaluate pair of SNPs -------------------------------------
    nocol = ncol(ds)      
    PAIRxy <<- matrix(nrow=nocol, ncol=nocol)
 
    for (i in 1:(nocol-1)) {
        for (j in (i+1):nocol) {
           total = cbind(cl, ds[,i],ds[,j])
           colnames(total) = c("CL", "SNP1", "SNP2")
           
           tmpcase = total[total[,1]==0,]
           tmpcontrol = total[total[,1]==1,]
        
           diff.value = xtabs(~SNP1+SNP2, tmpcase)-xtabs(~SNP1+SNP2, tmpcontrol)
           diff.value = sum(abs(diff.value))/nrow(total)
           
           PAIRxy[i,j] <<- diff.value ; PAIRxy[j,i] <<- diff.value
        }  
    }

    ## Step. 2 evaluate combinations of SNPs -----------------------------------
    
    # 2.1 make combinations of SNPs
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
    
    # 2.2 perform mdr --------------------------------------
    
    whole = cbind(cl,ds)
    scombi = t(combi[,c(tmp)])
    rs = mdr(whole , scombi , x=5 , ratio = 1)
    
    return(rs)
}  


