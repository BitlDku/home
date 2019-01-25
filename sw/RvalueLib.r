##############################################################################
## Rvale R(f)
## ds: data vector/matrix, cl:class vector, k:no of nearest neighbor,
## th: threshold
##############################################################################
Rvalue.f = function (ds, cl, k=7, th=3)
{
  library(FNN)               # for using get.knn
  rvalue = 0
  if (is.vector(ds)) {
    ds = cbind(ds)
  }
  if (is.vector(cl)) {
    cl = cbind(cl)
  }

  noRow = nrow(ds)           # no of Row

  knns = get.knn(ds, k)      # find k nearest neighbors for all instances
  knns.ndx <- data.frame(knns[["nn.index"]])

  overlap = c()
  for (i in 1:noRow) {
    thisClass = cl[i]
    error  = 0;
    for (j in 1:k) {
      if (cl[knns.ndx[i,j]] != thisClass) {
          error = error + 1
      }
    }
   
    if (error > th) {
      overlap[i] = 1          # current instance i is in correct class area
    } else {
      overlap[i] = 0          # current instance i is in overlapped area
    }
  }
  
  cntoverlap = sum(overlap)

  rvalue =  cntoverlap/noRow
  return(rvalue)

}

##############################################################################
## Rvale R(Ci)
## ds: data vector/matrix, cl:class vector, k:no of nearest neighbor,
## th: threshold
##############################################################################
Rvalue.Ci = function (ds, cl, k=7, th=4)
{

}

##############################################################################
## Rvale R(Ci, Cj)
## ds: data vector/matrix, cl:class vector, k:no of nearest neighbor,
## th: threshold
##############################################################################
Rvalue.CiCj = function (ds, cl, k=7, th=4)
{

}
