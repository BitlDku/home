##############################################################################
## Rvale R(f)
## ds: data vector/matrix, cl:class vector, k:no of nearest neighbor,
## th: threshold
##############################################################################
library(FNN)               # for using get.knn
library(stats)             # for using table

#############################################################################
Rvalue.f = function (ds, cl, k=5)
{
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

  inOverlap = c()

  diffNeighbor  = 0;
  for (i in 1:noRow) {
    thisClass = cl[i]
    for (j in 1:k) {
      if (cl[knns.ndx[i,j]] != thisClass) {
          diffNeighbor = diffNeighbor + 1
      }
    }
  }

  rvalue =  diffNeighbor/(noRow*k)

  return(1 - rvalue)

}
#############################################################################
Rvalue.f.attrs = function (gattrs, k=5)
{
  ds = GLOBAL.DS[,c(gattrs)]
  cl = GLOBAL.CL
  
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

  inOverlap = c()

  diffNeighbor  = 0;
  for (i in 1:noRow) {
    thisClass = cl[i]
    for (j in 1:k) {
      if (cl[knns.ndx[i,j]] != thisClass) {
          diffNeighbor = diffNeighbor + 1
      }
    }
  }

  rvalue =  diffNeighbor/(noRow*k)
  cat(length(gattrs), 1-rvalue, "----\n")
  
  return(1 - rvalue)

}

