#' Import feature importance from CARET object to FItable object
#'
#' Get feature importance from CART to FItable
#'
#' THis function get feature importance made by 'filterVarImp' function in CARET package and
#' copy to FItable.
#'
#' @param caret.obj A prediction model (Classification or regression)
#' @param FItable.obj Training dataset (data frame) that is used to building model
#' @return [list] updated FItable object
#' @examples
#' library("FIG")
#' # for regression
#' data("Boston", package = "MASS")
#' model <- lm(medv ~ ., data = Boston)
#' FIobj <- FItable(model, train=Boston, target.name="medv", grid=50,
#'                  task="regression", interaction_type="OH2")
#' print(FIobj)
#' FIgraph(gtype="S", FIobj=FIobj, task="regression", seed=101)
#'
#' # get feature importance from CARET
#' myImp <- filterVarImp(x=Boston[,-14], y=Boston[,14])  #  14th column is 'mdev'
#' FIobj.new <- import_caret_imp(myImp, FIobj)
#' FIgraph(gtype="S", FIobj=FIobj.new, task="regression")
#'
#' @export
#'
import_caret_imp <- function(caret.obj, FItable.obj) {

  type <- ifelse(ncol(FItable.obj$Fimp)==2, "regression", "classification")
  FItable.obj.new = FItable.obj

  if (type=="regression") {
    FItable.obj.new$Fimp[,2] <- NA
    for (i in 1:nrow(FItable.obj.new$Fimp)) {
      idx <- which(rownames(caret.obj) == FItable.obj.new$Fimp$feature[i])
      FItable.obj.new$Fimp[i,2] <- caret.obj[idx,1]
    }

  } else {
    FItable.obj.new$Fimp[,3] <- NA
    for (i in 1:nrow(FItable.obj.new$Fimp)) {
      idx <- which(rownames(caret.obj) == FItable.obj.new$Fimp$feature[i])
      if (FItable.obj.new$Fimp$class[i] == "all") {
        FItable.obj.new$Fimp[i,3] <- rowMeans(caret.obj[idx,])
      } else {
        FItable.obj.new$Fimp[i,3] <- caret.obj[idx,FItable.obj.new$Fimp$class[i]]
      }
    }
  }

  return(FItable.obj.new)

}

