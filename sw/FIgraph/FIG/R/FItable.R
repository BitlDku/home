#' Create a feature interaction table
#'
#' This function create a feature interaction table.
#'
#' @param model A prediction model (Classification or regression)
#' @param train Training dataset (data frame) that is used to building model
#' @param target.name Name of target label name in train dataset
#' @param grid Number of grid to calculate partial dependency function. Default is 50)
#' @param task Prediction task. "regression" (default) or "classification".
#' @param interaction_type Type of feature interaction algorithms. One of "OH2" (default), "FRIEDMAN","GREENWELL".
#' @param all.class When task is classification, Feature importance can be calculated for overall classed (all.class==TRUE) or for each class (all.class==FALSE).
#' @return [list] feature interaction & feature importance table
#' @examples
#' library("FIG")
#' # for regression
#' data("Boston", package = "MASS")
#' model <- lm(medv ~ ., data = Boston)
#' FIobj1 <- FItable(model, train=Boston, target.name="medv", grid=50,
#'                  task="regression", interaction_type="OH2")
#' print(FIobj1)
#'
#' # for classification
#' library(e1071)
#' model2 <- svm(Species~., data=iris)
#' FIobj2 <- FItable(model2, train=iris, target.name="Species", grid=50,
#'                  task="classification", interaction_type="OH2", all.class=F)
#' print(FIobj2)
#'
#' @export
FItable <- function(model, train, target.name, grid=50, task="regression",
                    interaction_type="OH2", all.class=TRUE){
  # error handling
  if (!(interaction_type %in% c("FRIEDMAN", "GREENWELL", "OH2"))) {
    print("interaction_type should be one of 'FRIEDMAN', 'GREENWELL', 'OH2' ")
    return(NULL)
  }
  if (!(task %in% c("regression", "classification"))) {
    print("task should be one of regression', 'classification'")
    return(NULL)
  }

  train2 <- train[,-which(names(train)==target.name)]
  cl <-  train[,target.name]

  ## Get pairwise interaction value ------------------------------------------
  if (interaction_type=="FRIEDMAN") {
    result <- Friedman(model, train, target.name, grid, task)
  } else if (interaction_type=="GREENWELL") {
    result <- Greenwell(model, train, target.name, grid, task)
  } else if (interaction_type=="OH2") {
    result <- Oh2(model, train, target.name,  all.class=all.class, grid, task)
  }

  ## Get feature importance --------------------------------------------------

  if (task=="regression") {

   # if (importance_type=="permute") {
      mod <- Predictor$new(model, data = train2, y = cl)
      set.seed(100)
      imp <- FeatureImp$new(mod, loss = "mae")
      myimp <- imp$results[,c(1,3)]
   # }
    names(myimp) <- c("feature", "importance")

  } else { # classification

    #if (importance_type=="permute") {
    if(all.class) {
      mod <- Predictor$new(model, data = train2, y = cl, type = "prob")
      set.seed(100)
      imp <- FeatureImp$new(mod, loss = "ce")
      myimp <- data.frame(rep("all", ncol(train2)),imp$results[,c(1,3)])
    } else {
      labs <- levels(cl)
      merge.myimp = NULL
      for (lb in labs) {
        set.seed(100)
        mod <- Predictor$new(model, data = train2, y = cl==lb, type = "prob", class=lb)
        imp <- FeatureImp$new(mod, loss = "ce")
        myimp <- data.frame(rep(lb, ncol(train2)),imp$results[,c(1,3)])
        merge.myimp <- rbind(merge.myimp, myimp)
      } # end for
      myimp <- merge.myimp
    } # end else if
    #} else  if (importance_type=="caret") {
    #  # caret varImp
    #

    #}



    names(myimp) <- c("class","feature", "importance")
  }

  return(list(Fint=result, Fimp=myimp))
}

