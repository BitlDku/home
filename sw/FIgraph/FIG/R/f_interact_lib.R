####################################################################################
# Feature Interaction
# based on prediction performance
# Author : Sejong Oh
# 2019.04.20
# modified :2021.07.01
# Note : This code is not optimized
# ####################################################################################

####################################################################################
## Shuple featue values
####################################################################################
shupple <- function(Feature, seed=100) {
  half <- as.integer(length(Feature)/2)

  set.seed(seed)
  shup <- sample(1:length(Feature), replace=F)
  return(Feature[shup])

}

####################################################################################
# Feature Interaction
# model  : prediction model
# train  : training dataset without label
# cl     : label info (should be factor type)
# F1     : name of first feature
# F2     : name of second feature
# clabel : NULL or specific class label when classification
# rep    : repeat # for prediction test
####################################################################################
f.interact <- function(model, train, cl, F1, F2, clabel=NULL, rep=10) {
  ds <- train
  outcome <- cl

  ifelse(is.factor(outcome), task.type <- "classification",
         task.type <- "regression")

  # calculate
  pred.all <- predict(model, ds)
  if (task.type=="classification") {      # classification
    if (is.null(clabel)) {
      PF.all <- mean(pred.all==outcome)
    } else {
      PF.all <- mean(pred.all[cl==clabel] == outcome[cl==clabel])
    }
  } else {                                # regression
    PF.all <- sqrt(sum((pred.all-outcome)^2)/length(outcome))  # rmse
  }

  PF.F1 <- c(); PF.F2 <- c(); PF.F12 <- c()
  for (k in 1:rep) {
    seed <- 100+k
    val.F1 <- ds[,F1]
    val.F1 <- shupple(val.F1, seed)
    new.ds <- ds
    new.ds[,F1] <- val.F1
    pred.F1 <- predict(model, new.ds)

    if (task.type=="classification") {  # classification
      if (is.null(clabel)) {
        PF.F1[k] <- mean(pred.F1==outcome)
      } else {
        PF.F1[k] <- mean(pred.F1[cl==clabel] == outcome[cl==clabel])
      }
    } else {                            # regression
      PF.F1[k] <- sqrt(sum((pred.F1-outcome)^2)/length(outcome))  # rmse
    }

    seed <- 200+k
    val.F2 <- ds[,F2]
    val.F2 <- shupple(val.F2, seed)
    new.ds <- ds
    new.ds[,F2] <- val.F2
    pred.F2 <- predict(model, new.ds)
    if (task.type=="classification") {  # classification
      if (is.null(clabel)) {
        PF.F2[k] <- mean(pred.F2==outcome)
      } else {
        PF.F2[k] <- mean(pred.F2[cl==clabel]==outcome[cl==clabel])
      }
    } else {                            # regression
      PF.F2[k] <- sqrt(sum((pred.F2-outcome)^2)/length(outcome))  # rmse
    }

    new.ds <- ds
    new.ds[,F1] <- val.F1
    new.ds[,F2] <- val.F2
    pred.F12 <- predict(model, new.ds)
    if (task.type=="classification") {  # classification
      if (is.null(clabel)) {
        PF.F12[k] <- mean(pred.F12==outcome)
      } else {
        PF.F12[k] <- mean(pred.F12[cl==clabel]==outcome[cl==clabel])
      }
    } else {                            # regression
      PF.F12[k] <- sqrt(sum((pred.F12-outcome)^2)/length(outcome))  # rmse
    }
  }

  PF.F1 <- mean(PF.F1)
  PF.F2 <- mean(PF.F2)
  PF.F12 <- mean(PF.F12)

  if(task.type=="classification") {
    FI.F1 <- PF.all - PF.F1
    FI.F2 <- PF.all - PF.F2
    FI.F12 <- PF.all - PF.F12
  } else {   # regression
    FI.F1 <- PF.F1 - PF.all
    FI.F2 <- PF.F2 - PF.all
    FI.F12 <- PF.F12 - PF.all
  }

  interact <- FI.F12 - max(FI.F1, FI.F2)

  return(list(task=task.type,feature=c(F1, F2), fi=interact,
              ass.all=PF.all, ass.F1=FI.F1,
              ass.F2=FI.F2,ass.F12=FI.F12))
}

