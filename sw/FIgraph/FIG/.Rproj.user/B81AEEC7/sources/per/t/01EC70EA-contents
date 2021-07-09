################################################################################
# Feature Interaction Graph
#
# Author : Sejong Oh
# 2021.07.01
# Note : This code is not optimized
################################################################################
#' @import igraph
#' @import vip
#' @import R6
#' @import checkmate
#' @import data.table
#' @import Metrics

# library("igraph")
# library("fastshap")
# library("vip")
# library("R6")
# library("checkmate")
# library("data.table")

# source("InterpretationMethod.R")     # from iml package
# source("MarginalGenerator.R")        # from iml package
# source("Predictor.R")                # from iml package
# source("find_y.R")                   # from iml package
# source("Data.R")                     # from iml package
# source("utils.R")                    # from iml package
# source("inferTask.R")                # from iml package
# source("create_predict_fun.R")       # from iml package
# source("FeatureImp.R")               # from iml package
# source("Shapley.R")                  # from iml package
# source("Interaction.R")              # from iml package (modified)
# source("f_interact_lib.R")


###############################################################################
# Create Feature interaction & Feature importance table
###############################################################################

## Feature interaction method
Friedman <- function(model, train, target.name, grid=50, task="regression") {
  mod <- NULL
  if (task == "regression") {
    mod <- Predictor$new(model, data = train[,-which(names(train) == target.name)])
  } else {  # classification
    mod <- Predictor$new(model, data = train[,-which(names(train) == target.name)],
                         type="prob", y=train[,target.name] )
  }
  ##
  cbn = combn(setdiff(names(train), target.name), 2)
  var1 <- c();  var2 <- c(); class <- c(); fi <- c()
  print("Calculate feature interaction ...")
  if (task=="regression") {
    for(i in 1:ncol(cbn)) {
      var1[i] <- cbn[1,i]
      var2[i] <- cbn[2,i]
      ia <- Interaction$new(mod, feature=var1[i], feature2=var2[i], grid.size=grid)
      qq=ia$results
      fi[i] = qq$.interaction
      cat(var1[i], var2[i], fi[i], class[i], "\n")
    }
    result <- data.frame(from=var1,to=var2, weight=fi)
    result$weight <- round(result$weight,3)
  } else {   # classification
    tmpdf <- NULL
    for(i in 1:ncol(cbn)) {
      var1[i] <- cbn[1,i]
      var2[i] <- cbn[2,i]
      ia <- Interaction$new(mod, feature=var1[i], feature2=var2[i], grid.size=grid)
      tmpdf = rbind(tmpdf, ia$results)
      print(ia$results); print(" ")
    }
    f1f2 <- as.data.frame(do.call(rbind, strsplit(tmpdf$.feature, ":")))
    result <- data.frame(class=tmpdf$.class, from=f1f2[,1],to=f1f2[,2], weight=tmpdf$.interaction)
    result <- result[order(result$class),]
  }
  result$weight <- round(result$weight,3)

  return(result)
}

## Feature interaction method ##############################################
Greenwell <- function(model, train, target.name, grid=50, task="regression") {

  train2 <- train[,-which(names(train)==target.name)]
  cbn = combn(setdiff(names(train), target.name), 2)
  var1 <- c();  var2 <- c(); class <- c(); fi <- c()
  print("Calculate feature interaction ...")

  if (task=="regression") {
    for(i in 1:ncol(cbn)) {
      var1[i] <- cbn[1,i]
      var2[i] <- cbn[2,i]
      ia <- vint(model, feature_names = c(var1[i], var2[i]), train=train2,
                 type="regression")
      fi[i] = data.frame(ia)[1,2]
      cat(var1[i], var2[i], fi[i], class[i], "\n")
    }
    result <- data.frame(from=var1,to=var2, weight=fi)
    result$weight <- round(result$weight,3)
  } else { # classification
    idx <-1
    for(c in unique(train[,target.name])){
      for(i in 1:ncol(cbn)) {
        class[idx] <- c
        var1[idx] <- cbn[1,i]
        var2[idx] <- cbn[2,i]
        ia <- vint(model, feature_names=c(var1[idx], var2[idx]), train=train2,
                   type="classification", pprob=T, which.class=c)
        fi[idx] <- data.frame(ia)[1,2]
        cat( class[idx],var1[idx], var2[idx], fi[idx], "\n")
        idx <- idx+1
      }
    }
    result <- data.frame(class=class, from=var1, to=var2, weight=fi)
    result$weight <- round(result$weight,3)
  }
  return(result)
}

## Feature interaction method ###############################################
Oh2 <- function(model, train, target.name, all.class=TRUE, grid=50, task="regression") {

  train2 <- train[,-which(names(train)==target.name)]
  cbn = combn(setdiff(names(train), target.name), 2)
  var1 <- c();  var2 <- c();  fi <- c(); class <- c();
  print("Calculate feature interaction ...")

  if (task=="regression") {
    cl= train[, target.name]
  } else {  # classification
    cl= factor(train[, target.name])
  }

  if (all.class) {
      for(i in 1:ncol(cbn)) {
        class[i] <- "all"
        var1[i] <- cbn[1,i]
        var2[i] <- cbn[2,i]
        ia <- f.interact(model, train=train2, cl= cl,
                         F1= var1[i],F2=var2[i], rep=10)
        fi[i] = ia$fi
        cat(class[i], var1[i], var2[i], fi[i],  "\n")
      }
  } else {  # FI separated by class
    idx <- 1
    labs <- levels(cl)
    for (lb in labs) {
      for(i in 1:ncol(cbn)) {
        class[idx] <- lb
        var1[idx] <- cbn[1,i]
        var2[idx] <- cbn[2,i]
        ia <- f.interact(model, train=train2, cl= cl,
                         F1= var1[i],F2=var2[i], clabel=lb, rep=10)
        fi[idx] = ia$fi
        cat(class[idx], var1[idx], var2[idx], fi[idx],  "\n")
        idx <- idx +1
      }
    }
  }
  result <- data.frame(class=class, from=var1,to=var2, weight=fi)
  result$weight <- round(result$weight,3)
  if (task=="regression") {
    result <- result[,-1] # remove class column
  }

  return(result)
}

###############################################################################
#
# task             : "regression" (default), "classification"
# interaction_type : "Friedman" (default), "Greenwell", "OH2"
# importance_type  : "permute" (default), "firm"
#
# FItable <- function(model, train, target.name, grid=50, task="regression",
#                     interaction_type="Friedman", all.class=TRUE){
#
#   train2 <- train[,-which(names(train)==target.name)]
#   cl <-  train[,target.name]
#
#   ## Get pairwise interaction value ------------------------------------------
#   if (interaction_type=="FRIEDMAN") {
#     result <- Friedman(model, train, target.name, grid, task)
#   } else if (interaction_type=="GREENWELL") {
#     result <- Greenwell(model, train, target.name, grid, task)
#   } else if (interaction_type=="OH2") {
#     result <- Oh2(model, train, target.name,  all.class=all.class, grid, task)
#   }
#
#   ## Get feature importance --------------------------------------------------
#
#   if (task=="regression") {
#
#     if (importance_type=="permute") {
#       mod <- Predictor$new(model, data = train2, y = cl)
#       set.seed(100)
#       imp <- FeatureImp$new(mod, loss = "mae")
#       myimp <- imp$results[,c(1,3)]
#     }
#     names(myimp) <- c("feature", "importance")
#
#   } else { # classification
#
#     #if (importance_type=="permute") {
#        if(all.class) {
#          mod <- Predictor$new(model, data = train2, y = cl, type = "prob")
#          set.seed(100)
#          imp <- FeatureImp$new(mod, loss = "ce")
#          myimp <- data.frame(rep("all", ncol(train2)),imp$results[,c(1,3)])
#        } else {
#          labs <- levels(cl)
#          merge.myimp = NULL
#          for (lb in labs) {
#            set.seed(100)
#            mod <- Predictor$new(model, data = train2, y = cl, type = "prob", class=lb)
#            imp <- FeatureImp$new(mod, loss = "ce")
#            myimp <- data.frame(rep(lb, ncol(train2)),imp$results[,c(1,3)])
#            merge.myimp <- rbind(merge.myimp, myimp)
#          } # end for
#          myimp <- merge.myimp
#        } # end else if
#     #} else  if (importance_type=="caret") {
#     #  # caret varImp
#     #
#
#     #}
#
#
#
#     names(myimp) <- c("class","feature", "importance")
#   }
#
#   return(list(Fint=result, Fimp=myimp))
# }


# ##############################################################################
# # Main function
# FIgraph <- function(th=NULL, gtype="S", FIobj, task="regression", class.name=NULL,
#                     show.edge.weight=TRUE) {
#
#   # Get Feature interation  & feature importance
#   result <- FIobj$Fint   # feature interaction table
#   myimp  <- FIobj$Fimp   # feature importance table
#
#   # filter low interaction
#   if(is.null(th)) {
#     th <- mean(result$weight)
#     if (!is.null(class.name)) th <- mean(result$weight[result$class==class.name])
#   }
#   result <- result[result$weight>=th,]
#   cat("th: ", th, "\n")
#
#   ## Draw interaction Graph #############################################
#
#   if (task=="regression") {
#     edges <- result
#     nodes <- data.frame(myimp)
#   } else {
#     edges <- result[result$class==class.name,-1]
#     nodes <- data.frame(myimp[myimp$class==class.name,-1])
#   }
#
#   net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
#
#   # node size
#   impsize <- (V(net)$importance - min(V(net)$importance)) /
#     (max(V(net)$importance)-min(V(net)$importance))
#   impsize <- ((impsize/1.2)+0.2) *50
#   V(net)$size <- impsize # V(net)$importance*9  # node size
#
#   # node color
#   pal1 <- heat.colors(nrow(nodes), alpha=1)
#   idx <- (V(net)$importance - min(V(net)$importance)) /
#     (max(V(net)$importance)-min(V(net)$importance))
#   idx <- nrow(nodes)-as.integer(idx*(nrow(nodes)-1)+1)+1
#   V(net)$color <- pal1[idx]
#
#
#   E(net)$label <- edges$weight     # edge weight
#   #E(net)$label.cex=1
#   #E(net)$label.font=2
#
#   # edge width
#   ewidth <- (E(net)$weight - min(E(net)$weight)) /
#     (max(E(net)$weight)-min(E(net)$weight))
#   ewidth <- (ewidth + 0.2)*7
#   E(net)$width <- ewidth # edges$weight*20  # edge width
#
#   # Draw graph
#
#   title <- ""
#   if (task=="classification") title <- paste0("class: ", class.name)
#
#   set.seed(101)
#   if(!show.edge.weight) {
#     elabel <- NA
#   } else {
#     elabel <- E(net)$value
#   }
#
#   if (gtype=="S") {  # static graph
#     plot(net, edge.arrow.size=.4, vertex.label=V(net)$nodes, edge.label=elabel,
#          main=title)
#   } else {            # interactive graph
#     tkplot(net, edge.label=elabel)
#   }
#
# }
#

################################################################################
# test
#
# rm(list=ls())
# setwd("D:/OneDrive - ?ܱ????б?/????/idea_feature_interaction_network/iml_oh")
#
# data("Boston", package = "MASS")
# library("randomForest")
# rf <- randomForest(medv ~ ., data = Boston, ntree = 50)
#
# FIobj1 <- FItable(rf, train=Boston, target.name="medv", grid=50, task="regression")
# FIgraph(th=0.1, gtype="S", FIobj=FIobj1, task="regression")
#
# library(e1071); library(pdp)
# data(pima)
# #sv <- svm(iris[,-5], iris$Species, probability = T)
# pima <- pima[complete.cases(pima),]
# sv <- svm(pima[,-9], pima$diabetes)
#
# FIobj2 <- FItable(sv, train=iris, target.name="Species", grid=50, task="classification",
#                   interaction_type="Greenwell", importance_type="shap")
# FIgraph(th=NULL, gtype="S", FIobj=FIobj2, task="classification", class.name="setosa",
#         show.edge.weight=T)
#
# source("f_interact_lib.R")
# FIobj3 <- FItable(sv, train=iris, target.name="Species", all.class=F, grid=50, task="classification",
#                   interaction_type="OH2")
# FIgraph(th=NULL, gtype="S", FIobj=FIobj2, task="classification", class.name=NULL,
#         show.edge.weight=T)
#
