####################################################################################
# Feature Interaction
# based on prediction performance
# Author : Sejong Oh
# 2019.04.20
# Note : This code is not optimised 
####################################################################################
library(ggplot2)
library(forcats)      # for desc
library(RColorBrewer)

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
## Feature Impact (importance)
####################################################################################
f.impact <- function(model, Feature, get.seed=100, rep=10) {
  trainData <- model$trainingData
  ds <- trainData[,-ncol(trainData)]
  outcome <- trainData[,ncol(trainData)]
  ifelse(is.factor(outcome), task.type <- "Classification",
         task.type <- "Regression") 
  
  # calculate 
  pred.all <- predict(model, ds)
  if (task.type=="Classification") {      # classification
    PF.all <- mean(pred.all==outcome) 
  } else {                                # regression
    PF.all <- sqrt(sum((pred.all-outcome)^2)/length(outcome))  # rmse
  }

  PF.F1 <- c()
  for (k in 1:rep) {
    seed <- get.seed+k
    val.F1 <- ds[,Feature]
    val.F1 <- shupple(val.F1, seed)
    new.ds <- ds
    new.ds[,Feature] <- val.F1 
    pred.F1 <- predict(model, new.ds)
    if (task.type=="Classification") {    # classification
      PF.F1[k] <- mean(pred.F1==outcome) 
    } else {                              # regression
      PF.F1[k] <- sqrt(sum((pred.F1-outcome)^2)/length(outcome))  # rmse
    }
  }
  PF.F1 <- mean(PF.F1) 
  if (task.type=="Classification") { # classification
    FI.F1 <- PF.all - PF.F1
  } else {                           # regression
    FI.F1 <- PF.F1 - PF.all
  }
  return(FI.F1)
}
  
####################################################################################
# Feature Intraction
# model : prediction model produced by caret package
# F1    : name of first feature
# F2    : name of second feature
# rep   : repeat # for prediction test 
####################################################################################
f.interact <- function(model, F1, F2, rep=10) {
  trainData <- model$trainingData
  ds <- trainData[,-ncol(trainData)]
  outcome <- trainData[,ncol(trainData)]
  ifelse(is.factor(outcome), task.type <- "Classification",
         task.type <- "Regression") 

  # calculate 
  pred.all <- predict(model, ds)
  if (task.type=="Classification") {      # classification
     PF.all <- mean(pred.all==outcome) 
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
     if (task.type=="Classification") {  # classification
       PF.F1[k] <- mean(pred.F1==outcome) 
     } else {                            # regression
       PF.F1[k] <- sqrt(sum((pred.F1-outcome)^2)/length(outcome))  # rmse
     }

     seed <- 200+k
     val.F2 <- ds[,F2]
     val.F2 <- shupple(val.F2, seed)
     new.ds <- ds
     new.ds[,F2] <- val.F2 
     pred.F2 <- predict(model, new.ds)
     if (task.type=="Classification") {  # classification
       PF.F2[k] <- mean(pred.F2==outcome) 
     } else {                            # regression
       PF.F2[k] <- sqrt(sum((pred.F2-outcome)^2)/length(outcome))  # rmse
     }
     
     new.ds <- ds
     new.ds[,F1] <- val.F1 
     new.ds[,F2] <- val.F2 
     pred.F12 <- predict(model, new.ds)
     if (task.type=="Classification") {  # classification
       PF.F12[k] <- mean(pred.F12==outcome) 
     } else {                            # regression
       PF.F12[k] <- sqrt(sum((pred.F12-outcome)^2)/length(outcome))  # rmse
     }
  }
  
  PF.F1 <- mean(PF.F1)
  PF.F2 <- mean(PF.F2)
  PF.F12 <- mean(PF.F12)
  
  if(task.type=="Classification") {
    FI.F1 <- PF.all - PF.F1
    FI.F2 <- PF.all - PF.F2
    FI.F12 <- PF.all - PF.F12
    interact <- (FI.F1+FI.F2) - FI.F12
  } else {
    FI.F1 <- PF.F1 - PF.all  
    FI.F2 <- PF.F2 - PF.all  
    FI.F12 <- PF.F12 - PF.all 
    interact <- FI.F12 - (FI.F1+FI.F2)  
  }
  
  
  return(list(task=task.type,feature=c(F1, F2), fi=interact, 
              ass.all=PF.all, ass.F1=FI.F1, 
              ass.F2=FI.F2,ass.F12=FI.F12))  
} 


####################################################################################
# Print Result of f.interact 
####################################################################################
print.fi <- function(result){
  cat("** feature interaction between [",result$feature[1],
      "],[",result$feature[2], "] \n")
  cat("Task type                      : ", result$task, "\n")
  cat("Performance of whole dataset   : ", result$ass.all, "\n")
  cat("Reduced performance by F1      : ", result$ass.F1, "\n")
  cat("Reduced performance by F2      : ", result$ass.F2, "\n")
  cat("Reduced performance by {F1,F2} : ", result$ass.F12, "\n")
  cat("Interaction between {F1,F2}    : ", result$fi, "\n")
}


####################################################################################
## Feature Interaction Plot between F1 and others
####################################################################################
fi.single.barplot <- function(model, F1) {
  print("This function takes some minutes ...")
  trainData <- model$trainingData
  outcome <- trainData[,ncol(trainData)]
  features <- names(trainData)
  features <- setdiff(features, c(F1, ".outcome"))
  fi <- c()
  for (i in 1:length(features)) {
    fi[i] <- f.interact(model, F1, features[i])$fi 
  }
  df <- data.frame(features, fi)
  
  mymax <- max(abs(fi))
  if (is.factor(outcome)) {  # classification
    if (mymax <= 0.05) {
      mylimit <- 0.05 
    } else if (mymax <= 0.1) {
      mylimit <- 0.1 
    } else if (mymax <= 0.5) { 
      mylimit <- 0.5 
    } else {
      mylimit <- 1.0 
    }
  } else {                  # regression
    mylimit <- mymax * 1.1
  }
  df$features <- factor(df$features, levels = rev(features))
  p<- ggplot(data=df, aes(x=features, y=fi)) +
      geom_bar(stat="identity")+
      scale_y_continuous(limits = c(-1*mylimit,mylimit))
  p<- p+labs(title=paste("Feature Interaction for [ ", F1, " ] \n",sep=""), 
             y="degree of interaction")
  p+coord_flip()
} 


####################################################################################
## Create Feature Interaction Matrix 
####################################################################################
fi.matrix <- function(model) {
  print("This function may take long time ...")
  trainData <- model$trainingData
  features <- names(trainData)
  features <- setdiff(features, c(".outcome"))
  cnt <- length(features)
  fi <- matrix(NA, nrow=cnt, ncol=cnt)
  
  for (i in 1:(cnt-1)) {
    for (j in (i+1):cnt) {
      tmp <- f.interact(model, features[i], features[j])
      fi[i,j] <- tmp$fi 
      fi[j,i] <- fi[i,j]
    }
    fi[i,i] <- f.impact(model, features[i])
    cat(i, "/", length(features), " end..\n")
  }
  fi[cnt,cnt] <- f.impact(model, features[cnt])
  colnames(fi) <- features 
  rownames(fi) <- features  
  cat(cnt, "/", length(features), " end..\n")
  return(fi)
}


####################################################################################
## Feature Interaction Grid Chart 
####################################################################################
fi.grid.chart <- function(model, title="") {
  # print("This function may take long time ...")
  trainData <- model$trainingData
  features <- names(trainData)
  features <- setdiff(features, c(".outcome"))
  cnt <- length(features)
##
  fi <- fi.matrix(model)   # calculate feature interaction matrix 
  
  fi2 <- data.frame(fi)
  rownames(fi2) <- features
  colnames(fi2) <- features
  
  df <- data.frame(
    rep(features,cnt),
    rep(features, each=cnt)
  ) 
  colnames(df)<- c("x.f","y.f")
  df$x.f <- as.character(df$x.f)
  df$y.f <- as.character(df$y.f)
  
  for (i in 1:nrow(df)) { 
    df$z[i] <- round(fi2[df$x.f[i],df$y.f[i]],4) 
  }
  
  #define tile color

  col.pos <- rev(brewer.pal(9,"Blues")[1:6]) 
  col.neg <- c("#FFFF99","#FFFFCC")           # yellow
  
  
  fi.na <- fi2
  for(i in 1:cnt) {
    fi.na[i,i] <- NA
  }
  
  mymax <- max(fi.na, na.rm=T) ; mymin <- min(fi.na, na.rm=T) 

  for (i in 1:nrow(df)) { 
    
    if (df$x.f[i] == df$y.f[i]) {
      col.val = "#B5B5B5"          # gray
    } else if (df$z[i] > 0) {
      tmp <- max(6-round((df$z[i]/mymax)*6,0),1)
      col.val = col.pos[tmp]  
    } else if (df$z[i] == 0){
      col.val = "#FFFFFF"          
    } else  {
      tmp <- max(2-round((df$z[i]/mymin)*2,0),1)
      col.val = col.neg[tmp]
    }
    df$color[i] <- col.val 
    
  }
  
  ggplot(df, aes(x=x.f, y=fct_rev(y.f)))+
    geom_tile(fill=df$color)+
    geom_text(data=df, aes(x.f, y.f, label = z), size=rel(3))+
    scale_x_discrete(position="top")+ 
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 315))+
    ggtitle(paste("Feature Interaction Table : ", title, sep=""))+
    xlab("") + ylab("") 
  
} 


####################################################################################
## Calculate Total Feature Interaction of each feature and 
## draw plot
####################################################################################
fi.total.plot <- function(fi.matrix) {
  
  for (i in 1:nrow(fi.matrix)) {
    fi.matrix[i,i] <- 0  
  }
  
  df <- NULL 
  for (i in 1:nrow(fi.matrix)) {
    f.value <- fi.matrix[,i] 
    fi.pos <- sum(f.value[f.value>0])
    fi.neg <- sum(f.value[f.value<0])
    fi.tot <- fi.pos + abs(fi.neg)
    
    df <- rbind(df, c("pos", colnames(fi.matrix)[i] ,fi.pos))
    df <- rbind(df, c("neg", colnames(fi.matrix)[i] ,fi.neg))
    df <- rbind(df, c("tot", colnames(fi.matrix)[i] ,fi.tot))
  }
  
  df <- data.frame(df[,1:2],as.numeric(df[,3]))

  names(df) <- c("type","feature","degree_of_interaction")

  ggplot(data=df, aes(x=feature, y=degree_of_interaction, fill=type)) +
    geom_bar(stat="identity", position=position_dodge())+
    ggtitle("Total interaction Plot")
  
}


####################################################################################
## Find best pair of feature interactions  
####################################################################################
fi.best <- function(fi.matrix, n=5, type="pos") {
  for (i in 1:nrow(fi.matrix)) {
    fi.matrix[i,i] <- 0  
  }

  cnt <- ncol(fi.matrix)
  tbl <- matrix(0, nrow=(sum(1:(cnt-1))), ncol=3)
  idx = 1
  for (i in 1:(cnt-1)) {
    for (j in (i+1):cnt) {
      tbl[idx,1] = i
      tbl[idx,2] = j
      tbl[idx,3] = fi.matrix[i,j]
      idx <- idx + 1
    }    
  }
 
  if (type=="pos") {
    tbl <- tbl[order(tbl[,3], decreasing=T),]
    tbl <- subset(tbl, tbl[,3]>0)
    g.title <- "Best Interacted features (Positive) \n"
  } else {
    tbl <- tbl[order(tbl[,3]),]
    tbl <- subset(tbl, tbl[,3]<0)
    tbl[,3] <- abs(tbl[,3]) 
    g.title <- "Best Interacted features (Negative) \n"
  }
  
  if (nrow(tbl)==0) {
    Print("There is no candidate feature.")
    return()
  } else {
    n <- min(n, nrow(tbl))
  }
  
  ds <- tbl[1:n,3]
  nms <-c()
  for(i in 1:n) {
    name.F1 <- substr(colnames(fi.matrix)[tbl[i,1]],1,5)
    name.F2 <- substr(colnames(fi.matrix)[tbl[i,2]],1,5)
    
    nms[i] <- paste(name.F1, "-", name.F2, sep="") 
  }
  
  gds <- data.frame(ds,nms)
  gds$nms <- factor(gds$nms, levels = gds$nms[order(-gds$ds)])  

  
  p<- ggplot(data=gds, aes(x=nms, y=ds)) +
    geom_bar(stat="identity")
  p<- p+labs(title=g.title, x="pair of features", y="degree of interaction")
  p
}
