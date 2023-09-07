library(parallel)
library(ranger)
library(xgboost)
library(caret)
library(pROC)

get_model_result <- function(ds, cl, test.idx){
  model <- xgboost(ds[-test.idx,], cl[-test.idx], 
                   nthread = 4, 
                   nrounds = 5, 
                   num_class=length(unique(cl)), 
                   objective = "multi:softprob", 
                   verbose = 0)
  
  prob <- predict(model, ds[test.idx,])
  prob <- prob[seq(2,length(prob),2)]
  
  res <- eval_metric(prob, cl[test.idx])
  return(res)
}

evaluate_model <- function(ds, cl, k=3, SEED=100, rep=1) {
  if(is.vector(ds)) {
    ds <- data.frame(ds)
  }
  ds <- as.matrix(ds)
  cl <- as.integer(factor(cl)) -1
  
  res.rep <- c()
  
  for(r in 1:rep){
    #k-fold cross validation
    set.seed(SEED+r)
    res <-avg <- std <- c()
    
    if (k==1) {
      set.seed(SEED+r)
      test.idx <- sample(1:nrow(ds), nrow(ds)*.3)
      res <- get_model_result(ds, cl, test.idx)
      
    } else {
      fold <- createFolds(cl, k, list = TRUE)
      for(i in 1:k) {
        test.idx <- fold[[i]] 
        res.cur <- get_model_result(ds, cl, test.idx)
        res <- rbind(res, res.cur)
      }
    }
    res.rep <- rbind(res.rep, res)
  }
  avg <- colMeans(res.rep, na.rm = T)
  std <- colSums((res.rep - avg)^2, na.rm = T)/nrow(res.rep)
  res <- list(mean = avg, sd =std)

  return(res)  
} 


eval_metric <- function(prob, y){
  y <- factor(y)
  pred <- ifelse(prob > 0.5, 1, 0)
  pred <- factor(pred, levels = levels(y))
  suppressMessages({
    cm <- confusionMatrix(pred, y)  
    AUC <- auc(roc(y, prob))  
  })
  
  res <- c(cm$overall['Accuracy'],
           cm$byClass['Sensitivity'],
           cm$byClass['Specificity'],
           cm$byClass['F1'],
           AUC)
  
  names(res) <- c('accuracy','sensitivity','specificity','fl','auc')
  
  return(res)
}


# Filtering
gradual_permutation_filtering <- function(ds, cl){
  cnt <- 50
  filtered.set <- colnames(ds)
  filtered.len <- 0
  
  for(ratio in seq(0.1,0.9,0.1)){ 
    rank <- c()
    ds.new <- cbind(cl, ds[,filtered.set])
    
    for(i in 1:cnt){
      set.seed(100+i)
      imp <- ranger(dependent.variable.name=colnames(ds.new)[1], 
                    data=ds.new, 
                    importance = 'permutation')
      rank <- rbind(rank, imp$variable.importance)
    }
    idx <- apply(rank, 2, function(x){sum(x > 0) >= cnt*ratio})
    neg <- names(idx)[!idx]
    
    rank <- colMeans(rank)
    rank[neg] <- 0
    idx <- order(rank, decreasing = T)
    ranking <- names(rank[idx])
    filtered.set <- ranking[rank[idx] > 0]
    cat(length(filtered.set),"of", ncol(ds)-1, "features are selected.\n")
    
    if(length(filtered.set) == filtered.len) break
    filtered.len <- length(filtered.set)
  }
  
  return(filtered.set)
}


##############3
fitness <- function(performance, theta, selected_n, N){
  theta*performance + (1-theta)*(1-log(selected_n, base = N))
}


heuristic_tribrid_search <- function (filtered, cl, i, N) {
  RATIO <- 0.5
  fold <- 5
  tradeoff <- 0.005
  theta <- log(2,base = N)/{tradeoff+log(2,base = N)}
  
  RANKED <- colnames(filtered)
  chosen <- RANKED[i]
  candidate <- setdiff(RANKED, chosen)
  
  
  # First conduct a forward search. 
  # If there is no performance improvement, execute the Consolation method.
  
  stop.flag <- F
  consol <- 0
  max.fit <- 0
  del.fs <- c()
  dummy <- rep(1, nrow(filtered))
  
  print("Forward Search")
  while(stop.flag == F) {
    if (length(candidate)==0) break 
    
    end <-  max(1, as.integer(length(candidate)* RATIO))
    
    # Forward search
    consol.flag <- T
    stop.flag <- T
    max.fs <- NULL
    
    for (k in 1:end) {
      performance <- evaluate_model(cbind(dummy, filtered[,c(chosen, candidate[k])]), cl, fold)$mean['auc']
      fit <- fitness(performance, theta, length(c(chosen, candidate[k])), N)
      
      if(fit > max.fit) {
        consol.flag <- F
        stop.flag <- F
        max.fit <- fit
        max.fs <- candidate[k]
      }
    }
    
    chosen <- c(chosen, max.fs)
    candidate <- setdiff(candidate, max.fs)
    
    if(consol.flag == T & length(chosen)>2){
      stop.flag <- T
      
      if(end > 1){
        consol.obj <- consolate(filtered, cl, chosen, candidate, max.fit, end, theta, fold)
        chosen <- consol.obj[[1]]
        candidate <- consol.obj[[2]]
        consol.fit <- consol.obj[[3]]
        consol <- consol + consol.obj[[4]]
        
        if (consol.fit > max.fit) {
          stop.flag <- F
          max.fit <- consol.fit 
          max.fs <- chosen
          cat("Consolation Match=", consol, max.fit,  "\n")
        } 
      } 
    }
    
    if(stop.flag == F) cat('Selected: ', chosen, 'fit = ', max.fit,'\n')
  }    
  
  ## backward elimination ########################################
  print("Backward Elimination")
  while(length(chosen) > 1){
    back.fs <- NA
    for (sel in 1:length(chosen)) {
      performance <- evaluate_model(cbind(dummy, filtered[, chosen[-sel]]), cl, k=5)$mean['auc']
      fit <- fitness(performance, theta, length(c(chosen, candidate[k])), N)
      if (fit >= max.fit) {
        max.fit <- fit
        back.fs <- chosen[sel]
        print(back.fs)
      }
    }
    if (!is.na(back.fs)) {
      chosen <- setdiff(chosen, back.fs)
      cat(back.fs, "is removed \n")
    }else break
  }
  
  cat('Selected: ', chosen, 'fit = ', max.fit,'\n')
  
  return(list(features=chosen, fit=max.fit))
}

consolate <- function(ds, cl, chosen, candidate, max.fit, end, theta, fold){
  #consolation match #########################################
  in.fs <- NA
  out.fs <- NA
  
  consol <- 0
  c.st <- ceiling(length(chosen)/2)
  c.end <- length(chosen)-1
  
  for (k in c.st:c.end) { 
    for (m in 0:end) {    
      if(m == 0) {
        tmp <- c(chosen[-k])
      }else{
        tmp <- c(chosen[-k], candidate[m])    
      }
      performance <- evaluate_model(ds[,tmp], cl, fold)$mean['auc']
      fit <- fitness(performance, theta, length(c(chosen, candidate[k])), ncol(ds))
      
      if(fit > max.fit) {
        max.fit <- fit
        in.fs <- ifelse(m != 0, candidate[m], "del")
        out.fs <- chosen[k]
        cat("in & out : ", m, ',', k,'\n')
      }
    }
  }
  
  if (!is.na(in.fs)) {
    if(in.fs == "del"){
      chosen <- setdiff(chosen, out.fs)
      candidate <- c(out.fs, candidate)   
      print("consol-del")
    }else{
      chosen <- c(setdiff(chosen, out.fs), in.fs)  
      candidate <- c(out.fs, setdiff(candidate, in.fs))    
    }
    consol <- 1
  }
  
  return(list(chosen, candidate, max.fit, consol))
}

heuristic_tribrid_search_parallel <- function(i, filtered, cl, N) {
  
  # For parallel computation
  numCores <- detectCores()
  clst <- makeCluster(numCores)
  
  clusterExport(clst, "createFolds")
  clusterExport(clst, "xgboost")
  clusterExport(clst, "confusionMatrix")
  clusterExport(clst, "auc")
  clusterExport(clst, "roc")
  
  clusterExport(clst, "get_model_result")
  clusterExport(clst, "evaluate_model")
  clusterExport(clst, "eval_metric")
  clusterExport(clst, "fitness")
  
  clusterExport(clst, "gradual_permutation_filtering")
  clusterExport(clst, "heuristic_tribrid_search")
  clusterExport(clst, "heuristic_tribrid_search_parallel")
  clusterExport(clst, "consolate")
  
  
  best.set <- NULL
  best.num <- ncol(filtered)
  best.acc <- 0
  best.start <- 0
  
  res <- heuristic_tribrid_search(filtered, cl, i, N) 
  
  stopCluster(clst)
  
  return(list("fit" = res$fit, "features" = res$features))
}


hybrid_feature_selection <- function (ds, cl) {
  #RANKED, theta, RATIO = 0.5, fold=5, start=1
  
  #gradual permutation filtering
  print("Now graudal permutation filtering  ...")
  rank <- gradual_permutation_filtering(ds, cl)  
  filtered <- ds[,rank]
  
  # Heuristic Tribrid Search: This section can be processed in parallel.
  print("Now Heuristic Tribrid Search  ...")
  N <- ncol(ds)
  num <- min(length(rank),10)
  #final_set <- heuristic_tribrid_search(filtered, cl, 1, N)
  candidate_set <- parLapply(clst, 1:num, heuristic_tribrid_search_parallel, filtered, cl, N)
  
  fit <- c()
  for(m in 1:num){fit[m] <- candidate_set[[m]]$fit}
  best.start <- which(fit == max(fit))[1]
  best.set <- candidate_set[[best.start]]$features
  best.fit <- candidate_set[[best.start]]$fit
  
  return(list(feature=best.set, fit=best.fit))
}

