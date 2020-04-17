#####################################################################################
# hybrid RFE  
# 2020-04-15
# R version 3.6.0
# Bio-Information Technology Lab.
# 
# This code references multiple svmRFE on
#       https://github.com/johncolby/SVM-RFE
#
######################################################################################

set.seed(12345)
lib = list('randomForest', 'e1071', 'caret', 'gbm')
unlist(lapply(lib, require, character.only = T))

if (!require('randomForest')) install.packages('randomForest')
if (!require('e1071')) install.packages('e1071')
if (!require('caret')) install.packages('caret')
if (!require('gbm')) install.packages('gbm')

library('randomForest')
library('e1071')
library('caret')
library('gbm')


#####################################################################################
##  kfold HybridRFE             
## 
##  input parameters:
## 
##    X       : dataset. first column should be factor type class label  
##    max.row : maximum number of rows for RandomForest modeling  
##    k       : number of fold for RFE experiment
##    halve.above : see https://github.com/johncolby/SVM-RFE
##    method  : ('s.sum'- simple sum function,  'w.sum'-weighted sum function
#####################################################################################

hybridRFE.kfold <- function(X, k=5, halve.above=500, max.row=1000, method = c('s.sum', 'w.sum')) {
    set.seed(1234)
    folds = caret::createFolds(y = X[,1], k = 5, list = T)

    X[,1] = as.numeric(as.factor(X[,1]))
    results = lapply(folds, hybridRFE.wrap, X, k, halve.above, max.row, method)

    top.features = WriteFeatures(results, X, save=F)
    rank = top.features$FeatureID

    return(rank)
}

#####################################################################################
##  hybridRFE.wrap                                                  
#####################################################################################

hybridRFE.wrap <- function(test.fold, X, k, halve.above, max.row, method) {
    # Wrapper to run hybridRFE function while omitting a given test fold
    train.data = X[-test.fold, ]
    test.data  = X[test.fold, ]
    
    # Rank the features
    features.ranked = hybridRFE(train.data, k, halve.above, max.row, method)
    
    return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
}


#####################################################################################
##  normalize Fun.                                                 
#####################################################################################
normalize <- function(x) { 
    z=x
    if(min(x)<max(x)){ 
        z=(x - min(x)) / (max(x) - min(x))
    }
    return(z)
}  


#####################################################################################
#  Run HybridRFE Function              
#####################################################################################
hybridRFE <- function(X, k, halve.above, max.row, method) {
    # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
    # X = input.ds
    n = ncol(X) - 1
    
    # Scale data up front so it doesn't have to be redone each pass
    cat('Scaling data...')
    X[, -1] = normalize(X[, -1])
    cat('Done!\n')
    flush.console()
    
    pb = txtProgressBar(1, n, 1, style=3)
    
    i.surviving = 1:n
    i.ranked    = n
    ranked.list = vector(length=n)
    
    # Recurse through all the features
    while(length(i.surviving) > 0) {
        if(k > 1) {
            # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)
            set.seed(100)
            # folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
            # folds = lapply(1:k, function(x) which(folds == x))
            folds = caret::createFolds(y = X[,1], k = k, list = T)   #####
            
            
            # Obtain weights for each training set
            w = lapply(folds, getWeights.hybrid, X[, c(1, 1+i.surviving)], max.row, method)
            w = do.call(rbind, w)
            
            # Normalize each weights vector
            w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
            
            
            # Compute ranking criteria
            # v    = w * w
            # vbar = apply(v, 2, mean) # ????�� ???? / ?????? ????
            # vsd  = apply(v, 2, sd) # ?????? ?л? 
            # c    = vbar / vsd
            
            c = colMeans(w)
            
        } else {
            # Only do 1 pass (i.e. regular SVM-RFE)
            w = getWeights.hybrid(NULL, X[, c(1, 1+i.surviving)], max.row, method)      # weight value for each feature
            # print(w)           # OH
            # c = w * w
            c = w
        }
        
        # Rank the features
        ranking = sort(c, index.return=T)$ix
        if(length(i.surviving) == 1) {
            ranking = 1
        }
        
        if(length(i.surviving) > halve.above) {
            # Cut features in half until less than halve.above
            nfeat = length(i.surviving)
            ncut  = round(nfeat / 2)
            n     = nfeat - ncut
            
            cat('Features halved from', nfeat, 'to', n, '\n')
            flush.console()
            
            pb = txtProgressBar(1, n, 1, style=3)
            
        } else ncut = 1
        
        # Update feature list
        ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
        i.ranked    = i.ranked - ncut
        i.surviving = i.surviving[-ranking[1:ncut]]
        
        if (i.ranked == 1) {
            ranked.list[1] <- i.surviving[1]
            break        # exit while loop
        }
        
        setTxtProgressBar(pb, n-length(i.surviving))
        flush.console()
    } # end while
    
    setTxtProgressBar(pb, n)
    flush.console()
    close(pb)
    
    return (ranked.list)
}



#####################################################################################
##  Calculate feature weight        
#####################################################################################

getWeights.hybrid <- function(test.fold, X, max.row, method) {
    # Fit a linear SVM model and obtain feature weights
    train.data = X
    # test.fold = folds[[2]]
    if(!is.null(test.fold)) train.data = X[-test.fold, ];  test.data = X[test.fold, ]
    
    svmModel = svm(train.data[, -1], train.data[, 1], cost=10, cachesize=500,
                   scale=F, type="C-classification", kernel="linear")
    
    svm.pred = predict(svmModel, test.data[,-1])
    svm.acc = mean(svm.pred == test.data[,1])
    
    w.svm = t(svmModel$coefs) %*% svmModel$SV
    w.svm = w.svm*w.svm
    w.svm <- (w.svm-min(w.svm))/(max(w.svm)-min(w.svm)) # normalize
    
    # [,drop = T]
    if(nrow(train.data)> max.row) {
        set.seed(1234)
        sel <- sample(nrow(train.data), max.row)
        train.data.RF <- train.data[sel,]
    } else {
        train.data.RF <- train.data
    }
    
    rfModel <- randomForest(train.data.RF[,-1], as.factor(train.data.RF[,1]), importance=TRUE, proximity=TRUE)
    
    rf.pred = predict(rfModel, test.data[,-1])
    rf.acc = mean(rf.pred == as.factor(test.data[,1]))
    
    w.rf <- rfModel$importance[,"MeanDecreaseAccuracy"]
    w.rf <- (w.rf-min(w.rf))/(max(w.rf)-min(w.rf)) # normalize
    
    
    #----
    names(train.data)[1] = 'group'
    names(test.data)[1] = 'group'
    train.data[,1] = factor(train.data[,1])
    test.data[,1] = factor(test.data[,1])
    
    gbmModel <- gbm(group ~ ., data=train.data,distribution = "multinomial")
    gbm.pred = predict(gbmModel, test.data[,-1], n.trees=100, type="response")
    labels = apply(gbm.pred, 1, which.max) - 1
    gbm.acc = mean(labels == test.data[,1])
    
    tmp <- summary(gbmModel, plotit=F)
    idx = c()
    for(i in 1:ncol(train.data[,-1])) {
        idx[i] = which(tmp$var == names(train.data[,-1])[i])
    }
    
    w.gbm <- tmp$rel.inf[idx]
    w.gbm <- (w.gbm-min(w.gbm))/(max(w.gbm)-min(w.gbm)) # normalize
    # ----
    
    if(method == 's.sum'){
        w.sum <- w.svm + w.rf + w.gbm
    } else {
        w.sum <- (w.svm * svm.acc) + (w.rf * rf.acc) + (w.gbm * gbm.acc)  
    }
    return(w.sum)
}

#####################################################################################
##  calculate mean weight from k-fold 
#####################################################################################

WriteFeatures <- function(results, input, save=F, file='features_ranked.txt') {
    # Compile feature rankings across multiple folds
    featureID = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$ix
    avg.rank  = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$x
    feature.name = colnames(input[, -1])[featureID]
    features.ranked = data.frame(FeatureName=feature.name, FeatureID=featureID, AvgRank=avg.rank)
    if(save==T) {
        write.table(features.ranked, file=file, quote=F, row.names=F)
    } else {
        features.ranked
    }
}

FeatSweep.wrap <- function(i, results, input) {
    # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
    svm.list = lapply(results, function(x) tune(svm,
                                                train.x      = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                                train.y      = input[x$train.data.ids, 1],
                                                validation.x = input[x$test.data.ids, 1+x$feature.ids[1:i]],
                                                validation.y = input[x$test.data.ids, 1],
                                                # Optimize SVM hyperparamters
                                                ranges       = tune(svm,
                                                                    train.x = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                                                    train.y = input[x$train.data.ids, 1],
                                                                    ranges  = list(gamma=2^(-12:0), cost=2^(-6:6)))$best.par,
                                                tunecontrol  = tune.control(sampling='fix'))$perf)
    
    error = mean(sapply(svm.list, function(x) x$error))
    return(list(svm.list=svm.list, error=error))
}

PlotErrors <- function(errors, errors2=NULL, no.info=0.5, ylim=range(c(errors, errors2), na.rm=T), xlab='Number of Features',  ylab='10x CV Error') {
    # Makes a plot of average generalization error vs. number of top features
    AddLine <- function(x, col='black') {
        lines(which(!is.na(errors)), na.omit(x), col=col)
        points(which.min(x), min(x, na.rm=T), col='red')
        text(which.min(x), min(x, na.rm=T), paste(which.min(x), '-', format(min(x, na.rm=T), dig=3)), pos=4, col='red', cex=0.75)
    }
    
    plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
    AddLine(errors)
    if(!is.null(errors2)) AddLine(errors2, 'gray30')
    abline(h=no.info, lty=3)
}

########################################################################
## Show feature ranking by dataset order
########################################################################
show.rank.dataset <- function(X, ranking) {
    names(ranking) = colnames(X[,-1])
    cat('** feature ranking by dataset order \n')
    print(ranking)
}

show.rank.imp <- function(X, ranking) {
    ss <- colnames(X[,-1])[order(ranking)]
    cat('** feature ranking by importance order \n')
    print(ss)
}







