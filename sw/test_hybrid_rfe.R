################################################################################################
## Test hybrid-RFE 
################################################################################################
setwd('C:\\Users\\mango\\OneDrive - 단국대학교\\논문지도\\rfe_feature_selection\\github_publish')

source('hybrid_RFE.R')

## dataset : iris ###################################################################

ds = iris[,c(5,1,2,3,4)]
head(ds)

ranking = hybridRFE.kfold(ds,k=5,halve.above=500,max.row=1000, method='s.sum') 

show.rank.dataset(ds, ranking)

show.rank.imp(ds, ranking) 

## dataset : Ionosphere (from mlbench) ###############################################
library(mlbench)
data(Sonar)
ds = Sonar[,c(61,1:60)]
head(ds)

ranking = hybridRFE.kfold(ds,k=5,halve.above=500,max.row=1000, method='s.sum') 

show.rank.dataset(ds, ranking)

show.rank.imp(ds, ranking) 

# SVM Test & graph
library(e1071)
features <- c(5,10,15,20,25,30,35,40,45,50,55,60)

X <- ds[,-1]
Y <- ds[,1]

acc = c()
i <- 1
for (no in features) {
  model <- svm(X[,1:no],Y)
  pred <- predict(model,X[,1:no])
  acc[i] <- mean(Y==pred) 
  i=i+1
}

plot(features, acc, main='Training Accuracy', xlab='# of features', ylab='accuracy',
     type='b', axes=T)


