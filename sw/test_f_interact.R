########################################################################
## install packages if not exist
#install.packages("mlbench")
#install.packages("caret")
#install.packages("caret")
#install.packages("ggplot2")
#install.packages("forcats")     
#install.packages("RColorBrewer")

setwd("D:/rworks")
source("f_interact_lib.R")
library(mlbench)
library(caret)

########################################################################
## Regression


## Prepare data
data("BostonHousing")
ds <- BostonHousing[,-14]
cl <- BostonHousing[,14]    # price of a house    

## Fit a model
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 2, returnResamp="all")

model.lm <- train(x=ds,y=cl, trControl=fitControl, method="lm") 
pred <- predict(model.lm, ds)
rmse <- RMSE(pred,cl)
rmse    # rmse of a whole model

## Printf feature impact of each feature
for(i in 1:ncol(ds)) {
  cat(names(ds)[i], "\t",f.impact(model.lm, i), "\n")
}

## Degree of feature interaction between two features
result <- f.interact(model.lm, "crim", "zn")
print.fi(result)

## Create a feature interaction matrix
f.matrix <- fi.matrix(model.lm)
print(f.matrix)

## Feature interaction between single feature and others
fi.single.barplot(model.lm, "lstat")


## Grid chart of two-way feature interaction
fi.grid.chart(model.lm, "BostonHousing")

## Total degree of interaction of each feature
fi.total.plot(f.matrix) 

## Best pair of features (Postive interaction)
fi.best(f.matrix, type="pos")

## Best pair of features (Negative interaction)
fi.best(f.matrix, type="neg")








