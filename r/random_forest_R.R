# script to execute random forest models, predicting hydrologic signatures based on catchment attribute datasets
# Note that originally was generating random forests in Python, but decided to use R packages used in other studies

library(tidyverse)
library(randomForest)
library(caret)
library(rpart)
library(rpart.plot)
library(sf)


#### Import datasets ####


#### Adjust datasets for different runs ####


#### Random forests, with cross validation for evaluation ####

# list of predictor variables (signatures) to loop through

sigs_list = c()

# Set seed for reproducibility
set.seed(42)

# Define repeated cross-validation with 10 folds and three repeats
# repeat_cv <- trainControl(method = 'repeatedcv', number = 10, repeats = 3, search = 'grid',
#                           allowParallel = TRUE)

repeat_cv <- trainControl(method = 'repeatedcv', number = 10, repeats = 3)

# hyperparameter testing, for mtry variable, which is the most influential
hyper_grid <- expand.grid(
  mtry = c(1:16)
)

# train random forest using cross-validation

forest <- train(
  # signature to predict
  BFI ~ .,
  # input attribute dataset, includes signature
  data = rf_caravan,
  # Random forest method
  method = 'rf',
  # metric to evaluate model performance
  metric = 'Rsquared',
  # metric = 'RMSE',
  # Number of trees
  # adding the repeated cross validation
  trControl = repeat_cv,
  ntree = 500,
  # hyperparameter testing
  # tuneGrid = hyper_grid,
  # return importance, want %IncMSE data
  importance = TRUE
)

## Extract the trained random forest model
rf_model <- forest$finalModel
print(forest)
plot(forest)


# var_importance = as.data.frame(rf_model[["importance"]])

# mse for regression
## scale... for permutation, the measures are divided by their standard errors??

var_importance = importance(rf_model, type = 1,scale = TRUE)

var_importance_df = as.data.frame(var_importance)
