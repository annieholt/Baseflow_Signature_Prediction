# script to execute random forest models, predicting hydrologic signatures based on catchment attribute datasets
# Note that originally was generating random forests in Python, but decided to use R packages used in other studies

library(tidyverse)
library(randomForest)
library(caret)
library(rpart)
library(rpart.plot)
library(sf)


#### Import datasets ####

sigs_all = read.csv('E:/SDSU_GEOG/Thesis/Data/Signatures/sigs_camels_v2.csv', colClasses = c(gauge_id = "character"))

# choosing a subset, based on recommendations from McMillan et al. 2022
# (McMillan et al. 2022 found Storage Fraction was not very reliable for large samples; Average Storage more reliable)
# this dataset also includes BFI90, as recommended by Gnann et al., 2021
# no MRC_num_segments, Spearmans_rho
sigs_final = sigs_all %>% 
  select(gauge_id, TotalRR, RR_Seasonality, EventRR, Recession_a_Seasonality,
         AverageStorage, RecessionParameters_a, RecessionParameters_b, RecessionParameters_c,
         First_Recession_Slope, Mid_Recession_Slope, EventRR_TotalRR_ratio,
         VariabilityIndex, BFI, BFI_90, BaseflowRecessionK) %>% 
  as.data.frame() %>% 
  rename(RecessionParameters_T0 = RecessionParameters_c)

# camels_attribs_v2 = read.csv("E:/SDSU_GEOG/Thesis/Data/RandomForest_R/inputs/camels_attribs_v2.csv")
camels_attribs_v3 = read.csv("E:/SDSU_GEOG/Thesis/Data/RandomForest_R/inputs/camels_attribs_v3.csv", colClasses = c(gauge_id = "character"))

# camels_new_attribs = read.csv("E:/SDSU_GEOG/Thesis/Data/RandomForest_R/inputs/camels_new_attribs.csv", colClasses = c(gauge_id = "character"))

# camels_caravan_attribs_v2 = read.csv("E:/SDSU_GEOG/Thesis/Data/RandomForest_R/inputs/camels_caravan_attribs_v2.csv", colClasses = c(gauge_id = "character"))
# camels_caravan_attribs_v3 = read.csv("E:/SDSU_GEOG/Thesis/Data/RandomForest_R/inputs/camels_caravan_attribs_v3.csv", colClasses = c(gauge_id = "character"))


#### Adjust datasets for different runs ####

rf_input_attribs = camels_attribs_v3


#### Random forests, with cross validation for evaluation ####

# list of predictor variables (signatures) to loop through

# sigs_list = c('EventRR', 'TotalRR', 'RR_Seasonality', 'Recession_a_Seasonality', 'AverageStorage',
#               'RecessionParameters_a', 'RecessionParameters_b', 'RecessionParameters_T0', 
#               'First_Recession_Slope', 'Mid_Recession_Slope','EventRR_TotalRR_ratio',
#               'VariabilityIndex', 'BaseflowRecessionK',
#               'BFI', 'BFI_90')

sigs_list = c('EventRR')

rf_out_var_importance <- list()
rf_out_r2 <- list()

for(sig in sigs_list){
  
  print(sig)
  
  num_var = length(rf_input_attribs) - 1
  
  rf_df = rf_input_attribs %>% 
    left_join(sigs_final %>% select(gauge_id, sig), by = "gauge_id") %>% 
    select(-gauge_id) %>% 
    drop_na()
  
  # set seed for reproducibility
  set.seed(42)
  
  # define repeated cross-validation with 10 folds and three repeats
  # allow for parameter tuning, for mtry grid; range through the total number of predictor variables
  hyper_grid <- expand.grid(
    mtry = c(1:num_var)
  )

  tenfold_cv <- trainControl(method = 'cv', number = 10, search = "grid", verboseIter = TRUE)
  
  # train random forest using cross-validation
  
  forest <- train(
    # signature to predict
    formula(paste(sig, "~ .")),
    # input attribute dataset, includes signature
    data = rf_df,
    # Random forest method
    method = 'rf',
    # metric to evaluate model performance
    metric = 'Rsquared',
    # metric = 'RMSE',
    # Number of trees
    # adding the repeated cross validation
    trControl = tenfold_cv,
    ntree = 500,
    # hyperparameter testing
    tuneGrid = hyper_grid,
    # return importance, want %IncMSE data
    importance = TRUE
  )
  
  ## extract the trained random forest model
  rf_model <- forest$finalModel
  print(forest)
  print(rf_model)
  
  # append r2 value 
  rf_out_r2[[sig]] =mean(rf_model$rsq)
  
  # extract variable importance
  var_importance = importance(rf_model, type = 1,scale = TRUE)
  var_importance_df = as.data.frame(var_importance)
  # Add column for predictor variable names
  var_importance_df <- rownames_to_column(var_importance_df, var = "predictor")
  # Add column for response variable name
  var_importance_df$signature <- sig
  
  # append to larger output list, variable importance
  rf_out_var_importance[[sig]] <- var_importance_df
  
  }

all_var_importance <- bind_rows(rf_out_var_importance)
all_r2 = bind_rows(rf_out_r2)
