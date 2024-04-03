# Script to review CAMELS dataset and random forest modeling

library(tidyverse)
library(randomForest)
library(caret)
library(rpart)
library(rpart.plot)
library(sf)


#### Hydrologic Signatures ####

sigs_c = read.csv('E:/SDSU_GEOG/Thesis/Data/Signatures/sigs_camels.csv', colClasses = c(gauge_id = "character"))

sigs_c_2 = sigs_c %>% 
  select(gauge_id, TotalRR, RR_Seasonality, EventRR, Recession_a_Seasonality,
         AverageStorage, RecessionParameters_a, RecessionParameters_b, RecessionParameters_c, MRC_num_segments,
         First_Recession_Slope, Mid_Recession_Slope, Spearmans_rho, EventRR_TotalRR_ratio,
         VariabilityIndex, BFI, BaseflowRecessionK) %>% 
  as.data.frame()


#### New Attributes ####

nwi_c = st_read('E:/SDSU_GEOG/Thesis/Data/NWI_outputs/Shapefiles/nwi_camels_metrics_ecoregions.shp') %>% 
  as.data.frame() %>% 
  select(gauge_id, shed_area, fresh, lake, other)

giws = st_read('E:/SDSU_GEOG/Thesis/Data/GIWs/giws_metrics.shp') %>% 
  as.data.frame() %>% 
  select(gauge_id, area_frac)

# age by major geologic unit
geol_c = st_read('E:/SDSU_GEOG/Thesis/Data/Geology_outputs/Shapefiles/sgmc_camels_metrics.shp') %>% 
  as.data.frame() %>% 
  # select(gauge_id, av_age) %>% 
  rename(geol_major_age_ma = av_age)

# average age of catchment (age weighted by area of geologic unit type)
geol_c_av = st_read('E:/SDSU_GEOG/Thesis/Data/Geology_outputs/Shapefiles/sgmc_camels_metrics_age_weighted.shp') %>% 
  as.data.frame() %>% 
  # select(gauge_id, av_age_w) %>% 
  rename(geol_av_age_ma = av_age_w)


new_attrib = nwi_c %>% 
  left_join(giws, by = "gauge_id") %>% 
  # some places have zero isolated wetlands, so replace NAs with zero since shapefiles weren't returned for those
  mutate(giw_frac= ifelse(is.na(area_frac), 0, area_frac)) %>% 
  mutate(non_giw_frac = fresh + lake - giw_frac) %>% 
  left_join(geol_c, by = "gauge_id") %>% 
  left_join(geol_c_av, by = "gauge_id") %>% 
  select(gauge_id, giw_frac, non_giw_frac, geol_major_age_ma, geol_av_age_ma)


#### CAMELS attributes ####

camels_path = "E:/SDSU_GEOG/Thesis/Data/CAMELS/camels-20230412T1401Z"

# camels_name = read.table(paste(camels_path, "/camels_name.txt", sep = ""), sep = ";", header = TRUE)
camels_hydro = read.table(paste(camels_path, "/camels_hydro.txt", sep = ""), sep = ";", header = TRUE)

# just area, mean slope, mean elevation
camels_topo = read.table(paste(camels_path, "/camels_topo.txt", sep = ""), sep = ";", header = TRUE)

camels_clim = read.table(paste(camels_path, "/camels_clim.txt", sep = ""), sep = ";", header = TRUE)

# remove second dominant geological class, as many catchments are one major geology
# follwoing Addor et al., 2018
camels_geol = read.table(paste(camels_path, "/camels_geol.txt", sep = ""), sep = ";", header = TRUE)

# remove soil porosity and conductivity, following Addor et al., 2018
# both are highly correlated with sand fraction
camels_soil = read.table(paste(camels_path, "/camels_soil.txt", sep = ""), sep = ";", header = TRUE)

# remove lai_diff and gvf_diff, following Addor et al., 2018
# both are highly correlated with leaf area index maximum
camels_vege = read.table(paste(camels_path, "/camels_vege.txt", sep = ""), sep = ";", header = TRUE)


# join all together, make sure gauge_id has leading zeros
camels_attribs = camels_topo %>% 
  select(-gauge_lat, -gauge_lon, -area_geospa_fabric) %>% 
  left_join(camels_clim, by = "gauge_id") %>% 
  left_join(camels_geol %>% select(-geol_2nd_class, -glim_2nd_class_frac), by = "gauge_id") %>% 
  left_join(camels_soil %>% select(-soil_porosity, -soil_conductivity), by = "gauge_id") %>% 
  left_join(camels_vege %>% select(-lai_diff, -gvf_diff), by = "gauge_id") %>% 
  mutate(gauge_id = str_pad(string = as.numeric(gauge_id), width = 8, side = 'left', pad = 0))


#### CARAVAN attributes ####

camels_caravan = read_csv("E:/SDSU_GEOG/Thesis/Data/Caravan/attributes/camels/attributes_caravan_camels.csv")

camels_hydroatlas = read_csv("E:/SDSU_GEOG/Thesis/Data/Caravan/attributes/camels/attributes_hydroatlas_camels.csv")

camels_caravan_other = read_csv("E:/SDSU_GEOG/Thesis/Data/Caravan/attributes/camels/attributes_other_camels.csv")


# add catchment area???? look at the other attribute table

# seleted attributes: 
# topography: mean elevation, average terrain slope, average stream gradient
# landcover: forest cover extent, glacial extent (there are also land use classes, but not including for now)
# soils: clay/silt/sand fraction in soils, organic carbon contect, karst area, soil erosion
# geology: NA (just have lith classes, which are a lot of categories; add perm/porosity??)

# climate: all caravan provides (including frac_snow!!)

camels_caravan_attribs = camels_hydroatlas %>% 
  select(gauge_id, ele_mt_smn, slp_dg_sav, sgr_dk_sav,
         for_pc_sse, gla_pc_sse,
         cly_pc_sav, slt_pc_sav, snd_pc_sav, soc_th_sav, kar_pc_sse, ero_kh_sav) %>% 
  left_join(camels_caravan, by = "gauge_id") %>% 
  left_join(camels_caravan_other %>% select(gauge_id, area), by = "gauge_id") %>% 
  mutate(gauge_id = gsub("^camels_", "", gauge_id))


camels_caravan_attribs_v2 = camels_hydroatlas %>% 
  select(gauge_id, ele_mt_smn, slp_dg_sav,
         for_pc_sse,
         cly_pc_sav, slt_pc_sav, snd_pc_sav, soc_th_sav, kar_pc_sse, swc_pc_syr) %>% 
  left_join(camels_caravan, by = "gauge_id") %>% 
  left_join(camels_caravan_other %>% select(gauge_id, area), by = "gauge_id") %>% 
  mutate(gauge_id = gsub("^camels_", "", gauge_id))

#### data for practicing RF modeling ####

# ## Get the iris dataset
# data("iris")
# ## View first few rows
# head(iris)


rf_df_final = sigs_c_2 %>% 
  left_join(new_attrib, by = "gauge_id") %>% 
  left_join(camels_attribs, by = "gauge_id")

# write.csv(rf_df_final, "E:/SDSU_GEOG/Thesis/Data/RandomForest/sigs_attributes_master.csv", row.names = FALSE)


rf_caravan_final = camels_caravan_attribs_v2 %>% 
  left_join(new_attrib, by = "gauge_id") %>% 
  left_join(sigs_c_2, by = "gauge_id") %>% 
  mutate(total_wet = giw_frac + non_giw_frac) %>% 
  filter(total_wet > 0.05)

# write.csv(rf_caravan_final, "E:/SDSU_GEOG/Thesis/Data/RandomForest/outputs/sigs_attributes_caravan_master_v2_wet.csv", row.names = FALSE)


rf_df_1 = camels_hydro %>% 
  select(gauge_id, baseflow_index) %>% 
  left_join(camels_attribs, by = "gauge_id")
  # drop_na() %>% 
  # select(-gauge_id) %>% 
  # select(-low_prec_timing, -high_prec_timing, -geol_1st_class, -glim_1st_class_frac, -dom_land_cover, -dom_land_cover_frac)



# write.csv(rf_df_1, "E:/SDSU_GEOG/Thesis/Data/RandomForest/BFI_attributes_test_data.csv")

rf_df = camels_hydro %>% 
  select(gauge_id, baseflow_index) %>% 
  left_join(camels_attribs, by = "gauge_id") %>% 
  drop_na() %>%
  select(-gauge_id) %>% 
  select(-low_prec_timing, -high_prec_timing, -geol_1st_class, -glim_1st_class_frac, -dom_land_cover, -dom_land_cover_frac)


# select(-dom_land_cover, -geol_1st_class)
  
#### RF AND PERMUTATION ####

X <- rf_df %>% select(-baseflow_index)
y <- as.factor(rf_df$baseflow_index)

# Set seed for reproducibility
set.seed(42)

# Split the data into training and test sets
set.seed(42)  
indices <- sample(1:nrow(X), size = 0.75 * nrow(X))
X_train <- X[indices, ]
X_test <- X[-indices, ]

# Convert y to character, subset, and convert back to factor
y_char <- as.character(y)
y_train <- as.factor(y_char[indices])
y_test <- as.factor(y_char[-indices])

# Create and train the Random Forest classifier
clf <- randomForest(x = X_train, y = y_train, ntree = 100, seed = 42)

# Predict on the test set
predictions <- predict(clf, newdata = X_test)

# Create a confusion matrix
conf_matrix <- confusionMatrix(predictions, y_test)

# Calculate accuracy on the test data
accuracy <- conf_matrix$overall["Accuracy"]

# Print the accuracy on the test data
print(paste("Baseline accuracy on test data:", round(accuracy, 2)))









#### CROSS-FOLD VALIDATION ####

## Set seed for reproducibility
set.seed(123)

## Define repeated cross validation with 5 folds and three repeats
repeat_cv <- trainControl(method='repeatedcv', number=5, repeats=3)

## Split the data so that we use 70% of it for training
# train_index <- createDataPartition(y=iris$Species, p=0.7, list=FALSE)
train_index <- createDataPartition(y=rf_df$baseflow_index, p=0.7, list=FALSE)

## Subset the data
training_set <- rf_df[train_index, ]
testing_set <- rf_df[-train_index, ]

## Set seed for reproducibility
set.seed(123)

## Train a random forest model
forest <- train(
  
  # Formula. We are using all variables to predict Species
  baseflow_index~., 
  
  # Source of data; remove the Species variable
  data=training_set, 
  
  # `rf` method for random forest
  method='rf', 
  
  # Add repeated cross validation as trControl
  trControl=repeat_cv,
  
  # Accuracy to measure the performance of the model
  metric='MAE')

## Print out the details about the model
forest$finalModel



## Get variable importance, and turn into a data frame
var_imp <- varImp(forest, scale=FALSE)$importance
var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)

## Create a plot of variable importance
var_imp %>%
  
  ## Sort the data by importance
  arrange(importance) %>%
  
  ## Create a ggplot object for aesthetic
  ggplot(aes(x=reorder(variables, importance), y=importance)) + 
  
  ## Plot the bar graph
  geom_bar(stat='identity') + 
  
  ## Flip the graph to make a horizontal bar plot
  coord_flip() + 
  
  ## Add x-axis label
  xlab('Variables') +
  
  ## Add a title
  labs(title='Random forest variable importance') + 
  
  ## Some layout for the plot
  theme_minimal() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 20), 
  )

# ## Generate predictions
# y_hats <- predict(
#   
#   ## Random forest object
#   object=forest, 
#   
#   ## Data to use for predictions; remove the Species
#   newdata=testing_set[, -1])
# 
# 
# ## Print the accuracy
# accuracy <- mean(y_hats == testing_set$baseflow_index)*100
# cat('Accuracy on testing data: ', round(accuracy, 2), '%',  sep='')
