import pandas
from sklearn.inspection import permutation_importance
from sklearn.datasets import load_breast_cancer
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.model_selection import GridSearchCV

# import input data for random forest modeling
# in this case, dataset of hydrologic signatures and catchment attributes

rf_data_df = pandas.read_csv("E:/SDSU_GEOG/Thesis/Data/RandomForest/sigs_attributes_caravan_master.csv")

rf_data_df_dropna = rf_data_df.dropna()

# create a dictionary of the response variables for subsequent modeling
# keys will be hydrologic signature variable names
# include a separate predictor dataframe

sig_list = ['EventRR', 'TotalRR', 'RR_Seasonality', 'Recession_a_Seasonality', 'AverageStorage','RecessionParameters_a',
            'RecessionParameters_b', 'RecessionParameters_c', 'MRC_num_segments', 'First_Recession_Slope',
            'Mid_Recession_Slope', 'Spearmans_rho', 'EventRR_TotalRR_ratio', 'VariabilityIndex', 'BaseflowRecessionK',
            'BFI']

attrib_df = rf_data_df_dropna.drop(['gauge_id', 'geol_major_age_ma', 'non_giw_frac',
                           'EventRR', 'TotalRR', 'RR_Seasonality', 'Recession_a_Seasonality', 'AverageStorage',
                           'RecessionParameters_a', 'RecessionParameters_b', 'RecessionParameters_c',
                           'MRC_num_segments', 'First_Recession_Slope', 'Mid_Recession_Slope',
                           'Spearmans_rho', 'EventRR_TotalRR_ratio', 'VariabilityIndex', 'BaseflowRecessionK', 'BFI'],
                          axis=1)


# empty dictionary for results
sig_dic = {}

# Loop over each response variable
for var in sig_list:
    # response variable extraction
    y = rf_data_df_dropna[var]
    # add the tuple to the dictionary with the hydrologic signature name as the key
    sig_dic[var] = y

# print(sig_dic)

# random forest modeling, optimizing number of trees
# includes 10 k fold cross validation
# test/train split of 80/20

# Assuming you have a DataFrame 'predictors_data' containing predictor variables
# and a dictionary 'response_data' containing response variables for each of the 15 variables

# Define parameter grid
param_grid = {'n_estimators': [50, 100, 150, 200, 250]}

# Initialize results dictionary
rf_results = {}

# Loop over each variable
for sig in sig_dic.keys():
    X = attrib_df  # Features
    y = sig_dic[sig]  # Target variable

    # Split the data into training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Perform grid search with cross-validation
    rf = RandomForestRegressor(random_state=42)
    grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=10, scoring='r2')
    grid_search.fit(X_train, y_train)

    # Store results
    rf_results[sig] = {'best_params': grid_search.best_params_, 'best_score': grid_search.best_score_}

# Print results
for variable, result in rf_results.items():
    print(f"Variable: {variable}")
    print(f"Best parameters: {result['best_params']}")
    print(f"Best score: {result['best_score']}")
    print()

print(rf_results)