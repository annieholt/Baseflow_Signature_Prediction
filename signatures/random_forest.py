import pandas
from sklearn.inspection import permutation_importance
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

# import input data for random forest modeling
# in this case, dataset of hydrologic signatures and catchment attributes

rf_data_df = pandas.read_csv("E:/SDSU_GEOG/Thesis/Data/RandomForest/sigs_attributes_caravan_master.csv")

rf_data_df_dropna = rf_data_df.dropna()

# rf_data_wet = dataframe[dataframe['Percentage'] > 70]

# create a dictionary of the response variables for subsequent modeling
# keys will be hydrologic signature variable names
# include a separate predictor dataframe

sig_list = ['EventRR', 'TotalRR', 'RR_Seasonality', 'Recession_a_Seasonality', 'AverageStorage','RecessionParameters_a',
            'RecessionParameters_b', 'RecessionParameters_c', 'MRC_num_segments', 'First_Recession_Slope',
            'Mid_Recession_Slope', 'Spearmans_rho', 'EventRR_TotalRR_ratio', 'VariabilityIndex', 'BaseflowRecessionK',
            'BFI']

# sig_list = ['EventRR', 'TotalRR']

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
# param_grid = {'n_estimators': [10, 50, 100, 150, 200, 250, 500, 1000]}
# 100 trees was generally enough to maximize accuracy but minimze computation
# (accuracy only varied about 1% with up to 1000 trees)
# accuracy varied a bit more with max_features... so letting log2 or sqrt vary for each signature
param_grid = {'n_estimators': [100], 'max_features': ['log2', 'sqrt']}

# Initialize results dictionary
rf_results = {}

# Loop over each variable
for sig in sig_dic.keys():
    X = attrib_df  # Features
    y = sig_dic[sig]  # Target variable

    # Split the data into training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Perform grid search with cross-validation
    rf = RandomForestRegressor(random_state=42, n_estimators=100)
    grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=10, scoring='r2')
    grid_search.fit(X_train, y_train)

    # Store results
    rf_results[sig] = {'best_params': grid_search.best_params_, 'best_score': grid_search.best_score_,
                       'best_estimator': grid_search.best_estimator_,
                       'mean_test_score':  grid_search.cv_results_['mean_test_score']}


rf_final_results = {}

# Based on best hyperparameters, generate final random forest models for each signature
# Loop over each variable
for sig in sig_dic.keys():
    X = attrib_df  # Features
    y = sig_dic[sig]  # Target variable

    # Split the data into training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    max_features = rf_results[sig]['best_params']['max_features']
    n_estimators = rf_results[sig]['best_params']['n_estimators']

    rf = RandomForestRegressor(n_estimators=n_estimators, max_features=max_features, random_state=42)
    rf.fit(X_train, y_train)
    # print(f"R-squared on test data: {rf.score(X_test, y_test):.2f}")c

    # Perform 10-fold cross-validation
    cv_scores = cross_val_score(rf, X, y, cv=KFold(n_splits=10, shuffle=True, random_state=42))

    # mean squared error
    y_pred = rf.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    # print("Mean Squared Error (MSE):", mse)

    # store results
    rf_final_results[sig] = {'mean R-squared from cross-validation': cv_scores.mean(),
                             'R-squared on test data': rf.score(X_test, y_test),
                             'mean squared error': mse}

    # permutation importance
    # NOTE: this is different then increase in mean squared error
    # so later, try MSE... and decide if perturb the variable or remove it entirely
    
    # on test data (prob also need to check on train data to test for overfitting??)
    var_importance = permutation_importance(rf, X_test, y_test, n_repeats=10, random_state=42)

    sorted_importances_idx = var_importance.importances_mean.argsort()
    importances = pandas.DataFrame(
        var_importance.importances[sorted_importances_idx].T,
        columns=X.columns[sorted_importances_idx],
    )
    ax = importances.plot.box(vert=False, whis=10)
    ax.set_title("Permutation Importances (test set)")
    ax.axvline(x=0, color="k", linestyle="--")
    ax.set_xlabel("Decrease in accuracy score")
    ax.figure.tight_layout()
    plt.show()


print(rf_results)
print(rf_final_results)


# Extract response variables and corresponding cross-validated accuracy scores from the results dictionary
response_variables = list(rf_results.keys())


# Create a separate plot for each response variable
for sig in response_variables:
    # Extract the cross-validated accuracy scores and corresponding hyperparameter combinations
    # for the current response variable
    scores = rf_results[sig]['mean_test_score']
    param_values = rf_results[sig]['best_params'].values()
    # param_combinations = [str(params) for params in zip(*param_values)]

    # Create a plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(range(len(scores)), scores, marker='o', linestyle='-')
    ax.set_xticks(range(len(scores)))
    # ax.set_xticklabels(param_combinations, rotation=45, ha='right')
    ax.set_xlabel('Parameter Combination')
    ax.set_ylabel('Accuracy Score')
    ax.set_title(f'Accuracy Scores for Response Variable {sig}')
    plt.tight_layout()
    plt.show()

print(rf_results)