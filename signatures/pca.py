# Potential correlations among baseflow signatures will initially be identified and summarized
# using Spearman rank correlation coefficients with the SciPy Python module and the Variance Inflation Factor (VIF)
# method using the statmodels Python module

from signatures.calculate_sigs_camels import *
import scipy
import statsmodels
import pandas
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import plotly.express as px
from statsmodels.stats.outliers_influence import variance_inflation_factor

# import data
pi_path = 'C:/Users/holta/Documents/Repositories/baseflow_prediction/signatures/camels_gw_sigs_pi.obj'
tossh_path = 'C:/Users/holta/Documents/Repositories/TOSSH'
sig_df = get_sig_df(tossh_path, pi_path)
nan_count = sig_df.isna().sum()
print(nan_count)

sig_df_drop = sig_df[["TotalRR", "EventRR", "RR_Seasonality","Recession_a_Seasonality","AverageStorage",
                      "MRC_num_segments", "BFI", "BaseflowRecessionK","First_Recession_Slope",
                      "Mid_Recession_Slope", "Spearmans_rho", "EventRR_TotalRR_ratio", "VariabilityIndex"]].dropna()

# Mid_Recession_Slope  and Variability Index have several NAS


# reduce to the five recommended signatures
sig_df_5 = sig_df[["TotalRR", "Recession_a_Seasonality", "AverageStorage", "BFI", "BaseflowRecessionK"]]
# Count the NaN values in multiple rows, don't want nas
# nan_count = sig_df_5.isna().sum()
# print(nan_count)


# PCA on signatures

# set variance to be explained??
variance = 0.9

pca = PCA(n_components=variance, svd_solver="full")
standardized_df = StandardScaler().fit_transform(sig_df_drop)

# Calculate the components
principal_components = pca.fit_transform(np.array(standardized_df))
print("Explained variance of the components (sorted in ascending order):")
print(pca.explained_variance_ratio_)

# Make it a dataframe again
principal_df = pd.DataFrame(data=principal_components, index=sig_df_drop.index)
# Give the columns more meaningful names
principal_df.columns = ["PC " + str(i) for i in range(1, len(principal_df.columns) + 1)]
print(principal_df)


# look at components during debugging (pca, components), see influence of each variable to the components

# visualize

total_var = pca.explained_variance_ratio_.sum() * 100

fig = px.scatter_3d(
    principal_components, x=0, y=1, z=2,
    title=f'Total Explained Variance: {total_var:.2f}%',
    labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
)
fig.show()



# features = ["TotalRR", "Recession_a_Seasonality", "AverageStorage", "BFI", "BaseflowRecessionK"]
# loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
#
# fig = px.scatter(principal_components, x=0, y=1)
# for i, feature in enumerate(features):
#     fig.add_annotation(
#         ax=0, ay=0,
#         axref="x", ayref="y",
#         x=loadings[i, 0],
#         y=loadings[i, 1],
#         showarrow=True,
#         arrowsize=2,
#         arrowhead=2,
#         xanchor="right",
#         yanchor="top"
#     )
#     fig.add_annotation(
#         x=loadings[i, 0],
#         y=loadings[i, 1],
#         ax=0, ay=0,
#         xanchor="center",
#         yanchor="bottom",
#         text=feature,
#         yshift=5,
#     )
# fig.show()





# # import libraries
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import spearmanr
# # set figure size
# plt.figure(figsize=(10,7))
# # Generate a mask to onlyshow the bottom triangle
# mask = np.triu(np.ones_like(sig_df.corr(), dtype=bool))
# # generate heatmap
# # ANNOTATE WITH CORRELATION SIGNIFICANCE... good have high correlation but not be signficant?
# sns.heatmap(sig_df.corr(), annot=True, mask=mask, vmin=-1, vmax=1, cmap='Blues')
# plt.title('Correlation Coefficient Of Predictors')
# plt.show()



# Variance Inflation Factor (VIF)
# # method using the statmodels Python module

# look for multicollinearity, when independent variables have a high correlation
# takes each feature and calculates regression with others
# greater VIF is greater correlation (above 5 is high)


# sig_vif = pandas.DataFrame
# # sig_vif["signature"] = sig_df.columns.tolist()
#
# # calculating VIF for each signature
# sig_vif["VIF"] = [variance_inflation_factor(sig_df.values, i) for i in range(len(sig_df.columns))]
#
# print(sig_vif)


# notes
# spatial autocorrelation of signatures and catchment characteristics?
# we are assuming the variables are indepedenent, but these variables are geographically related so they aren't
# closer catchments could have similar values, biasing the correlations

# type of variables? continuous? what value ranges? are there NAs?

# type of correlation...Pearsons or Spearmans Rank?...will be biased with spatial autocorrelation
# how do these work? linear relationships?

# PCA to create less correlated dataset (doesn't reduce the data, just weights?)
# use those groups, run clustering analysis on groups

