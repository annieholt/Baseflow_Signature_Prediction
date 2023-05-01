# Workflow to generate groundwater signatures for CAMELS catchments
# runs MATLAB TOSSH Toolbox, then reformats output into one final dataframe
# Dataframe contains results for all signatures, and gauge location information

# Note that this workflow was run like Addor et al., 2018:
# time period 1 Oct 1989 to 30 Sep 2009; 671 CAMELS catchments

# Import packages
import matlab.engine
import pandas as pd
import pickle
import os


# Function to run MATLAB codes from the TOSSH toolbox in Python
def run_tossh(tossh_path):
    # create MATLAB engine
    eng = matlab.engine.start_matlab()
    # add paths to necessary MATLAB scripts
    # scripts were forked from GitHub Gnann et al., 2021
    # several scripts were updated for personal use... basically needed to feed CAMELS data into groundwater workflow
    path = eng.genpath(tossh_path)
    eng.addpath(path, nargout=0)
    # camels_data = eng.loadCAMELSstruct()

    # results_daily = eng.CAMELS_dataprep()
    # datasets for Q, t, P, PET
    # Q_mat: streamflow [mm/timestep] matrix (cell array)
    # t_mat: time [Matlab datenum] matrix (cell array)
    # P_mat: precipitation [mm/timestep] matrix (cell array)

    # calculate groundwater signature set
    # this takes about a minute
    output = eng.CAMELS_groundwater_2()

    return output

# function to load/generate data, and then extract into dataframes for analyses
# requires paths, for pickle file storing outputs and for tossh toolbox repo

# returns dataframe, for now just 5 recommended signatures, can later do all 15??
def get_sig_df(tossh_path, pi_path):

    # create empty object for tossh results
    tossh_results = None

    # if file doesn't exist, run tossh and assign to results
    # otherwise import pickle file and assign to results

    if not os.path.exists(pi_path):
        # run tossh function, which calls matlab tools
        tossh_results = run_tossh(tossh_path)
        # then save into pickle
        file_pi = open(pi_path, 'wb')
        pickle.dump(tossh_results, file_pi)
    else:
        # load pickle file
        file_pi = open(pi_path, 'rb')
        tossh_results = pickle.load(file_pi)

    # now put results into dataframes
    id_df = pd.DataFrame(tossh_results['gauge_id'])
    id_df.rename(columns={0: 'gauge_id'}, inplace=True)
    id_df.reset_index(inplace=True)

    lat_df = pd.DataFrame(tossh_results['gauge_lat'])
    lat_df.rename(columns={0: 'gauge_lat'}, inplace=True)
    lat_df.reset_index(inplace=True)

    lon_df = pd.DataFrame(tossh_results['gauge_lon'])
    lon_df.rename(columns={0: 'gauge_lon'}, inplace=True)
    lon_df.reset_index(inplace=True)

    sig_df = pd.DataFrame(tossh_results['sigs'])
    sig_df.reset_index(inplace=True)

    # join all dataframes into one final dataframe, based on index/gage number
    sig_df_final = id_df.merge(lat_df, on='index').merge(lon_df, on='index').merge(sig_df, on='index')
    # print(sig_df_final.dtypes)

    # drop decimals, make gauge id a string, add leading zero so always 8 digits long
    sig_df_final['gauge_id'] = sig_df_final['gauge_id'].astype(int).astype(str).str.zfill(8)
    print(sig_df_final.dtypes)

    # no StorageFraction, RecessionParameters for now, didn't deal with the formatting
    sig_df_2 = sig_df_final[["gauge_id", "gauge_lat", "gauge_lon", "TotalRR", "EventRR", "RR_Seasonality",
                             "Recession_a_Seasonality",
                             "AverageStorage", "MRC_num_segments", "BFI", "BaseflowRecessionK",
                             "First_Recession_Slope", "Mid_Recession_Slope", "Spearmans_rho", "EventRR_TotalRR_ratio",
                             "VariabilityIndex"]]

    # select 5 most recommended signatures by McMillan et al., 2022
    # Total RR, Average Storage, Recession a seasonality, BFI, Baseflow Recession K
    # sig_df_2 = sig_df_final[["TotalRR", "Recession_a_Seasonality", "AverageStorage", "BFI", "BaseflowRecessionK"]]
    # print(sig_df_2)

    return sig_df_2




# note that some sigs have multiple output values:
# Storage fraction: S_fraction: ratio between active and total storage capacity [-]
# %   S_active: active storage capacity [mm]
# %   S_total: total storage capacity [mm]

# Recession Parameters
# matrix with parameters alpha, beta (=1 for
# %       exponential fit) for each recession segment
