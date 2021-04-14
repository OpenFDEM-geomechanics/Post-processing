# /////////////////////////////////////////////////////////////// #
# Python Script initially created on 2021-01-07
# Compiled by Aly @ Grasselli Geomechanics Group, UofT, 2021
# Created using PyCharm
# /////////////////////////////////////////////////////////////// #

'''
GRAPHICAL REPRESENTATION.
    Inputs:
        - Post Processing directory
        - Start Time Step
        - End Time Step
        - Frequency of evaluation.
        - Name of graph
    Processing:
        - Loads seismicclistering csv
        - Manipulates data as DataFrames
    Output:
        - Time step vs. Event Energy colored based on Failure Mode
        - Histogram of Magnitudes
        - b-value curve of Cumulative and Individual event counts based on Magnitude.
'''

import sys
import os
import csv
import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import operator
# import formatting_codes
import pandas as pd
import math

from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
# import statsmodels.api as sm

''' Default MATPLOTLIB Fonts '''
# plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams["figure.figsize"] = [5, 5]
matplotlib.rcParams['font.family'] = ['arial']
matplotlib.rcParams['font.size'] = 8


def load_process_file(seismicclustering_csv, history_csv):
    print("Processing b Values")

    # Load csv
    df = pd.read_csv(seismicclustering_csv)
    stress_strain = pd.read_csv(history_csv)
    # Calculate and add Magnitude series
    df['Magnitude'] = (float(2) / float(3)) * (np.log10(df['Event Energy'] / 1.0e0) - 4.8)
    # =============================================================================================================================
    Descending_M = df.sort_values('Magnitude', ascending=False).reset_index(drop=True)

    # Redefine colors based on Crack Type
    un_crack = df.CrackType.unique()
    print("There are %s Unique Crack Types" % len(un_crack))
    for i in un_crack:
        print("\t" + i)
    colors = {'intragranular': 'red', 'intergranular_Material': 'blue', 'intergranular': 'blue',
              'intergranular_DFN': 'green'}
    # =================================Adding Axial stress and strain (and volumetric) importing to AE dataframe ======================

    stress, ax_strain, vol_strain = [], [], []
    df = df.sort_values('Start Step', ascending=True).reset_index(drop=True)
    for ts_2 in df['Start Step']:
        for ts_1 in stress_strain['Time Step (-)']:
            if ts_1 == ts_2:
                ts_indx1 = stress_strain[stress_strain['Time Step (-)'] == ts_1].index.item()
                stress.append(stress_strain['Axial stress (MPa)'][ts_indx1])
                ax_strain.append(stress_strain['Strain y from platen (%)'][ts_indx1])
                vol_strain.append(stress_strain['Volumetric strain (%)'][ts_indx1])
    df.loc[:, 'Axial Stress'] = stress
    df.loc[:, 'Axial Strain'] = ax_strain
    df.loc[:, 'Volumetric Strain'] = vol_strain

    # =============================================================================================================================
    ## FIRST FIGURE ##
    # Plot Time vs. Event Energy colored by Failre Mode
    ax = df.plot.scatter(x='Start Step', y='Event Energy', c='Failure Mode', colormap='viridis')
    # Plot Time vs. Event Energy colored by Crack Type

    # Setting up the figure
    ax.set_yscale('log')
    ax.set_ylim(1)
    # ax.set_xlim(0, 4000)
    ax.set_ylabel("Event Energy")
    ax.set_xlabel("Time Step of Event Occurrence")
    plt.tight_layout()
    plt.show()

    ## SECOND FIGURE ##
    # Create bins
    bin_size = 0.1
    bins = np.arange(df['Magnitude'].min() - bin_size, df['Magnitude'].max() + bin_size, bin_size)
    # Draw histogram based on bins and hte Magnitude column
    ax1 = df.hist(column='Magnitude', bins=bins, grid=False, )
    # Setting up the figure
    ax1 = ax1[0]
    for x in ax1:
        x.set_ylim(0, 350)
        # x.set_xlim(-7, 2)
        x.set_title("")
        x.set_xlabel("Magnitude")
        x.set_ylabel("Frequency")
    plt.tight_layout()
    plt.show()

    ##

    ## SECOND FIGURE ##
    # Setup a secondary DataFarme compirsing of the Magnitude
    stat_df = pd.DataFrame(df['Magnitude'])
    stat_df.columns = ['Magnitude']

    # Cut the dataframes into bins
    stat_df["bins"] = pd.cut(stat_df["Magnitude"], bins=bins)
    # Add a column that represents the mid value of the bin
    stat_df["bin_centres"] = stat_df["bins"].apply(lambda x: x.mid)

    ##

    # Group the DataFrame based on the bin mid-value to obtain a frequency
    statsxx_df = stat_df \
        .groupby('bin_centres') \
        ['bin_centres'] \
        .agg('count') \
        .pipe(pd.DataFrame) \
        .rename(columns={'bin_centres': 'frequency'})
    # Reset Index to calculate the cumulative from the highest magnitude
    statsxx_df = statsxx_df.sort_values('bin_centres', ascending=False).reset_index()

    # Add column for Frequency Cumulative Sum
    statsxx_df['total'] = statsxx_df['frequency'].cumsum()

    ## THIRD FIGURE ##
    fig = plt.figure()
    # Plot Scatter of Cumulative Events
    plt.scatter(statsxx_df["bin_centres"], statsxx_df['total'], marker='.', color='blue',
                label='Cumulative count (binned)')
    plt.scatter(Descending_M['Magnitude'], Descending_M['Magnitude'].index + 1, marker='.', color='green',
                label='Cumulative count (all AE events)')
    # Plot Scatter of Individual Events
    plt.scatter(statsxx_df["bin_centres"], statsxx_df["frequency"], marker='.', color='red',
                label='Non-cumulative count')
    # Setting up the figure

    # ================================================= b-value calculation using Maximum Likelihood Estimation Method (MLE)  ===================================================
    # Mc = -1.35
    Mc_indx = statsxx_df[statsxx_df['frequency'] == statsxx_df['frequency'].max()].index.item()
    print(Mc_indx)
    Mc = statsxx_df['bin_centres'][Mc_indx]

    Mbin = bin_size
    total_events = []
    Mlarger = []
    # for i in Descending_M["Magnitude"]:

    for i in statsxx_df["bin_centres"]:
        if i >= Mc:
            indx = statsxx_df[statsxx_df["bin_centres"] == i].index.item()
            events_in_bin = statsxx_df['frequency'][indx]
            total_events.append(events_in_bin)
            total_in_bin = i * events_in_bin  # sum of event magnitudes in a single bin
            Mlarger.append(total_in_bin)

    Mavg = sum(Mlarger) / sum(total_events)
    b_value = np.log10(math.e) / (Mavg - (Mc - 0.5 * Mbin))
    b_value_uncertainty = b_value / len(total_events)  # Aki (1965) relation is used

    print("b-value using MLE method is:", b_value)
    print("Mc is:", Mc)
    print("Uncertainty associated with b-value:", b_value_uncertainty)

    Cum_F_Mc = float(statsxx_df['total'][Mc_indx])

    C = np.log10(Cum_F_Mc) + b_value * Mc
    M_bins_list = np.array(statsxx_df['bin_centres'].tolist())
    Y = []
    for i in M_bins_list:
        yi = pow(10, -i * b_value + C)
        Y.append(yi)

    plt.plot(M_bins_list, Y, color='orange', label="b-value (MLE)")
    plt.scatter(Mc, Cum_F_Mc, color="purple", s=50, marker="s", label="Magnitude of Completeness (Mc)")

    # ================================================ b-value calculation using least squares regression method LSM ============================
    Mags = []
    Freqs = []
    for i in statsxx_df["bin_centres"]:
        if i >= Mc:
            Mags.append(i)
            indx = statsxx_df[statsxx_df["bin_centres"] == i].index.item()
            Freqs.append(statsxx_df['total'][indx])

    Mags = np.array(Mags)
    Freqs = np.array(Freqs)

    Mags = np.array(Mags).reshape((len(Mags), 1))

    y = np.log10(Freqs)  # Apply natural log function to the target
    model = LinearRegression()
    model.fit(Mags, y)
    y_pred = np.power(10, model.predict(Mags))
    plt.plot(Mags, y_pred, 'k--', label="b-value (LSM) R^2 = %1.4f" % model.score(Mags, y))
    print('R^2', model.score(Mags, y))
    print('b-value using LSM is:', abs(model.coef_))

    # ========================================================================================================================================

    plt.yscale('log')
    plt.ylim([0.8, 10000])
    # plt.autoscale(enable=True, axis='y')
    plt.xlabel('Magnitude')
    plt.ylabel('Event Count')
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ========================================== Stress and strain vs. AE ======================================
    # ============ Now binning is based on time steps not magnitude ==================
    bin_size_2 = 50
    bins_2 = np.arange(df['Start Step'].min() - bin_size_2, df['Start Step'].max() + bin_size_2, bin_size_2)

    # Setup a secondary DataFarme compirsing of the timesteps
    ts_df = pd.DataFrame(df['Start Step'])
    ts_df.columns = ['Start Step']

    # Cut the dataframes into bins
    ts_df["bins"] = pd.cut(ts_df["Start Step"], bins=bins_2)
    # Add a column that represents the mid value of the bin
    ts_df["bin_centres"] = ts_df["bins"].apply(lambda x: x.mid)
    stresses_in_bin = []
    print(ts_df)
    '''
    for i in df['Start Step']:
        for ii in ts_df["bins"]:
            if i in ii:
                stress_indx = df[df['Start Step'] in ii].index.item()
                stresses_in_bin.append(df['Axial Stress'][stress_indx])
                print(stresses_in_bin)
    '''


Equiv_D = [1.65, 2.55, 3.1, 12]
b_values = [1.8048242882381038, 1.1999387767041843, 1.2005900452786096, 0.771860483503401]
errors = [0.11280151801488149, 0.014457093695231136, 0.0521995671860265, 0.024898725274303257]
plt.errorbar(Equiv_D, b_values, yerr=errors, fmt='_', color='red')
plt.xlabel('Mean grain diameter, mm')
plt.ylabel('b-value and associated error')
plt.show()

if __name__ == "__main__":
    try:
        # seismicclustering_csv = "E:\Google Drive\PhD\Projects\CIV1429 Project\AE_postprocessing\Cleavage_41.csv"
        # history_csv = "E:\Google Drive\PhD\Projects\CIV1429 Project\AE_postprocessing\Tugrul_3_history.csv"

        # name_csv = "/external/ARMA_JP/00MPa/000/post_processing/seismicclustering.csv"
        seismicclustering_csv = "/external/Ejection_Model/0002.50kPa/post_processing/seismicclustering.csv"
        history_csv = "/external/Ejection_Model/0002.50kPa/post_processing/history.csv"
        load_process_file(seismicclustering_csv, history_csv)
    except KeyboardInterrupt:
        exit("TERMINATED BY USER")

