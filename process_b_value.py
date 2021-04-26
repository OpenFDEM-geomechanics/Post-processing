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
import formatting_codes
import pandas as pd
import math


''' Default MATPLOTLIB Fonts '''
# plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams["figure.figsize"] = [5, 5]
matplotlib.rcParams['font.family'] = ['arial']
matplotlib.rcParams['font.size'] = 8


def load_process_file(name_csv):

    print("Processing b Values")

    # Load csv
    df_seismic = pd.read_csv(name_csv)
    x_max = df_seismic['End Step'].max()

    # Calculate and add Magnitude series
    df_seismic['Magnitude'] = (float(2) / float(3)) * (np.log10(df_seismic['Event Energy']) - 4.8)

    # Redefine colors based on Crack Type
    un_crack = df_seismic.CrackType.unique()
    print("There are %s Unique Crack Types" % len(un_crack))
    for i in un_crack:
        print("\t" + i)
    colors = {'intragranular': 'red', 'intergranular_Material': 'blue', 'intergranular': 'blue', 'intergranular_df_seismicN': 'green'}

    ## FIRST FIGURE ##
    # Plot Time vs. Event Energy colored by Failre Mode
    ax = df_seismic.plot.scatter(x='Start Step', y='Event Energy', c='Failure Mode', colormap='viridis')
    # Plot Time vs. Event Energy colored by Crack Type
    # ax = df_seismic.plot.scatter(x='Start Step', y='Event Energy', c = df_seismic['CrackType'].apply(lambda x: colors[x]))

    # Setting up the figure
    ax.set_yscale('log')
    ax.set_ylim(1)
    ax.set_xlim(0, x_max)
    ax.set_ylabel("Event Energy")
    ax.set_xlabel("Time Step of Event Occurrence")
    plt.tight_layout()
    plt.show()

    ## SECOND FIGURE ##
    # Create bins
    bins = np.arange(df_seismic['Magnitude'].min(), df_seismic['Magnitude'].max(), 0.1)
    # Draw histogram based on bins and hte Magnitude column
    ax1 = df_seismic.hist(column='Magnitude', bins=bins, grid=False, )
    # Setting up the figure
    ax1 = ax1[0]
    for x in ax1:
        x.set_ylim(0, x_max)
        x.set_xlim(-7, 2)
        x.set_title("")
        x.set_xlabel("Magnitude")
        x.set_ylabel("Frequency")
    plt.tight_layout()
    plt.show()


    ## SECOND FIGURE ##
    # Setup a secondary DataFrame comprising of the Magnitude
    stat_df_seismic = pd.DataFrame(df_seismic['Magnitude'])
    stat_df_seismic.columns = ['Magnitude']
    # Cut the dataframes into bins
    stat_df_seismic["bins"] = pd.cut(stat_df_seismic["Magnitude"], bins=bins)
    # Add a column that represents the mid value of the bin
    stat_df_seismic["bin_centres"] = stat_df_seismic["bins"].apply(lambda x: x.mid)


    # Group the DataFrame based on the bin mid-value to obtain a frequency
    statsxx_df_seismic = stat_df_seismic \
    .groupby('bin_centres') \
    ['bin_centres'] \
    .agg('count') \
    .pipe(pd.DataFrame) \
    .rename(columns = {'bin_centres': 'frequency'})
    # Reset Index
    # Add column for Frequency Cumulative Sum
    statsxx_df_seismic = statsxx_df_seismic.sort_values('bin_centres', ascending=False).reset_index()
    statsxx_df_seismic['total'] = statsxx_df_seismic['frequency'].cumsum()

    ## THIRD FIGURE ##
    fig = plt.figure()
    # Plot Scatter of Cumulative Events
    plt.scatter(statsxx_df_seismic["bin_centres"], statsxx_df_seismic["total"], marker='.', color='blue', label='Cumulative')
    # Plot Scatter of Individual Events
    plt.scatter(statsxx_df_seismic["bin_centres"], statsxx_df_seismic["frequency"], marker='.', color='orange', label='Individual')
    # Setting up the figure
    plt.yscale('log')
    plt.xlabel('Magnitude')
    plt.ylabel('Event Count')
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    try:
        folder = "/external/Size_7/post_processing"
        name_csv = os.path.join(folder, "seismicclustering.csv")
        # ss_name_csv = "/external/Ejection_Model/0002.50kPa/post_processing/history.csv"
        # name_csv = "/hdd/home/aly/Desktop/Dropbox/20170600_AA_JP_GBM_in_FDEM/20170712_MODELS/20170712_statistics/05MPa/seismicclustering.csv"
        load_process_file(name_csv)
    except KeyboardInterrupt:
        exit("TERMINATED BY USER")
