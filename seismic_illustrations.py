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
matplotlib.rcParams['font.size'] = 12


def load_process_file(folder, binnings=10, big_bins=4):
    name_csv = os.path.join(folder, "seismicclustering.csv")
    ss_name_csv = os.path.join(folder, "history.csv")


    # Load seismicclustering file
    df_seismic = pd.read_csv(name_csv)

    # Load history file
    df_stress_strain = pd.read_csv(ss_name_csv)

    # Calculate and add Magnitude series
    df_seismic['Magnitude'] = (float(2) / float(3)) * (np.log10(df_seismic['Event Energy']) - 4.8)

    # Concatenate the crack count into the stressstrain DataFrame
    seismic_count = df_seismic['Start Step'].value_counts().sort_index(ascending=True)
    seismic_count.name = "SeismicCluster"
    df_combined = pd.concat([df_stress_strain, seismic_count],  axis=1, join_axes=[df_stress_strain.index])
    df_combined['SeismicCluster'].replace(np.nan, 0, inplace=True)

    # Assign CrackMode
    # create a list of our conditions
    conditions = [
        (df_seismic['Failure Mode'] >= 1) & (df_seismic['Failure Mode'] < 1.5),
        (df_seismic['Failure Mode'] > 1.5) & (df_seismic['Failure Mode'] <= 2.0),
        (df_seismic['Failure Mode'] == 3.0),
        (df_seismic['Failure Mode'] == 1.5)
    ]
    # create a list of the values we want to assign for each condition
    values = ['Tensile Dominant', 'Shear Dominant', 'Mixed Mode', 'Mixed Mode']

    # create a new column and use np.select to assign values to it using our lists as arguments
    df_seismic['CrackMode'] = np.select(conditions, values)

    # Assign CrackMode
    # create a list of our conditions
    conditions_type = [
        (df_seismic['CrackType'] == 'intergranular_Material'),
        (df_seismic['CrackType'] == 'intergranular_DFN'),
        (df_seismic['CrackType'] == 'intragranular')
    ]
    # create a list of the values we want to assign for each condition
    values_type = ['Intergranular', 'Intergranular','Intragranular']

    # create a new column and use np.select to assign values to it using our lists as arguments
    df_seismic['CrackType'] = np.select(conditions_type, values_type)


    # Get max time step
    x_max = df_seismic['End Step'].max()


    '''
    Plot stress-strain vs. Seismic Data
    '''

    fig1 = plt.figure(figsize=[12.8, 7.2])
    ax1 = fig1.add_subplot(1, 1, 1)

    ax1.plot(df_combined['Strain y from platen (%)'], df_combined['Axial stress (MPa)'], linestyle="-")

    ax1.set_xlabel("Strain (%)")
    ax1.set_xlim(xmin=0)
    ax1.set_ylabel("Stress (MPa)")
    ax1.set_ylim(ymin=0)
    #
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel('Cumulative AE Event Count', color='red')
    ax2.plot(df_combined['Strain y from platen (%)'], df_combined['SeismicCluster'].cumsum(), linestyle="-", color='red')
    ax2.set_ylim(ymin=0)

    ind_count = np.diff(df_combined['SeismicCluster'].cumsum())
    ind_count = [0] + ind_count.tolist()

    b= [sum(df_combined['Strain y from platen (%)'][i:i+binnings])/binnings for i in range(0, int(x_max), binnings)]
    bx = [sum(ind_count[i:i+binnings])/binnings for i in range(0, int(x_max), binnings)]
    width = (b[0] - b[1]) / 2

    ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    # right, left, top, bottom

    ax3.set_ylim(ymin=0)
    ax3.set_ylabel('Individual AE Event Count', color='green')
    ax3.bar(b, bx, color='green', width=width, alpha=0.5)
    ax3.set_ylim(ymin=0, ymax=max(bx) * 1.1)
    ax3.spines['right'].set_position(('outward', 45))

    plt.tight_layout()
    plt.show()


    '''
    Plot stress-strain vs. Crack Patterns
    '''

    filter_on = ['Group Type', 'CrackType', 'CrackMode']
    # filter_on = ['CrackMode']
    for filter_criteria in filter_on:
        fig2 = plt.figure(figsize=[12.8, 7.2])
        ax1 = fig2.add_subplot(1, 1, 1)

        ax1.plot(df_combined['Strain y from platen (%)'], df_combined['Axial stress (MPa)'], linestyle="-")

        ax1.set_xlabel("Strain (%)")
        ax1.set_xlim(xmin=0)
        ax1.set_ylabel("Stress (MPa)")
        ax1.set_ylim(ymin=0)

        un_crack = df_seismic[filter_criteria].unique()
        ax4 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        ax4.set_ylabel('Cumulative Crack Count', color='black')
        print(df_seismic[filter_criteria].value_counts())
        for i in un_crack:
            crack_filter = df_seismic[df_seismic[filter_criteria]==i]
            crack_count = crack_filter['Start Step'].value_counts().sort_index(ascending=True)
            crack_count.name = i
            df_combined = pd.concat([df_combined, crack_count], axis=1, join_axes=[df_combined.index])
            df_combined[i].replace(np.nan, 0, inplace=True)

        for idx, var in enumerate(un_crack):
            by = [df_combined[var].cumsum()[i:i + binnings].iloc[-1] for i in range(0, int(x_max), binnings)]


            if idx == 0:
                # print(var, by)
                ax4.bar(b, by, width=width, label=var, alpha=0.5)
                ll = by
            else:

                ax4.bar(b, by, bottom=ll, width=width, label=var, alpha=0.5)
                ll = [x + y for x, y in zip(ll, by)]

        ax4.set_ylim(ymin=0)
        ax4.legend(title=filter_criteria)
        plt.tight_layout()
        plt.savefig(os.path.join(folder, 'stress_strain_' + filter_criteria + '.pdf'))
        plt.savefig(os.path.join(folder, 'stress_strain_' + filter_criteria + '.svg'))
        plt.savefig(os.path.join(folder, 'stress_strain_' + filter_criteria + '.png'))

        plt.show()


    bb = [sum(df_combined['Time Step (-)'][i:i+big_bins])/big_bins for i in range(0, int(x_max), big_bins)]
    width = (bb[0] - bb[1]) / 2

    fig11 = plt.figure(figsize=[12.8, 7.2])

    ax_sub = fig11.add_subplot(221)  # General Plot
    ax_sub.bar(df_seismic['Group Type'].unique(), df_seismic['Group Type'].value_counts(),
               alpha = 0.5)
    ax_sub.set_ylabel('Frequency')
    ax_sub.set_xlabel('Material Contact Group')

    filter_dict = {222:'CrackMode', 223:'CrackType', 224:'Group Type'}
    for key, filter_criteria in filter_dict.items():
        ax_sub = fig11.add_subplot(key)  # General Plot

        # ax1.plot(df_combined['Strain y from platen (%)'], df_combined['Axial stress (MPa)'], linestyle="-")

        ax_sub.set_xlabel("Time Step (-)")
        ax_sub.set_xlim(xmin=0, xmax=x_max)
        ax_sub.set_ylabel("Cumulative Crack Count")

        un_crack = df_seismic[filter_criteria].unique()

        for idx, var in enumerate(un_crack):
            by = [df_combined[var].cumsum()[i:i + big_bins].iloc[-1] for i in range(0, int(x_max), big_bins)]

            if idx == 0:
                ax_sub.bar(bb, by, width=width, label=var, alpha=0.5)
                ll = by
            else:
                ax_sub.bar(bb, by, bottom=ll, width=width, label=var, alpha=0.5)
                ll = [x + y for x, y in zip(ll, by)]

        ax_sub.set_ylim(ymin=0)
        ax_sub.legend(title=filter_criteria)
    fig11.tight_layout()
    plt.savefig(os.path.join(folder, "histogram.svg"))
    plt.savefig(os.path.join(folder, "histogram.pdf"))
    plt.savefig(os.path.join(folder, "histogram.png"), dpi=100)
    fig11.show()



if __name__ == "__main__":
    try:
        post_processing_folder_name = "/external/Size_7/post_processing"
        load_process_file(post_processing_folder_name, 10, 100)

    except KeyboardInterrupt:
        exit("TERMINATED BY USER")
