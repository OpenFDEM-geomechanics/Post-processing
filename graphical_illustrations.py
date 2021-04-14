# /////////////////////////////////////////////////////////////// #
# !python2.7
# -*- coding: utf-8 -*-
# Python Script initially created on 2019-04-17 
# Compiled by Aly @ Grasselli Geomechanics Group, UofT, 2019 
# Created using PyCharm 
# Current Version - Dated Apr 23, 2018
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
        - Loads the temporary dictionary files
        - Processes the data
        - Generate the graphs based on input
    Output:
        - ax1 = Histogram chart of broken Phase Interfaces
        - ax2 = Stacked Bar Chart different broken Phase Interfaces based on Time Step
        - ax3 = Stacked Bar Chart of Crack Type (IntER/IntRA granular cracks)
        - ax4 = Stacked Bar Chart for Fracture Mode Type (Tensile Dominant / Shear Dominant / Mixed Mode)
'''

import sys
import os
import csv
import json
import matplotlib
import matplotlib.pyplot as plt
import numpy
import operator
import formatting_codes
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

failure_mode = {1: "Mode 1", 2: "Mode 2", 3: "Simultaneous"}
crack_dir = {"intergranular": 1, "intragranular": 2, "intergranular / intragranular": 3}


def plot_bar_from_counter(post_processing, begin, end, bin_freq, graph_name, ):

    ''' Default MATPLOTLIB Fonts '''
    # plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams["figure.figsize"] = [5, 5]
    matplotlib.rcParams['font.family'] = ['arial']
    matplotlib.rcParams['font.size'] = 8

    # Temporary files
    temp_counter = os.path.join(post_processing, 'temp_counter.csv')  # AX1
    temp_failures = os.path.join(post_processing, 'temp_failures.csv')  # AX2
    temp_types = os.path.join(post_processing, 'temp_types.csv')  # AX3
    temp_break = os.path.join(post_processing, 'temp_break.csv')  # AX4

    #Load data from files
    # AX1 Data
    with open(temp_counter, mode='r') as infile:
        reader = csv.reader(infile)
        counter = {rows[0]: int(rows[1]) for rows in reader}
    # print(type(counter), counter)

    # AX2 Data
    with open(temp_failures) as f:
        fail_series = json.load(f)
    # print("FAIL SERIES", fail_series)

    # AX3 Data
    with open(temp_types) as f:
        fail_type_series = json.load(f)
    # print("FAIL SERIES TYPE", fail_type_series)

    # AX4 Data
    with open(temp_break) as f:
        break_series = json.load(f)
    # print("BREAK SERIES TYPE", break_series)


    # Setup divisions
    alist = range(int(begin), int(end), int(bin_freq))

    try:
        from cycler import cycler
    except AttributeError:
        print("matplotlib may not be the latest version. Rebuild/Update matplotlib and cycler, to try to resolve this.")
        exit("Can not generate graphs. Terminating")
    fig = plt.figure(figsize=[12.8, 7.2])
    if graph_name == "main_graph":
        # if ax1 is None:
        fig = plt.figure(figsize=[12.8, 7.2])
        ax1 = fig.add_subplot(221) # General Plot
        ax2 = fig.add_subplot(222) # Fracture Group No.
        ax3 = fig.add_subplot(223) # Crack Type
        ax4 = fig.add_subplot(224) # Fracture Mode Type

        ''' AX1 '''

        frequencies = counter.values()
        names = counter.keys()

        x_coordinates = numpy.arange(len(counter))
        ax1.bar(x_coordinates, frequencies, align='center', alpha= 0.5)

        ax1.xaxis.set_major_locator(plt.FixedLocator(x_coordinates))
        ax1.xaxis.set_major_formatter(plt.FixedFormatter(names))

        ax1.set_ylabel("Frequency")
        ax1.set_xlabel("Material Contact")

        ''' AX2 '''
        N = [len(x) for x in fail_series.values()][0]

        # If only 1 contact type found. Increase by 1 to avoid DivisionError.
        if len(fail_series.keys()) == 1:
            colors = 2
        else:
            colors = len(fail_series.keys())
        cmap = plt.get_cmap('viridis') # Refer color_maps Python documentation for more maps!
        ax2.set_prop_cycle(cycler('color', [cmap(1. * (i) / (colors - 1)) for i in range(colors)]))

        bottom, x = numpy.zeros(N), range(N)

        for key in fail_series:
            ax2.bar(x, fail_series[key], alpha= 0.5, label=key, bottom=bottom)
            bottom += numpy.array(fail_series[key])

        ax2.set_xticks(numpy.arange(N))
        ax2.set_xticklabels(numpy.array(alist), rotation=45)
        ax2.set_ylabel("Frequency")
        ax2.set_xlabel("Time Step")
        ax2.legend()

        ''' AX3 '''

        if len(fail_type_series.keys()) == 1:
            colors_ax3 = 2
        else:
            colors_ax3 = len(fail_type_series.keys())
        ax3.set_prop_cycle(cycler('color', [cmap(1. * (i) / (colors_ax3 - 1)) for i in range(colors_ax3)]))
        P = [len(x) for x in fail_type_series.values()][0]

        bottom, x = numpy.zeros(P), range(P)

        for key in fail_type_series:
            ax3.bar(x, fail_type_series[key], alpha=0.5, label=key, bottom=bottom)
            bottom += numpy.array(fail_type_series[key])

        ax3.set_xticks(numpy.arange(P))
        ax3.set_xticklabels(numpy.array(alist), rotation=45)
        ax3.set_ylabel("Frequency")
        ax3.set_xlabel("Time Step")
        ax3.legend()

        ''' AX4 '''

        colors_ax4 = len(failure_mode.keys()) + 1
        ax4.set_prop_cycle(cycler('color', [cmap(1. * (i) / (colors_ax4 - 1)) for i in range(colors_ax4)]))

        M = [len(x) for x in break_series.values()][0]
        bottom, x = numpy.zeros(M), range(M)

        for key in break_series:
            ax4.bar(x, break_series[key], alpha= 0.5, label=key, bottom=bottom)
            bottom += numpy.array(break_series[key])

        ax4.set_xticks(numpy.arange(M))
        ax4.set_xticklabels(numpy.array(alist), rotation=45)
        ax4.set_ylabel("Frequency")
        ax4.set_xlabel("Time Step")
        ax4.legend()

        plt.suptitle("Time Window %s to %s" % (begin, end))
        plt.savefig(os.path.join(post_processing, "histogram.svg"))
        plt.savefig(os.path.join(post_processing, "histogram.pdf"))
        plt.savefig(os.path.join(post_processing, "histogram.png"), dpi=100)

    plt.close('all')

'''
PLOT CURVES.
    Inputs:
        - File containing stress/strain information
        - Simulation type
    Output:
        - Stress/Strain & Stress/Volumetric Stain Curve in various formats.
        - Stress/Strain Curve in various formats.
'''

def single_graph(post_processing, sim_type):

    # Variables
    axial_strain, stress, event_count, volumetric_strain = [], [], [], []
    # Name of file to look for
    names = ["history.csv"]
    filename = os.path.join(post_processing, names[0])
    # Configure name of file (Name of sub_folder)
    pdf_name = os.path.basename(os.path.dirname(post_processing))
    # Read file to obtain necessary data
    with open(filename) as history_file:
        reader = csv.DictReader(history_file, delimiter=',')
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        for row in reader:
            if sim_type == "UCS Simulation":
                axial_strain.append(float(row['strain y from platen (%)']))
                stress.append(float(row['axial stress (mpa)']))
                event_count.append(float(row['event count']))
                volumetric_strain.append(float(row['volume (mm3)']))
            else:
                axial_strain.append(float(row['y displacement(mm)']))
                stress.append(float(row['indirect tensile stress (mpa)']))
    history_file.close()

    ind_count = numpy.diff(event_count)
    ind_count = [0] + ind_count.tolist()
    print(ind_count)
    print(type(ind_count))

    '''
    DRAW STRESS/STRAIN CURVE
    '''

    fig2 = plt.figure(figsize=[12.8, 7.2])
    ax1 = fig2.add_subplot(1, 1, 1)

    ax1.plot(axial_strain, stress, linestyle="-", label=pdf_name)
    ax1.set_title(pdf_name)
    ax1.set_xlabel("Strain (%)")
    ax1.set_xlim(xmin=0)
    ax1.set_ylabel("Stress (MPa)")
    ax1.set_ylim(ymin=0)

    # Save in the relevant sub-directory directory
    plt.tight_layout()
    fig2.savefig(os.path.join(post_processing, "stress_strain.pdf"))
    fig2.savefig(os.path.join(post_processing, "stress_strain.svg"))
    fig2.savefig(os.path.join(post_processing, "stress_strain.png"), dpi=100)

    '''
    DRAW STRESS/STRAIN CURVE
    '''

    fig2 = plt.figure(figsize=[12.8, 7.2])
    ax1 = fig2.add_subplot(1, 1, 1)

    ax1.plot(axial_strain, stress, linestyle="-", label=pdf_name)
    ax1.set_title(pdf_name)
    ax1.set_xlabel("Strain (%)")
    ax1.set_xlim(xmin=0)
    ax1.set_ylabel("Stress (MPa)")
    ax1.set_ylim(ymin=0)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel('Cumulative AE Event Count', color='red')
    ax2.plot(axial_strain, event_count, color='red')
    ax2.set_ylim(ymin=0)

    binnings = 10
    b= [sum(axial_strain[i:i+binnings])/binnings for i in range(0, len(axial_strain), binnings)]
    bx = [sum(ind_count[i:i+binnings])/binnings for i in range(0, len(ind_count), binnings)]

    ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    # right, left, top, bottom

    ax3.set_ylim(ymin=0)
    ax3.set_ylabel('Individual AE Event Count', color='green')
    # ax3.plot(axial_strain, ind_count, color='green')
    ax3.bar(b, bx, color='green', width=0.005, alpha=0.5)
    ax3.set_ylim(ymin=0, ymax=max(bx)*1.1)
    ax3.spines['right'].set_position(('outward', 60))

    # Save in the relevant sub-directory directory
    fig2.savefig(os.path.join(post_processing, "stress_strain_with_count.pdf"))
    fig2.savefig(os.path.join(post_processing, "stress_strain_with_count.svg"))
    fig2.savefig(os.path.join(post_processing, "stress_strain_with_count.png"), dpi=100)


    '''
    DRAW STRESS/VOLUMETRIC STRAIN CURVE
        - Only for UCS Simulations
    '''

    if sim_type == "UCS Simulation":
        # slice the volumetric strain
        initial_volume = float(volumetric_strain[0])  # Initialise Initial Volume
        volumetric_strain_calculated = [((i - initial_volume) * 100 / initial_volume) for i in
                                        volumetric_strain]  # Calculate volumetric strain (dV / V0)
        volumetric_strain_limit = float(max(axial_strain)) * 1.1  # Limit of slice at 110% of max Axial Strain
        volumetric_strain_sliced = [(i * -1) for i in volumetric_strain_calculated if
                                    i <= volumetric_strain_limit]  # Slice volumetric strain
        stress_sliced = stress[:len(volumetric_strain_sliced):]  # Slice stress

        # Check is Simulation results OR EXP results
        fig1 = plt.figure(figsize=[12.8, 7.2])
        ax = fig1.add_subplot(1, 1, 1)

        ax.plot(axial_strain, stress, linestyle="-", label=pdf_name)
        ax.plot(volumetric_strain_sliced, stress_sliced, linestyle=":", c='red')
        # Manipulate Figure
        ax.set_title(pdf_name)
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_xlabel("Strain (%)")
        ax.set_ylabel("Stress (MPa)")

        # Save in the relevant sub-directory directory
        fig1.savefig(os.path.join(post_processing, "stress_vol_strain.pdf"))
        fig1.savefig(os.path.join(post_processing, "stress_vol_strain.svg"))
        fig1.savefig(os.path.join(post_processing, "stress_vol_strain.png"), dpi=100)

    try:
    # Obtain the Max_Stress in each file and its relative index in the list
        max_value_index, max_values_stress =  max(enumerate(stress), key=operator.itemgetter(1)) # Enumerate to get index of max value in list
        elastic_modulus = ((stress[int(max_value_index) / 2]) / (axial_strain[int(max_value_index) / 2]) / 10 )
        print(formatting_codes.green_text("Simulation Results"))
        print("\tMax values at Output No. %s\n\tMax Stress %10.2f\n\tElastic Modulus %10.2f" % (max_value_index, max_values_stress, elastic_modulus))
    except ZeroDivisionError:
        print(formatting_codes.red_text("There is a problem in the history.csv file\nThe stress/strain values are extremely low."))

# Function chooser
func_arg = {"-histogram": plot_bar_from_counter, "-graphs": single_graph}

if __name__ == "__main__":
    try:
        if sys.argv[1] == "-histogram":
            func_arg[sys.argv[1]](sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],)
        elif sys.argv[1] == "-graphs":
            func_arg[sys.argv[1]](sys.argv[2], sys.argv[3])
    except KeyboardInterrupt:
        exit("TERMINATED BY USER")