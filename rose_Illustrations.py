# /////////////////////////////////////////////////////////////// #
# !python2.7
# -*- coding: utf-8 -*-
# Python Script initially created on 2019-04-17 
# Compiled by Aly @ Grasselli Geomechanics Group, UofT, 2019 
# Created using PyCharm 
# Current Version - Dated Apr 23, 2018
# /////////////////////////////////////////////////////////////// #

import formatting_codes
import sys
import os
import math
from windrose import WindroseAxes
# from windrose_edit import WindroseAxes2
import matplotlib.pyplot as plt
import matplotlib
import numpy
import subprocess

'''
Function for 2D Rose Diagrams
    Inputs:
        - Angle of Crack Element
        - Failure Mode
        - Length of Crack Element
    Returns:
        OR - A) 2D Rosette of Initial DFN, if exists
        OR - B) 2D Rosette of Mesh
        OR - C) 2D Rosette of Crack Elements as simulation progresses. 
'''



'''
Function to find max/min bin for cracks
    Inputs:
        - Reads all the information of the Angle of the crack elements
    Returns:
        - Max and Min bins of the cracks through the entire simulation
        - Interval for the major axis division
'''


def illus_var(dam):
    global max_range, min_range, inter
    if (max(dam) - min(dam)) > 1:
        max_range = int(max(dam) + 1)
        min_range = int(min(dam))
        inter = 1
    else:
        max_range = (max(dam))
        min_range = (min(dam) - 0.1)
        inter = 0.1
    return min_range, max_range, inter

def rose_illustration(post_processing, width, height, no_of_lines, reset_scale, fname, fname1, ax=None, dx=0, dy=0):
    sub_post_processing = os.path.join(post_processing, "Simulation_Rose_Plot")
    temp_rose_angle = os.path.join(post_processing, 'temp_rose_angle.csv')
    temp_str_damage = os.path.join(post_processing, 'temp_str_damage.csv')
    # reset_scale = float(reset_scale)
    rose_angle, damage = [], []
    with open(temp_rose_angle, "r") as f:
        for line in f:
            rose_angle.append(float(line.strip()))
    f.close()
    with open(temp_str_damage, "r") as f:
        for line in f:
            damage.append(float(line.strip()))
    f.close()
    # print(type(rose_angle), len(rose_angle))
    # print(type(damage), len(damage))
    # exit()
    illus_var(damage)
    sub_post_processing = os.path.join(post_processing, "Simulation_Rose_Plot")

    ''' Initialise Figure '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='windrose')
    viridis = plt.get_cmap('viridis')

    # Plot Sector -90 to 90
    # Set locations of X-Axis
    ax.set_thetagrids([90, 45, 0, -45, -90], labels=['0$^o$', '', '90$^o$', '', '0$^o$'])
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    ax.set_theta_zero_location("N")
    # Location of Radial Labels
    ax.set_rlabel_position(270)

    # Format the grid/plot
    plt.grid(alpha=0.5)
    plt.tight_layout()

    ''' Generate appropriate y-axis division'''
    # Create bins and counters to find mode value. Bins of 10 degrees.
    counts, bins = numpy.histogram(rose_angle, bins=17, range=(5, 175))
    divisions = (int(str(1).ljust(len(str(max(counts))) - 1, '0')))
    count = (int(math.ceil((max(counts) + divisions) / divisions)) * divisions)
    grid_labels = numpy.linspace(0, count, 11, endpoint=True, dtype=int)

    ''' A) 2D Rosette of Initial DFN, if exists '''
    if fname == "2D_initial_damage_intensity_rosette":

        resultant_force_vector = math.sqrt((dx * dx) + (dy * dy))
        resultant_force_angle = (math.degrees(math.atan2(dy, dx)))
        r_text = "%.2f $mm$ at %.2f$^{o}$" % (resultant_force_vector, resultant_force_angle)

        ax.bar(rose_angle, damage, normed=False, opening=0.8, edgecolor='white', nsector=36, cmap=viridis,
               bins=numpy.arange(min_range, max_range, inter))

        # Set Legend And Legend Properties
        legend = ax.set_legend(title="Length", bbox_to_anchor=(0.5, 0.05), loc='center', ncol=4)
        legend.get_title().set_fontsize('10')
        legend.get_title().set_weight('bold')
        plt.setp(plt.gca().get_legend().get_texts(), fontsize='10')  # legend 'list' fontsize

        # Set locations of Y-Axis & Label
        ax.set_yticks(grid_labels)
        ax.tick_params(axis='y', pad=7.5, rotation=45)
        ax.set_yticklabels(grid_labels, verticalalignment="top", horizontalalignment='center',
                           size=8)  # works ONLY for on matplotlib 2.2.3+
        ax.text(0.75, 0.15, 'Frequency', horizontalalignment='center', verticalalignment='center', size=10, weight="bold",
                transform=ax.transAxes)

        # Figure Information
        width = float(width)
        height = float(height)
        plt.suptitle("Fracture Intensity $P_{21}$ %.2f $mm^{-1}$\nResultant Vector %s" % (sum(damage) / (width * height), r_text),
                     fontsize=14, y = 0.9, size = 10, weight="bold")
        plt.savefig(os.path.join(post_processing, fname), dpi=600, bbox_inches = 'tight', pad_inches = 0)

        ''' B) 2D Rosette of Mesh '''
    elif fname == "2D_mesh_rosette":
        ax.bar(rose_angle, damage, normed=False, opening=0.8, edgecolor='white', nsector=36, cmap=viridis,
               bins=numpy.arange(min_range, max_range, inter))

        # Set Legend And Legend Properties
        legend = ax.set_legend(title="Length", bbox_to_anchor=(0.5, 0.05), loc='center', ncol=4)
        legend.get_title().set_fontsize('10')
        legend.get_title().set_weight('bold')
        plt.setp(plt.gca().get_legend().get_texts(), fontsize='10')  # legend 'list' fontsize

        # Set locations of Y-Axis & Label
        ax.set_yticks(grid_labels)
        ax.tick_params(axis='y', pad=7.5, rotation=45)
        ax.set_yticklabels(grid_labels, verticalalignment="top", horizontalalignment='center',
                           size=8)  # works ONLY for on matplotlib 2.2.3+
        ax.text(0.75, 0.15, 'Frequency', horizontalalignment='center', verticalalignment='center', size=10,
                weight="bold", transform=ax.transAxes)

        # Figure Information
        plt.suptitle("Mesh Plot - %s lines" % no_of_lines, fontsize=16, y = 0.85, size = 10, weight="bold")
        plt.savefig(os.path.join(post_processing, fname), dpi=600)

        ''' C) 2D Rosette of Crack Elements as simulation progresses. '''
    else:
        fname = float(fname)
        reset_scale = int(reset_scale)
        # Process a Paraview screenshot
        # screenshot(post_processing, fname)

        ''' Generate appropriate y-axis division '''
        divisions = (int(str(1).ljust(len(str(reset_scale)) - 1, '0')))
        count = (int(math.ceil((reset_scale + divisions) / divisions)) * divisions)
        grid_labels = numpy.linspace(0, count, 11, endpoint=True, dtype=int)

        # print (reset_scale, divisions, count)
        ax.bar(rose_angle, damage, normed=False, opening=0.85, edgecolor='white', nsector=36, cmap=viridis,
                   bins=[1., 1.5, 2.01])

        # Set Legend And Legend Properties
        legend = ax.set_legend(labels=['Tensile Dominant', 'Shear Dominant', 'Mixed Mode'], title="Failure Mode",
                               bbox_to_anchor=(0.5, 0.05), loc='center', ncol=3)
        legend.get_title().set_fontsize('10')
        legend.get_title().set_weight('bold')
        plt.setp(plt.gca().get_legend().get_texts(), fontsize='10')  # legend 'list' fontsize

        # Set locations of Y-Axis & Label
        ax.set_yticks(grid_labels)
        ax.tick_params(axis='y', pad=7.5, rotation=45)
        ax.set_yticklabels(grid_labels, verticalalignment="top", horizontalalignment='center',
                           size=8)  # works ONLY for on matplotlib 2.2.3+
        ax.text(0.75, 0.15, 'Frequency', horizontalalignment='center', verticalalignment='center', size=10, weight="bold",
                transform=ax.transAxes)

        # Figure Information
        plt.suptitle("Time Step %s - %s cracks" % (str("%0.4d" % int(fname)), len(damage)), fontsize=16, y=0.9, size=10, weight="bold")

        # Create subfolder
        if not os.path.exists(sub_post_processing):  # Check to see if the folder exists
            os.makedirs(sub_post_processing)  # if not then makes the folder
        # Save figure to folder
        rose_name = os.path.join(sub_post_processing,
                                 ("Time_Step_" + str("%0.4d" % int(fname)) + "_2D_final_damage_rosette.png"))
        plt.savefig(rose_name, dpi=600)

        # Start merge of this figure and the Paraview Screen shot
        # Pass to "python" for the PIL module
        # if fname == alist[-1]:
        #     fname1 = "FinalFrame"
        # else:
        #     fname1 = "NIL"
        # Screen shot folder
        Screen = os.path.join(post_processing, "Screenshots")
        # Create subfolder
        if not os.path.exists(Screen):  # Check to see if the folder exists
            os.makedirs(Screen)  # if not then makes the folder
        # Save Screen shot to folder
        ffname = os.path.join(Screen, (str("%0.4d" % int(fname)) + ".png"))
        # command = ("python /hdd/home/aly/Desktop/Dropbox/Python_Codes/irazu_post_processing/merge_images.py %s %s %s %s %s" % (str(post_processing) , str(rose_name) , str(ffname) , str(fname) , str(fname1)))
        import merge_images
        merge_images.merge_images(str(post_processing) , str(rose_name) , str(ffname) , str(fname) , str(fname1))
        # print command
        # subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        # import merge_images
        # merge_images.merge_images(post_processing, rose_name, ffname, fname, fname1)

    plt.close('all')

    '''
    # In case convert locks up your system
    # https://github.com/phw/peek/issues/112
    locate this file => /etc/ImageMagick-6/policy.xml
    < policy domain = "resource" name = "memory" value = "2GiB" / >
    < policy domain = "resource" name = "disk" value = "1GiB" / >
    '''
    # At final timestep
    # Create an animated gif from the stitched images
    if fname == "FinalFrame":
        # print(red_text("\nCreating Animation"))
        os.system("convert -delay 100 -dispose Background " + str(sub_post_processing) + "/*.png -loop 1 " + str(
            sub_post_processing) + "/animation.gif")


if __name__ == "__main__":
    try:
        import ast
        # a = ast.literal_eval(sys.argv[1])
        # b = ast.literal_eval(sys.argv[6])
        # c = ast.literal_eval(sys.argv[7])
        rose_illustration(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5], sys.argv[6],sys.argv[7],sys.argv[8])
        # rose_illustration(analysis_directories)
    except KeyboardInterrupt:
        exit("TERMINATED BY USER")