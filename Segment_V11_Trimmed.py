# /////////////////////////////////////////////////////////////// #
# !python2.7
# -*- coding: utf-8 -*-
# Python Script initially created on 15/03/17
# Compiled by Aly @ Grasselli's Geomechanics Group, UofT, 2017
# Base Script by Omid Mahdabi
# Created using PyCharm 
# Current Version 11 - Dated Mar 21, 2018
# /////////////////////////////////////////////////////////////// #

# THINGS TO IMPROVE
    # 2) Integrate Qi Matlab Code for "seismic tool" <+> Maybe work on returning the file is CSV for Matlab (V#07)? # Needs to be verified
    # 4) Monitor relative movement of the broken joints => Each Element is not tracked
    # 5) If the initial model has DFN's. Subtract that from the total.

# UPDATES
# V01
    # -  Seismic Clustering Added.
    # -  Graphical representation of crack type.
# V02
    # -  Graphical representation of crack type - Enhanced
    # -  Input routine enhanced
    # -  Graphical display detaches as a Child Window
# V03
    # -  Added Graphical Display (Failure Mode Type)
    # -  Code Enhancements & General Cleanup
    # -  Percentage Materials displayed (takes into account the platens)
#V04
    # -  Locate Broken Joint along DFN
    # -  Added Failure_Mode dictionary & Changed algorithm
    # -  Added crack type dictionary
    # -  Time stamp of each output streamlined
    # -  Color scale for ax3 and ax4 amended
    # -  Range of Color Scheme rectified
    # -  dfn_line_list has set to allow only unique values
    # -  pvpython cycler error highlighted.
#VO5
    # -  seismicclustering function updated
        # 1) To read the #11 (12) ID in the reader
        # 2) Return the time step as an integer so it can do the max/min operation.
    # -  Allow error free execution if seismic data is not available
#V06
    # - Visualization code runs separately from the rest of the code to overcome pvpython/python limitation on Glass computers.
#V07
    # - Optimized Timer Display
    # - Optimized the pvpython/python Visualization problem
    # - Better Validations and Terminates if no broken joints
    # - Checks for missing modules prior to execution to save on processing time
#V08
    # - Table name inserted instead of Array Number
    # - Identifies 2D or 3D Simulation
    # - Based on 2D or 3D, skip variable added to iterate over points of the broken joints
    # - Batch processing added
    # - Return 'Area' instead of length

#V09
    # - Returns the angle of the crack (orientation)
    # - WindRose (Radar) Diagram of the crack orientation binned based on failuremode. Bins can be changed to parameter of choice.

#V11
    # - os.walk() added to collect all entire directory tree.
    # - gather information on the file being processed
    # - create a "post_processing" folder that saves all the files to it.

try:
    paraview.simple
except:
    from paraview.simple import *
import itertools, operator, csv, argparse, subprocess
import numpy, os, re, sys, time, math, gc

import scipy.spatial as ss
import paraview.vtk as vtk
from collections import Counter
from shapely.geometry import Polygon
from windrose import WindroseAxes
import matplotlib.pyplot as plt
from operator import itemgetter
try:
    from cycler import cycler
except AttributeError:
    print("pvpython and matplotlib may not be linked properly. Try running \n$ sudo apt-get install paraview-python \nElse rebuild matplotlib and cycler, to try to resolve this.")

# Attempt to overcome affinity problem related to numpy (OPENBLAS_MAIN_FREE)
# by forcing the processes to be binned across CPUS's
# I can verify that all the cores are being utilized on my machine while the code is running by using the htop console program. In some cases, Python modules like Numpy, scipy, etc. may limit processes in Python to running on only one core on Linux machines, which defeats the purpose of writing concurrent code in this case. (See for example this discussion.) To fix this, we can import Python's os module to reset the task affinity in our code:
import multiprocessing
os.system("taskset -cp %d" % (os.getpid()))

# TRY Block for illustrations
# named_libs = [('matplotlib', 'mpl')]
# for (name, short) in named_libs:
#     try:
#         lib = __import__(name)
#     except:
#         var_matplotlib = name
#     else:
#         var_matplotlib = ''
#         globals()[short] = lib
#         import matplotlib as mpl
#         # mpl.use("TkAgg") # Force backend to TkAgg to avoid verbose error
#         from windrose import WindroseAxes
#         # from matplotlib import pyplot as plt
#         import matplotlib.cm as cm

# paraview.simple._DisableFirstRenderCameraReset()


# /// DICTIONARIES /// #

crack_dir = {"intergranular": 1, "intragranular": 2, "intergranular / intragranular": 3}
failure_mode = {1: "Mode 1", 2: "Mode 2", 3: "Simultaneous"}

numbers = re.compile(r'(\d+)')

# /// TIMER FUNCTION /// #

def calc_timer_values(end_time):
    minutes, sec = divmod(end_time, 60)
    if end_time < 60:
        return ("\033[1m%.2f seconds\033[0m" % end_time)
    else:
        return ("\033[1m%d minutes and %d seconds\033[0m." % (minutes, sec))


# /// ADMINISTRATIVE AND SORTING OF FILES IN FOLDER /// #


def red_text(val):  # RED Bold text
    tex = "\033[1;31m%s\033[0m" % val
    return tex


def green_text(val):  # GREEN Bold text
    tex = "\033[1;92m%s\033[0m" % val
    return tex


def bold_text(val):  # Bold text
    tex = "\033[1m%s\033[0m" % val
    return tex

# Strip the numerical portion of the file
def numericalSort(value):
    parts = numbers.split(value)  # Split the numerical part of the file
    parts[1::2] = map(int, parts[1::2])  # Return the numerical portion of the file
    return parts

# Attain the list of files in the directory then sort them based on their numerical value
def findOutputFiles(dir_path, file_extension, output_type):
    os.chdir(dir_path) # cd to the base directory path
    list_of_files = []
    # Get the list of files/directories
    files = os.listdir('.')
    files.sort()
    for f in files:
    # Filter only the files
        if os.path.isfile(f):
        # Clean the list to only include files with the given extension that contain output_type
            if f.endswith(file_extension):
                if output_type in f:
                    list_of_files.append(os.path.join(dir_path, f))
    # Sort list of files by their numerical value
    list_of_files = sorted(list_of_files, key=numericalSort)
    return list_of_files

''' 
Print iterations progress
Source: https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
'''
def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=50):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(bar_length * (100. / bar_length) * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '/' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()

# /// Find all the subdirectories of a given path /// #
# Inputs:
# - Starting directory (full path is required)
# Returns:
# - A list of subdirectories

def findSubdirectories(dir_path):
    sub_dirs = []
    for root, dirs, files in os.walk(dir_path):
        for dir_name in dirs:
            # print (os.path.join(root,name))
            sub_dirs.append(os.path.join(root,dir_name))
    sub_dirs = sorted(list(set(sub_dirs))) # Sort directories alphabetically in ascending order
    print("Found \033[1m%s\033[0m sub-directories" % len(sub_dirs))
    return sub_dirs

# /// Create post_processing folder /// #

def create_post_processing(dir_path):
    global post_processing
    post_processing = os.path.join((dir_path), "post_processing")
    if not os.path.exists(post_processing): # Check to see if the folder exists
        os.makedirs(post_processing) # if not then makes the folder
    return post_processing

# /// PROCESS PARAVIEW OUTPUTS /// #
    # Inputs:
        # Directory where output files are located

def processOutputPropID(dir_path, file_extension, specimen_center, platen_cells_prop_id, platen_points_prop_id, output_file_name):
    global list_of_files, list_of_files_broken, list_of_files_dfn, list_of_files_seismic
    start_prop = time.time() # Processing Start Time
    #  List BASIC Files => Fatal Error
    list_of_files = findOutputFiles(dir_path, file_extension, "_basic_")
    #  If no output files were found, warn the user and abort
    if list_of_files is None or not list_of_files:
        print(bold_text("\nERROR: ") +
            "No output files matching extension '%s' containing '%s' were found in %s\n%s\n" % (file_extension, "_basic_", dir_path, red_text("ABORTING")))
        return
    # List BROKEN JOINTS Files => Fatal Error
    list_of_files_broken = findOutputFiles(dir_path, file_extension, "_broken_joint_")
    #  If no output files were found, warn the user and abort
    if list_of_files_broken is None or not list_of_files_broken:
        print(bold_text("\nERROR: ") +
            "No output files matching extension '%s' containing '%s' were found in %s\n%s\n" % (file_extension, "_broken_joint_", dir_path, red_text("ABORTING")))
        return
    # List SEISMIC Files => Warning Error
    global list_of_files_seismic
    list_of_files_seismic = findOutputFiles(dir_path, file_extension, "_seismic_event_")
    #  If no output files were found, warn the user and continue
    if list_of_files_seismic is None or not list_of_files_seismic:
        print(bold_text("\nERROR: ") +
            "No output files matching extension '%s' containing '%s' were found in %s\nProceeding\n"
            % (file_extension, "_seismic_event_", dir_path))
        list_of_files_seismic = None
    # List DFN Files => Warning Error
    list_of_files_dfn = findOutputFiles(dir_path, file_extension, "_dfn_property_")
    #  If no output files were found, warn the user and continue
    if list_of_files_dfn is None or not list_of_files_dfn:
        print(bold_text("\nERROR: ") +
            "No output files matching extension '%s' containing '%s' were found in %s\n\033[1;31;0mWill not identify DFN and/or intragranular cracks\033[0m\nProceeding\n"
            % (file_extension, "_dfn_property_", dir_path))

    print("Processing data files in:" + red_text(dir_path))
    print("Found %d output files (=time steps)" % len(list_of_files))

    # return list_of_files, list_of_files_broken, list_of_files_dfn, list_of_files_seismic


    global output, output_broken, output_dfn, output_seismic, cells, points, broken_joint_info, data_object, pv_grid_reader, pv_grid_reader_broken, pv_grid_reader_seismic, pv_grid_reader_dfn, seismic_point_info

    # Read the vtu files using an XMLUnstructuredGridReader
    # print list_of_files
    pv_grid_reader = XMLUnstructuredGridReader(FileName=list_of_files)
    pv_grid_reader_broken = XMLUnstructuredGridReader(FileName=list_of_files_broken)
    if list_of_files_seismic:
        pv_grid_reader_seismic = XMLUnstructuredGridReader(FileName=list_of_files_seismic)
    if list_of_files_dfn:
        pv_grid_reader_dfn = XMLUnstructuredGridReader(FileName=list_of_files_dfn)

    # Show the grid
    Show(pv_grid_reader)
    Show(pv_grid_reader_broken)
    if list_of_files_seismic:
        Show(pv_grid_reader_seismic)
    if list_of_files_dfn:
        Show(pv_grid_reader_dfn)

    # Define a new cell array for the output
    cellArray = vtk.vtkCellArray()

    # Get the Client Side Data Object of the grid reader
    data_object = pv_grid_reader.GetClientSideObject()
    data_object_broken = pv_grid_reader_broken.GetClientSideObject()
    if list_of_files_seismic:
        data_object_seismic = pv_grid_reader_seismic.GetClientSideObject()
    if list_of_files_dfn:
        data_object_dfn = pv_grid_reader_dfn.GetClientSideObject()

    # Get the output
    output = data_object.GetOutput()
    output_broken = data_object_broken.GetOutput()
    if list_of_files_seismic:
     output_seismic = data_object_seismic.GetOutput()
    if list_of_files_dfn:
        output_dfn = data_object_dfn.GetOutput()

    # Get the required information of points and cells
    cells = output.GetCellData()
    points = output.GetPointData()
    broken_joint_info = output_broken.GetCellData()
    if list_of_files_seismic:
        seismic_point_info = output_seismic.GetPointData()
    if list_of_files_dfn:
        dfn_point_info = output_dfn.GetPointData()
        dfn_joint_info = output_dfn.GetCellData()

    do_strain_gauge_analysis = []
    top_platen_cellArray = vtk.vtkCellArray()
    bottom_platen_cellArray = vtk.vtkCellArray()

    # Get the available time steps
    # Set global variables
    global subId, pcoords, w, node_skip, time_steps, begin, end, mineral_bound_2, angle_deg, list_azimuth, list_dip

    time_steps = pv_grid_reader.TimestepValues
    cell_pr_txt = 'property_id'
    point_pr_txt = 'property_id'
    if cells.GetArray('property_id') is None:
        cell_pr_txt = 'material property ID'
    if points.GetArray('property_id') is None:
        point_pr_txt = 'boundary condition ID'
    # print 'time steps', time_steps

    # Check the current view time
    view = GetActiveView()
    view.ResetCamera()  # Reset the camera to display the entire model.

    # Check 2D (3 Points - Triangle) or 3D (4 Points - Tetrahedral) Simulations from the cell vertex.
    check_dim_cell = output.GetCell(0)
    number_of_points_per_cell = float(check_dim_cell.GetNumberOfPoints())
    if number_of_points_per_cell == 3: #2D
        print (green_text("2D Simulation"))
        global node_skip, skip
        node_skip = 4
        skip = 4
    else: #3D
        print (green_text("3D Simulation"))
        node_skip = 6
        exit("3D Simulation not supported")

    # Variables

    # global subId, pcoords, w, node_skip
    subId = vtk.mutable(0)  # Dummy for FindCell
    pcoords = [0.0, 0.0, 0.0]  # Dummy for FindCell
    w = [0.0, 0.0, 0.0]  # Dummy for FindCell
    idlist = vtk.vtkIdList() # Dummy for idList
    error_array, cluster, failure_modes = [], [], []  # For QAQC Error & Cluster CrackType
    mineral_bound, combos, breakdown, fail_count, crack_types, crack_count = [], [], [], [], [], []
    list_azimuth, list_dip, list_points = [], [], []

    # Prompt for user input first
    # global alist
    # begin, end = 150, 155
    global alist
    answer, answer1, answer2 = "N", "Y", "Y"
    bin_freq = 50
    begin = int(0)
    end = int(max(time_steps))
    # user_data()
    begin, end = int(begin), int(end) # initialise user inputs
    alist = range(int(begin), int(end), int(bin_freq))
    alist.pop(0)  # Remove first element of the list
    alist.append(int(end) - 1)  # Add the value "end" to the list


    ''' Initialization Data '''
    # Reads ONLY 1 time step and returns the dictionary of possible combinations
    # for t in time_steps:
    t = time_steps[0]
    view.ViewTime = t
    Render()
    sample_bound_x, sample_bound_y, sample_bound_z = [], [], []
    start_initial = time.time()
    # if t == 0:
    print(green_text("Starting Initialization"))

    dfn_line_list = []  # reset every time step to match dfn lines against respective broken joints
    cluster, failure_modes, crack_types = [], [], []  # reset for noncumulative graphs

    # Combinations of the the various mineral boundaries in the simulation results
    # Also return % composition within the model
    for i in range(0, output.GetNumberOfCells()):  # Create a list of the mineral types in the model
        mineral_bound.append(int(cells.GetArray('material property ID').GetTuple1(i)))
    mineral_bound_2 = range((min(mineral_bound)),
                            (max(mineral_bound)) + 1)  # return the max and min in that list
    const = Counter(mineral_bound)  # counts the frequency of each mineral type
    list_combo = [x for x in itertools.combinations_with_replacement(mineral_bound_2, 2)]  # Creates a list of all possible combinations
    for key, val in enumerate(list_combo):  # Enumerates the listcombo making it "Group #")
        combos.append("Group %d" % key)
        combos.append(val)
    global material_combinations
    material_combinations = dict(
        itertools.izip_longest(*[iter(combos)] * 2, fillvalue=""))  # combines, creating a dictionary
    for key, val in const.iteritems():
        per = "{0:.2f}".format(float(val) / float(output.GetNumberOfCells()) * 100)
        print ("Element ID %s - Percentage Composition %s %%" % (bold_text(key), bold_text(per)))
    print ("Please be aware of the presence of the Platens in the %'s")

    # If specimen center is not defined as input
    # Calculates the Entire model bounds
    if specimen_center is None:
        # print specimen_center
        model_bounds = output.GetBounds()  # xmin, xmax, ymin, ymax, zmin, zmax
        x_center = (model_bounds[0] + model_bounds[1]) / 2
        y_center = (model_bounds[2] + model_bounds[3]) / 2
        z_center = (model_bounds[4] + model_bounds[5]) / 2
        specimen_center = [x_center, y_center, z_center]
    else:
        model_bounds = output.GetBounds()  # xmin, xmax, ymin, ymax, zmin, zmax

    # If Strain gauges on, identify location using specimen center
    pv, ph = [],[]
    if do_strain_gauge_analysis:
        pv.append([specimen_center[0] + gauge_width / 2.0,
                   specimen_center[1] - gauge_length / 2.0,
                   0.0])
        pv.append([specimen_center[0] - gauge_width / 2.0,
                   specimen_center[1] - gauge_length / 2.0,
                   0.0])
        pv.append([specimen_center[0] + gauge_width / 2.0,
                   specimen_center[1] + gauge_length / 2.0,
                   0.0])
        pv.append([specimen_center[0] - gauge_width / 2.0,
                   specimen_center[1] + gauge_length / 2.0,
                   0.0])
        ph.append([specimen_center[0] + gauge_length / 2.0,
                   specimen_center[1] + gauge_width / 2.0,
                   0.0])
        ph.append([specimen_center[0] + gauge_length / 2.0,
                   specimen_center[1] - gauge_width / 2.0,
                   0.0])
        ph.append([specimen_center[0] - gauge_length / 2.0,
                   specimen_center[1] + gauge_width / 2.0,
                   0.0])
        ph.append([specimen_center[0] - gauge_length / 2.0,
                   specimen_center[1] - gauge_width / 2.0,
                   0.0])
    else:
        gauge_length = gauge_width = 0

    # If Platen Property ID not defined.
    if platen_points_prop_id is None:
        platen_points_prop_id = []
        for i in xrange(output.GetNumberOfCells()):
            # Get the Cell Object & Point Data
            # cell = output.GetCell(i)
            point_data = output.GetPointData()
            # Get the property_id of the current Cell
            cell_pr = int(cells.GetArray(cell_pr_txt).GetTuple1(i))
            # Get the id list of cell having the prescribed element property_id
            idlist = vtk.vtkIdList()
            id_list = output.GetCellPoints(i, idlist)
            cellidlist = ([idlist.GetId(k) for k in range(idlist.GetNumberOfIds())])  # Nodes of the CellID
            # Check if the property_id is equal to the prescribed element property_id
            if (cell_pr == platen_cells_prop_id):
                # Lookup each point of the cell and return the boundary condition ID.
                for j in range(0, len(cellidlist)):
                    # print point_data.GetArray('boundary condition ID')
                    b_cond = point_data.GetArray('boundary condition ID').GetTuple1(cellidlist[j])
                    platen_points_prop_id.append(int(b_cond))
            else:
                # Make a list of all X, Y, and Z of the points not in the element property ID
                for j in range(0, len(cellidlist)):
                    all_xs, all_ys, all_zs = output.GetPoint(cellidlist[j])
                sample_bound_x.append(all_xs)
                sample_bound_y.append(all_ys)
                sample_bound_z.append(all_zs)
        global width, height, thickness
        width = max(sample_bound_x) - min(sample_bound_x)
        height =  max(sample_bound_y) - min(sample_bound_y)
        thickness = max(sample_bound_z) - min(sample_bound_z)

    # Return unique elements (i.e., unique boundary id of platens)
    platen_points_prop_id = list(set(platen_points_prop_id))
    global sample_x_min, sample_x_max, sample_y_min, sample_y_max
    sample_x_min, sample_x_max  = min(sample_bound_x), max(sample_bound_x)
    sample_y_min,sample_y_max  = min(sample_bound_y), max(sample_bound_y)

    modelextend, cellextent = [], []
    top_platen_cells, top_platen_cell_ids  = [], []
    bottom_platen_cells, bottom_platen_cell_ids = [], []
    for i in xrange(output.GetNumberOfCells()):
        # Get the Cell object
        cell = output.GetCell(i)
        # Get the property_id of the current Cell
        cell_pr = int(cells.GetArray(cell_pr_txt).GetTuple1(i))
        # Check if the property_id is equal to the prescribed element property_id
        # Gather Platen information
        if (cell_pr == platen_cells_prop_id):
            # Get the first point of the cell
            # print cell.GetPointId(0)
            point_0_id = cell.GetPointId(0)
            point_0_pr = int(points.GetArray(point_pr_txt).GetTuple1(point_0_id))
            # Top Platen
            if (point_0_pr == platen_points_prop_id[0]):
                top_platen_cells.append(cell)
                top_platen_cell_ids.append(i)
                top_platen_cellArray.InsertNextCell(cell)
            # Bottom Platen
            elif (point_0_pr == platen_points_prop_id[1]):
                bottom_platen_cells.append(cell)
                bottom_platen_cell_ids.append(i)
                bottom_platen_cellArray.InsertNextCell(cell)
        # Gather information
        else:
            idlist = vtk.vtkIdList()
            output.GetCellPoints(i, idlist)
            cellidlist = ([idlist.GetId(k) for k in range(idlist.GetNumberOfIds())])
            for j in range(0, len(cellidlist)):
                # In 2D Read the co-ordinates of all the cells
                if node_skip == 4:
                    xs, ys, zs = output.GetPoint(cellidlist[j])
                    # Identifies if co-ordinate is on the edge.
                    if xs in [sample_x_min, sample_x_max] or ys in [sample_y_min, sample_y_max]:
                        modelextend.append(cellidlist[j])  # "Point ID on boundary"
                        cellextent.append(i)  # "cells on boundary"
                else:
                    # In 3D identifies ALL points in sample based on boundary condition.
                    if point_data.GetArray("boundary condition ID").GetTuple1(
                            i) == -1:  # Based on boundary condition (UCS). Needs to be reconsidered
                        xs, ys, zs = output.GetPoint(cellidlist[j])
                        modelextend.append(cellidlist[j])  # "Point ID on boundary"
                        cellextent.append(i)  # "cells on boundary"

    area_calculation(output, modelextend, specimen_center)

    '''
    DRAW THE ROSETTE FOR THE MESH
    - Platens excluded.
    '''
    dfn_output_types, initial_dfn_line_list = [], []
    areas = []
    create_post_processing(dir_path) # Create post-processing folder
    overall_points, rose_angle, damage = [], [], []
    print red_text("\nCreating Rosette (Excluding Platens)")
    for i in range(0, output.GetNumberOfCells()):
        # print(cells.GetArray('material property ID').GetTuple1(i), platen_cells_prop_id)
        if cells.GetArray('material property ID').GetTuple1(i) != float(platen_cells_prop_id):
            dfn_output_types.append(cells.GetArray('material property ID').GetTuple1(i))
            # Triangle points rotate anti-clockwise
            # Get Points
            point_a = output.GetPoint((i * 3) + 0)  # Get Point A of the Triangle
            point_b = output.GetPoint((i * 3) + 1)  # Get Point B of the Triangle
            point_c = output.GetPoint((i * 3) + 2)  # Get Point C of the Triangle
            # Use points to construct Triangles
            overall_points.append([point_a, point_b])  # Construct line #1 of Triangle
            overall_points.append([point_b, point_c])  # Construct line #2 of Triangle
            overall_points.append([point_c, point_a])  # Construct line #3 of Triangle

            # Use the points to calculate the area of the triangle
            tri_area = area(point_a, point_b, point_c)
            areas.append(tri_area)
            # print point_a, point_b, point_c, tri_area

    # Cleanup for lines that coincide (triangle edges)
    fset = set(frozenset(x) for x in overall_points)  # remove duplicate items from nested list.

    global new_overall_points

    new_overall_points = [list(x) for x in fset]  # return list of unique items

    # Calculate edge length and angle for each line
    for i in new_overall_points:
        ax, ay, az = list(i[0])
        bx, by, bz = list(i[1])
        angle_deg = (math.degrees(math.atan2(by - ay, bx - ax)))  # Calculate the slope of the line in degrees
        cal_len = math.sqrt((by - ay) ** 2 + (bx - ax) ** 2)  # Calculate the length of the line
        rose_angle.append(angle_deg)
        damage.append(cal_len)
        # print ax, ay, az, bx, by, bz, angle_deg, cal_len  # Uncomment to see line X Y Z Angle Length
    # dfn_output_types = list(set(dfn_output_types))
    # print("Mesh: Max Length %.2f \tMin Length %.2f\tAverage %.2f" % (max(damage), min(damage), reduce(lambda x, y: x + y, damage) / len(damage)))

    # illus_var(damage)
    rose_illustration(rose_angle, damage, "2D_mesh_rosette")
    print("Mesh Element Length:\tMax %.2f\tMin %.2f\tAverage %.2f" % (max(damage), min(damage), sum(damage) / len(damage)))
    print("Mesh Element Area:\tMax %.2f\tMin %.2f\tAverage %.2f" % (max(areas), min(areas), sum(areas) / len(areas)))
    global mesh_avg_length, mesh_avg_area
    mesh_avg_area = sum(areas) / len(areas)
    mesh_avg_length = sum(damage) / len(damage)



    # exit(500)

    '''
    DRAW THE ROSETTE FOR DFN's in the Model
    - Platens excluded.
    - Model Extents excluded.
    '''
    dfn_output_types, initial_dfn_line_list = [], []

    # illustration(output, cells, output, cells, output, cells, platen_cells_prop_id,
    #              dfn_output_types, "Mesh_rosette.pdf")

    if list_of_files_dfn is None or not list_of_files_dfn:
        print(red_text("No initial DFN found"))
    else:
        for i in range(0, output_dfn.GetNumberOfCells(), 2):
            dfn_output_types.append(dfn_joint_info.GetArray('property ID').GetTuple1(i))
        dfn_output_types = list(set(dfn_output_types))
        if len(dfn_output_types) == 1:
            for i in range(0, output_dfn.GetNumberOfPoints()):
                initial_dfn_line_list.append(output_dfn.GetPoint(i))  # Build initial dfn list based on number of points.
            initial_dfn_line_list = list(set(initial_dfn_line_list))
            # create_post_processing(dir_path)
            # out_name = "2D_initial_fracture_intensity_rosette.pdf"
            print("DFN Plot")
            illustration(output, cells, output_broken, broken_joint_info, output_dfn, dfn_joint_info, platen_cells_prop_id, dfn_output_types, "2D_initial_fracture_intensity_rosette.pdf")
        else:
            for i in range(0, output_broken.GetNumberOfCells(), 2):
                if broken_joint_info.GetArray('failure mode').GetTuple1(i) == 0:
                    ax, ay, az = point_a = output_broken.GetPoint(i * node_skip)  # Get Point A of the line
                    bx, by, bz = point_b = output_broken.GetPoint(
                        (i * node_skip) + node_skip)  # Get Point B of the line
                    if ax not in [sample_x_min, sample_x_max] and bx not in [sample_x_min,
                                                                             sample_x_max] and ay not in [
                        sample_y_min, sample_y_max] and by not in [sample_y_min, sample_y_max]:
                        cell_1 = output.FindCell(point_a, None, 0, 1e-4, subId, pcoords,
                                                 w)  # Lookup the cell of the point
                        cell_2 = output.FindCell(point_b, None, 0, 1e-4, subId, pcoords,
                                                 w)  # Lookup the cell of the point
                        material_id_cell_1 = cells.GetArray('material property ID').GetTuple1(
                            cell_1)  # Lookup the material ID of the cell
                        material_id_cell_2 = cells.GetArray('material property ID').GetTuple1(
                            cell_2)  # Lookup the material ID of the cell
                        if material_id_cell_1 != platen_cells_prop_id or material_id_cell_2 != platen_cells_prop_id:
                            initial_dfn_line_list.append(
                                output_dfn.GetPoint(i))  # Build initial dfn list based on number of points.
                            initial_dfn_line_list = list(set(initial_dfn_line_list))

            if list_of_files_dfn is None or not list_of_files_dfn or not initial_dfn_line_list:
                print("\033[1;31;mNo initial DFN found\033[0m")
            else:
                # create_post_processing(dir_path)
                illustration(output, cells, output_broken, broken_joint_info, output_dfn, dfn_joint_info, platen_cells_prop_id, dfn_output_types, "2D_initial_fracture_intensity_rosette.pdf")

    # Display some important information on Terminal
    # temp = sys.stdout
    # sys.stdout = open(os.path.join(post_processing, 'model_statistics.csv'), 'wb')
    print ("\nNo. of Points in Model: %s \nNo. of Mesh elements: %s" % (bold_text(output.GetNumberOfPoints()), bold_text(output.GetNumberOfCells())))
    print ("Model Bounds: xmin = %.2f, xmax = %.2f, ymin = %.2f, ymax = %.2f, zmin = %.2f, zmax = %.2f" % tuple(
        model_bounds))
    print ("Sample Dimensions: X = %.2f, Y = %.2f, Z = %.2f" % (width, height, thickness))
    print ("Specimen Center: X = %.2f, Y = %.2f, Z = %.2f" % tuple(specimen_center))
    if not platen_points_prop_id:
        print red_text("Unknown platen boundary conditions. Should be one of these values %s" % const.keys())
    elif len(platen_points_prop_id) == 1 or -1 in platen_points_prop_id:
        print red_text("Platen Boundary Conditions: %s" % bold_text(platen_points_prop_id)[0:])
        print red_text("Please recheck platen element property.")
    else:
        print ("Platen Boundary Conditions: %s" % bold_text(platen_points_prop_id)[0:])
    print ("No. of points on sample boundary - %s " % len(modelextend))
    print ("Initial Area/Volume: %.2f" % (init_area))
    # sys.stdout.close()
    # sys.stdout = temp
    print ("\nInitialization Complete: %s" % bold_text(calc_timer_values(time.time() - start_initial)))


    # exit(6)

    processOutputUCS(dir_path, list_of_files, width, height, specimen_center, platen_cells_prop_id, platen_points_prop_id, gauge_length, gauge_width, do_strain_gauge_analysis, modelextend, output_file_name)
    # plt.show()


'''
Find the points close to a centre line & sort them clockwise
    Inputs:
        - Data object containing the output (main point/cell data)
        - Coordinates of the center of the specimen [x, y, z]
        - Axis of loading: x=0; y=1; z=2
        - Diameter (width) of the sample in mm
        - Buffer size to find the closest point (should be smaller than element size).
    Returns:
        - A list of points IDs and coordinates [i, x, y, z] sorted clockwise
'''


# def find_sort_points(data_object_output, centre, axis_of_loading, diameter, buffer_size):
#     coordinates = []  # [point_id, x, y, z]
#     # Store unique coordinates to avoid adding duplicate nodes
#     uniqure_coords = []  # [(x, y, z)]
#     # Find the points close enough to the defined centre line (centre)
#     for i in range(data_object_output.GetNumberOfPoints()):
#         tmp_coord = data_object_output.GetPoint(i)
#         coord = tmp_coord[:3]
#         # Add the points if they're close enough to the centre
#         if (math.fabs(coord[axis_of_loading] - centre[axis_of_loading]) < buffer_size):
#             # Get the distance of this point to the center
#             # and only consider points that are farthest from the centre: exterior points and not interior
#             x, y, z = coord
#             dist = math.sqrt(math.pow(x - centre[0], 2) + \
#                              math.pow(y - centre[1], 2) + \
#                              math.pow(z - centre[2], 2))
#
#             # Only consider points that are farther than (radius - buffer size) away from the centre
#             if dist >= (diameter / 2.0 - buffer_size):
#                 # Only add the point coords if it was not already added
#                 if coord not in uniqure_coords:
#                     uniqure_coords.append(coord)
#                     coordinates.append([i] + list(coord))
#                     # print coord
#
#     # Sort the coordinates clockwise, for instance, for loading along y,
#     # sort based on x and z coordinates - the (x, z) of the centroid of the cube
#     # (hence the minus)
#     if (coordinates is not None):
#         coordinates.sort(key=lambda c: math.atan2(c[(axis_of_loading + 1) % 3 + 1] -
#                                                   centre[(axis_of_loading + 1) % 3],
#                                                   c[(axis_of_loading + 2) % 3 + 1] -
#                                                   centre[(axis_of_loading + 2) % 3]))
#     else:
#         print "No points found!"
#
#     return coordinates

# def processOutputUCS(dir_path, file_extension, output_type,
                     # platen_cells_prop_id, platen_points_prop_id,
                     # diameter, height, specimen_center,
                     # gauge_length, gauge_width, do_strain_gauge_analysis,
                     # axis_of_loading, specimen_shape):
def processOutputUCS(dir_path, list_of_files, diameter, height, specimen_center, platen_cells_prop_id, platen_points_prop_id, gauge_length, gauge_width, do_strain_gauge_analysis, modelextend, output_file_name):
    # By default this is a 2D model and analysis
    print(green_text("Processing Started"))
    start_initial_process = time.time()

    num_dimensions = 2

    # list_of_files = findOutputFiles(dir_path, file_extension, output_type)

    print "Processing UCS output files in:", red_text(dir_path)

    # Open a csv file to "append" the data to

    datafile = open(os.path.join(post_processing, output_file_name), 'w')
    row = csv.writer(datafile, delimiter=',')

    # Open a csv file to "append" the data to
    brokenjoint_file = 'brokenjoints.csv'
    broken_datafile = open(os.path.join(post_processing, brokenjoint_file), 'w')
    brokenjoint_row = csv.writer(broken_datafile, delimiter=',')

    source_file = 'sourcejoints.csv'
    source_datafile = open(os.path.join(post_processing, source_file), 'w')
    source_row = csv.writer(source_datafile, delimiter=',')

    # Data titles for crack type and number of broken elements
    br_row_titles = ['Time Step', "Broken Joint ID", "Cell_1", "Cell_2", "Validation_cell_seismic", "Material ID_1",
                  "Material ID_2", "Group Type", "Crack Type", "Failure Mode", "Broken Joint Area", "event energy",
                  "kinetic energy at failure", "kinetic energy at yielding", "Angle"]
    brokenjoint_row.writerow(br_row_titles)

    # Data titles for crack type and number of seismic post-processing
    source_row_titles = ["Time Step", "Node_0_X", "Node_0_Y", "Node_0_Z", "Node_1_X", "Node_1_Y", "Node_1_Z", "Node_2_X", "Node_2_Y", "Node_2_Z", "Node_3_X", "Node_3_Y", "Node_3_Z", "Event Energy", "Event Time Step", "Failure Mode", "Failure Time Step", "Kinetic Energy at Failuer", "Kinetic Energy at Yielding", "Yielding Time Step", "Center_X", "Center Y", "Center Z"]
    source_row.writerow(source_row_titles)

    # Data titles and data of mechanical properties
    row_titles = ["Time Step (-)", "Platen displacement y (mm)", "Strain y from platen (%)",
                  "Strain x (%)", "Strain y (%)", "Perimeter (mm)", "Lateral strain (%)",
                  "Platen Force (kN)", "Axial stress (MPa)",
                  "Volumetric strain (%)",
                  "Material property ID of platens = " + str(platen_cells_prop_id),
                  "Boundary condition ID of platen = " + str(platen_points_prop_id),
                  "Specimen diameter = " + str(diameter),
                  "Specimen height = " + str(height),
                  "Specimen centre = " + str(specimen_center),
                  "Strain gauge length = " + str(gauge_length),
                  "Strain gauge width = " + str(gauge_width)]

    # Perform strain gauge analysis only if the relevant data was provided as inputs

    if specimen_center is None or not specimen_center or not do_strain_gauge_analysis:
        do_strain_gauge_analysis = False
        # Remove row titles related to strain gauges
        row_titles.remove('Strain x (%)')
        row_titles.remove('Strain y (%)')
        row_titles.remove("Specimen centre = " + str(specimen_center))
        row_titles.remove("Strain gauge length = " + str(gauge_length))
        row_titles.remove("Strain gauge width = " + str(gauge_width))

    row.writerow(row_titles)

    # Read the vtu files using an XMLUnstructuredGridReader
    # pv_grid_reader = XMLUnstructuredGridReader(FileName=list_of_files)

    # Show the grid
    Show(pv_grid_reader)

    # Define a new cell array for the output
    cellArray = vtk.vtkCellArray()
    cellTypes = []

    # Get the Client Side Data Object of the grid reader
    # data_object = pv_grid_reader.GetClientSideObject()
    # Get the output
    # output = data_object.GetOutput()


    # Get the required points and cells
    # cells = output.GetCellData()
    # points = output.GetPointData()

    # Older vtu files used a generic 'property_id' while newer files user more
    # descriptive names. We first try to parse the old name for old files.
    # If this fails, which means we're using a more recent vtu file, then we get the new name/field.
    cell_pr_txt = 'property_id'
    point_pr_txt = 'property_id'
    if cells.GetArray('property_id') is None:
        cell_pr_txt = 'material property ID'
    if points.GetArray('property_id') is None:
        point_pr_txt = 'boundary condition ID'
    # print 'property_id of element 0: ', cells.GetArray('property_id').GetTuple1(0)
    # print 'displacement of point 0:',  points.GetArray('displacement').GetTuple3(0)

    # List of points for the vertical (v) and horizontal (h) strain gauges
    pv = []
    ph = []
    if do_strain_gauge_analysis:
        pv.append([specimen_center[0] + gauge_width / 2.0,
                   specimen_center[1] - gauge_length / 2.0,
                   0.0])
        pv.append([specimen_center[0] - gauge_width / 2.0,
                   specimen_center[1] - gauge_length / 2.0,
                   0.0])
        pv.append([specimen_center[0] + gauge_width / 2.0,
                   specimen_center[1] + gauge_length / 2.0,
                   0.0])
        pv.append([specimen_center[0] - gauge_width / 2.0,
                   specimen_center[1] + gauge_length / 2.0,
                   0.0])
        ph.append([specimen_center[0] + gauge_length / 2.0,
                   specimen_center[1] + gauge_width / 2.0,
                   0.0])
        ph.append([specimen_center[0] + gauge_length / 2.0,
                   specimen_center[1] - gauge_width / 2.0,
                   0.0])
        ph.append([specimen_center[0] - gauge_length / 2.0,
                   specimen_center[1] + gauge_width / 2.0,
                   0.0])
        ph.append([specimen_center[0] - gauge_length / 2.0,
                   specimen_center[1] - gauge_width / 2.0,
                   0.0])
    # Dummy variables for FindCell
    p = [0.0, 0.0, 0.0]
    subId = vtk.mutable(0)
    pcoords = [0.0, 0.0, 0.0]
    w = [0.0, 0.0, 0.0]

    # Get the available time steps
    time_steps = pv_grid_reader.TimestepValues
    # print 'time steps', time_steps

    # Check the current view time
    view = GetActiveView()

    # Initial sample volume
    initial_volume = init_area

    area = init_area
    initial_perimeter = init_area

    # Show progress bar
    print_progress(0, len(list_of_files), prefix='Progress:', suffix='Complete')

    # Loop over time steps and perform calculations
    # This portion needs to go to pool...
    # global alist
    error_array, cluster, failure_modes = [], [], []  # For QAQC Error & Cluster CrackType
    mineral_bound, combos, breakdown, fail_count, crack_types, crack_count = [], [], [], [], [], []
    # alist = range(0, int(max(time_steps)))
    # alist.pop(0)  # Remove first element of the list
    # alist.append(int(end) - 1)  # Add the value "end" to the list
    fail_modes, fail_angles = [], []  # reset every time step to match dfn lines against respective broken joints
    for t in time_steps:
        gc_counter = t
        # global areas
        # Update/render the view
        view.ViewTime = t
        Render()

        ###----------------------------------
        ###      VOLUMETRIC STRAIN ANAYLSIS
        ###----------------------------------

        # boundary_coor, xa, ya, za, sorted_x, sorted_y, sorted_z = [], [], [], [], [], [], []
        coordinates = []

        for j in modelextend:
            xs, ys, zs = output.GetPoint(j)[0], output.GetPoint(j)[1], output.GetPoint(j)[2]
            # print(modelextend, )
            # print (xs, ys, zs)
            coordinates.append([xs, ys, zs])
        coordinates.sort(key=lambda c: (math.atan2(c[1] - specimen_center[1], c[0] - specimen_center[0]),
                                        math.sqrt(
                                            (c[1] - specimen_center[1]) ** 2 + (
                                                    c[0] - specimen_center[0]) ** 2)))
        polygon = Polygon(coordinates)  # create polygon based on Point ID on boundary
        areas = polygon.area  # calculate area of the polygon
        # Volumetric strain calculation for the entire samples
        volume = areas - init_area
        if volume is not None and initial_volume >= 0.0:
            volumetric_strain = (volume - initial_volume) / initial_volume * 100.0
            # print t, volume, volumetric_strain

        ###----------------------------------
        ###      STRAIN GAUGE ANAYLSIS
        ###----------------------------------
        if do_strain_gauge_analysis:
            # Find the cells inside of which strain gauge points are located. We do this only at time step zero.
            if t == 0.0:
                cv = []
                ch = []
                for ps in xrange(0, len(pv)):
                    cv.append(int(output.FindCell(pv[ps], None, 0, 1e-3, subId, pcoords, w)))
                    ch.append(int(output.FindCell(ph[ps], None, 0, 1e-3, subId, pcoords, w)))
                    cellArray.InsertNextCell(output.GetCell(cv[ps]))
                    cellArray.InsertNextCell(output.GetCell(ch[ps]))
                    # SetCells of the output with only the selected analyzed cells
                # output.SetCells(vtk.VTK_TRIANGLE, cellArray)
                # print cv, ch

            strain_v_bottom = [0.0, 0.0, 0.0]
            strain_v_top = [0.0, 0.0, 0.0]
            displacement_y = 0.0
            displacement_x = 0.0
            # Loop over cells of strain gauges
            for ps in xrange(0, len(cv)):
                """ Vertical strain gauge """
                # Get the first point of the cell
                cell = output.GetCell(cv[ps])
                for j in xrange(cell.GetNumberOfPoints()):
                    # Get the displacement tuple for each point
                    point_id = cell.GetPointId(j)
                    displacement = points.GetArray('displacement').GetTuple3(point_id)
                    # Vertical contraction is assumed positive
                    # Bottom cells of vertical strain gauge
                    if ps < 2:
                        displacement_y += displacement[1]
                        # print 'bottom', cv[ps], displacement[1]
                    # Top cells of vertical strain gauge
                    else:
                        displacement_y -= displacement[1]
                        # print 'top', cv[ps], displacement[1]

                """ Horizontal strain gauge """
                # Get the first point of the cell
                cell = output.GetCell(ch[ps])
                for j in xrange(cell.GetNumberOfPoints()):
                    # Get the displacement tuple for each point
                    point_id = cell.GetPointId(j)
                    displacement = points.GetArray('displacement').GetTuple3(point_id)
                    # Horizontal expansion is assumed positive
                    # Right cells of horizontal strain gauge
                    if ps < 2:
                        displacement_x -= displacement[0]
                    # Left cells of horizontal strain gauge
                    else:
                        displacement_x += displacement[0]

            displacement_y = displacement_y / 6.0
            displacement_x = displacement_x / 6.0
            # print 'displacement_x & _y at time step', t, displacement_x, displacement_y

            # Calculate strains in percentage (%)
            strain_x = displacement_x / gauge_length * 100.0
            strain_y = displacement_y / gauge_length * 100.0
            # print 'strain_x & _y at time step', t, strain_x, strain_y


        ###----------------------------------
        ###        PLATEN ANAYLSIS
        ###----------------------------------
        top_platen_cellArray = vtk.vtkCellArray()
        bottom_platen_cellArray = vtk.vtkCellArray()

        # Initialize the array of platen forces and displacements with zero (x,y,z)
        top_platen_force = [0.0, 0.0, 0.0]
        bottom_platen_force = [0.0, 0.0, 0.0]
        avg_top_platen_disp = [0.0, 0.0, 0.0]
        avg_bottom_platen_disp = [0.0, 0.0, 0.0]

        # Loop over all cells and find the platen cells based on their element properties. We do this only at time step zero.
        if t == 0.0:
            top_platen_cells = []
            bottom_platen_cells = []
            top_platen_cell_ids = []
            bottom_platen_cell_ids = []
            for i in xrange(output.GetNumberOfCells()):
                # Get the Cell object
                cell = output.GetCell(i)
                # Get the property_id of the current Cell
                cell_pr = int(cells.GetArray(cell_pr_txt).GetTuple1(i))
                # Check if the property_id is equal to the prescribed element property_id
                if (cell_pr == platen_cells_prop_id):
                    # Get the first point of the cell
                    point_0_id = cell.GetPointId(0)
                    point_0_pr = int(points.GetArray(point_pr_txt).GetTuple1(point_0_id))
                    # Top Platen
                    if (point_0_pr == platen_points_prop_id[0]):
                        top_platen_cells.append(cell)
                        top_platen_cell_ids.append(i)
                        # top_platen_cellArray.InsertNextCell(cell)
                    # Bottom Platen
                    elif (point_0_pr == platen_points_prop_id[1]):
                        bottom_platen_cells.append(cell)
                        bottom_platen_cell_ids.append(i)
                        # bottom_platen_cellArray.InsertNextCell(cell)
        # print 'total no. of platen cells: ', len(top_platen_cells), len(bottom_platen_cells)
        # print bottom_platen_cell_ids

        # Get the average displacements and sum of forces
        top_platen_force, avg_top_platen_disp = calculatePlatenForceAndDisplacement(output,
                                                                                    points,
                                                                                    top_platen_cell_ids)
        bottom_platen_force, avg_bottom_platen_disp = calculatePlatenForceAndDisplacement(output,
                                                                                          points,
                                                                                          bottom_platen_cell_ids)
        # print 'top, bottom (Y)', avg_top_platen_disp[1], avg_bottom_platen_disp[1]

        ###----------------------------------
        ###        CRACK CLUSTERING
        ###----------------------------------

        # /// CHARACTERISING/CLUSTERING THE TYPE OF CRACK TAKING PLACE - INTRACRACK/INTERCRACK/GBC /// #
        # The broken joints are mapped based on the failure mode being threshold between 0 amd 4
        # As they occur in pairs, range is advanced in a step of 2. With i and i+1 being the opposing broken joints.
        # The skip is made in increments of 4 (2D crack rectangle) and 6 (3D crack triangle)
        # Match is made on the coordinates of the brokenjoint element.
        idlist = vtk.vtkIdList()  # Dummy for idList
        dfn_line_list = [] # reset every time step to match dfn lines against respective broken joints

        for i in range(0, output_broken.GetNumberOfCells(), 2):
            # print("CRACK clustering")
            # point_a, point_b, point_c = [], [], []
            if 0 < broken_joint_info.GetArray('failure mode').GetTuple1(i) < 4:
                # Co-relating Broken Joints TO Basic Cell #
                # // THIS IS 2D CALCULATIONS! // #
                ax, ay, az = point_a = output_broken.GetPoint(i * node_skip)
                bx, by, bz = point_b = output_broken.GetPoint((i * node_skip) + node_skip)
                ## ROSE DIAGRAM CALCULATIONS
                angle_deg = (math.degrees(math.atan2(by - ay, bx - ax)))  # Calculate the slope of the line in degrees

                if angle_deg < 0:
                    angle_deg = (angle_deg + 180) # All to Positive
                if angle_deg > 90:
                    angle_deg = (angle_deg - 90) # Rotate to 0 to 90
                else:
                    angle_deg = (angle_deg + 270) # Rotate to 270 to 360
                cal_len = math.sqrt((by - ay) ** 2 + (bx - ax) ** 2)  # Calculate the length of the line
                ## FIND CELL ID USING POINTS (from basic)
                cell_1 = output.FindCell(point_a, None, 0, 1e-4, subId, pcoords, w)
                cell_2 = output.FindCell(point_b, None, 0, 1e-4, subId, pcoords, w)
                if t in alist[:]:
                    fail_angles.append(angle_deg)

                # if node_skip == 4: ## ROSE DIAGRAM CALCULATIONS
                #     angle_deg = (math.degrees(math.atan2(point_b[1]-point_a[1], point_b[0]-point_a[0])))
                #     if angle_deg < 0:
                #         angle_deg = (angle_deg + 180) # All to Positive
                #     if angle_deg > 90:
                #         angle_deg = (angle_deg - 90) # Rotate to 0 to 90
                #     else:
                #         angle_deg = angle_deg + 270 # Rotate to 270 to 360
                #
                #  Check to see if Seismic files are output,
                #  If yes, get teh cell ID assocaited with the event for QAQC
                #  SEISMIC DATA MANIPULATION #

                if list_of_files_seismic:
                    seismic_point = seismic_point_info.GetArray('node 0').GetTuple3(i / 2)
                    seismic_point_coordinates = [seismic_point[0], seismic_point[1], seismic_point[2]]
                    cell_seismic = int(output.FindCell(seismic_point_coordinates, None, 0, 1e-4, subId, pcoords, w))

                # QAQC CHECK #

                if list_of_files_seismic:
                    if cell_1 == cell_2 and cell_seismic not in (cell_1, cell_2):
                        error_array.append((t, i, cell_1, cell_2))
                else:
                    if cell_1 == cell_2:
                        error_array.append((t, i, cell_1, cell_2))

                # POPULATE DATA #
                    # Cluster Failure Mode
                    # Cluster Crack Type
                try:
                    material_id_cell_1 = cells.GetArray('material property ID').GetTuple1(cell_1)
                    material_id_cell_2 = cells.GetArray('material property ID').GetTuple1(cell_2)
                    fail_mode = broken_joint_info.GetArray('failure mode').GetTuple1(i)  # Get Array "FAILURE MODE" of SPECIFIC CELL
                    if t in alist[:]:
                        fail_modes.append(fail_mode)
                except AttributeError:
                    print "Minor Error"
                    material_id_cell_1 = 0
                    material_id_cell_2 = 1
                    # fail_mode = 0

                # Clustering Failure Mode (material ids) - sorted in ascending order to avoid repetition
                if t in alist[:]:
                    cluster.append(tuple(sorted([int(material_id_cell_1), int(material_id_cell_2)])))
                for key, values in material_combinations.items():
                    if values == tuple(sorted([int(material_id_cell_1), int(material_id_cell_2)])):
                        group_type = key

                # Clustering Failure Modes (Mode 1/2/Mix)
                # looks for Failure Mode in Dictionary and returns type ELSE "Mix Mode"
                if t in alist[:]:
                    failure_modes.append(failure_mode.get(fail_mode, "Mix Mode"))

                # Identify the type of Crack
                # Between different material => Intergranular
                # Same material AND on DFN => Intergranular

                if list_of_files_dfn is None or not list_of_files_dfn:
                    crack_type = "intragranular"
                else:
                    if material_id_cell_1 != material_id_cell_2:
                        crack_type = "intergranular"
                    else:
                        # Point List of Broken Joint
                        output_broken.GetCellPoints(i, idlist)
                        broken_joint_idList_1 = ([idlist.GetId(k) for k in range(idlist.GetNumberOfIds())])
                        output_broken.GetCellPoints(i + 1, idlist)
                        broken_joint_idList_2 = ([idlist.GetId(k) for k in range(idlist.GetNumberOfIds())])
                        # Looking for X/Y/Z of Broken Joint in DFN

                        if dfn_line_list:
                            for cor_idList_1, cor_idList_2 in zip(range(0, len(broken_joint_idList_1)), range(0, len(broken_joint_idList_2))):
                                xyz_idList_1 = output_broken.GetPoint(broken_joint_idList_1[cor_idList_1])
                                xyz_idList_2 = output_broken.GetPoint(broken_joint_idList_2[cor_idList_2])
                                if xyz_idList_1 in dfn_line_list and xyz_idList_2 in dfn_line_list:
                                    crack_type = "intergranular"
                                else:
                                    crack_type = "intragranular"
                        else:
                            # Since DFN is not defined. There is no possibility of knowing the location of the crack with respect to the similar grains.
                            crack_type = "intergranular / intragranular"

                if t in alist[:]:
                    crack_types.append(crack_type)

                # Write to File

                if list_of_files_seismic:
                    if cell_seismic == cell_1 or cell_seismic == cell_2:
                        br_row = [t]
                        br_row.extend([i, int(cell_1), int(cell_2), int(cell_seismic), material_id_cell_1, material_id_cell_2, group_type, crack_type, fail_mode, broken_joint_info.GetArray('area').GetTuple1(i), seismic_point_info.GetArray('event energy').GetTuple1(i / 2), seismic_point_info.GetArray('kinetic energy at failure').GetTuple1(i / 2), seismic_point_info.GetArray('kinetic energy at yielding').GetTuple1(i / 2), angle_deg])

                        brokenjoint_row.writerow(br_row)

                        ###
                        x_tot, y_tot, z_tot = 0, 0, 0
                        for n in range(0, node_skip):
                            total_x = seismic_point_info.GetArray(n).GetTuple3(i / (node_skip / 2))[0]
                            total_y = seismic_point_info.GetArray(n).GetTuple3(i / (node_skip / 2))[1]
                            total_z = seismic_point_info.GetArray(n).GetTuple3(i / (node_skip / 2))[2]
                            x_tot += total_x
                            y_tot += total_y
                            z_tot += total_z
                        x_avg = x_tot / node_skip
                        y_avg = y_tot / node_skip
                        z_avg = z_tot / node_skip
                        source_rowData = [t]
                        source_rowData.extend(
                            [seismic_point_info.GetArray('node 0').GetTuple3(i / 2)[0],
                             seismic_point_info.GetArray('node 0').GetTuple3(i / 2)[1],
                             seismic_point_info.GetArray('node 0').GetTuple3(i / 2)[2],
                             seismic_point_info.GetArray('node 1').GetTuple3(i / 2)[0],
                             seismic_point_info.GetArray('node 1').GetTuple3(i / 2)[1],
                             seismic_point_info.GetArray('node 1').GetTuple3(i / 2)[2],
                             seismic_point_info.GetArray('node 2').GetTuple3(i / 2)[0],
                             seismic_point_info.GetArray('node 2').GetTuple3(i / 2)[1],
                             seismic_point_info.GetArray('node 2').GetTuple3(i / 2)[2],
                             seismic_point_info.GetArray('node 3').GetTuple3(i / 2)[0],
                             seismic_point_info.GetArray('node 3').GetTuple3(i / 2)[1],
                             seismic_point_info.GetArray('node 3').GetTuple3(i / 2)[2],
                             seismic_point_info.GetArray('event energy').GetTuple1(i / 2),
                             seismic_point_info.GetArray('event time step').GetTuple1(i / 2),
                             seismic_point_info.GetArray('failure mode').GetTuple1(i / 2),
                             seismic_point_info.GetArray('failure time step').GetTuple1(i / 2),
                             seismic_point_info.GetArray('kinetic energy at failure').GetTuple1(i / 2),
                             seismic_point_info.GetArray('kinetic energy at yielding').GetTuple1(i / 2),
                             seismic_point_info.GetArray('yielding time step').GetTuple1(i / 2),
                             x_avg, y_avg, z_avg])
                        source_row.writerow(source_rowData)
                        ###
                    else:
                        print("ERROR in %d: Check Cell contact at # %s # %s # %d" % (t, cell_1, cell_2, cell_seismic))
                else:
                    # All Seismic info removed as no Seismic files found.
                    if cell_2 != cell_1:
                        rowData = [t]
                        rowData.extend([i, int(cell_1), int(cell_2), "NA", material_id_cell_1, material_id_cell_2, group_type, crack_type, fail_mode, broken_joint_info.GetArray('area').GetTuple1(i), "NA", "NA", "NA", angle_deg])
                        brokenjoint_row.writerow(rowData)
                    else:
                        print("ERROR in %d: Check Cell contact at # %s # %s" % (t, cell_1, cell_2))

        # Create bins for graphical output

        if t in alist[:]:
            breakdown.append(Counter(cluster))
            fail_count.append(Counter(failure_modes))
            crack_count.append(Counter(crack_types))
            #FAILURE MODE! ROSE ANGLE! Assign here!
            if crack_types:
                rose_illustration(fail_angles, fail_modes, t)
                # plot_bar_from_counter(Counter(cluster), breakdown, fail_count, crack_count, fail_modes, fail_angles, t)
                fail_angles, fail_modes = [], []



        ###----------------------------------
        ###        FILE OUTPUT
        ###----------------------------------
        # Convert forces from microN to kN and get the average forces & displacements
        avg_platen_force = [0.0, 0.0, 0.0]
        avg_platen_disp = [0.0, 0.0, 0.0]
        for i in range(0, 3):
            # microN to kN
            avg_platen_force[i] = 0.5 * (abs(top_platen_force[i]) + abs(bottom_platen_force[i])) / 1.0e9
            avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

        rowData = [t]

        axis_of_loading = 1
        # stress in MPa (force in kN & area in mm^2)
        stress = avg_platen_force[axis_of_loading] / width * 1.0e3
        # strain calculated from platen displacements
        strain_from_platen = avg_platen_disp[axis_of_loading] / height * 100.0

        if do_strain_gauge_analysis:
            rowData.extend([avg_platen_disp[axis_of_loading], strain_from_platen, strain_x, strain_y,
                            areas, volume,
                            avg_platen_force[axis_of_loading], stress])
        else:
            rowData.extend([avg_platen_disp[axis_of_loading], strain_from_platen,
                            areas, volume,
                            avg_platen_force[axis_of_loading], stress])

        if volume is not None and initial_volume >= 0.0:
            rowData.extend([volumetric_strain])
        # print "\nRow data:", rowData

        # Append the data to the csv file
        row.writerow(rowData)

        # Update progress bar
        print_progress(int(t) + 1, len(list_of_files), prefix='Progress:', suffix='Complete')

    datafile.close()
    broken_datafile.close()
    source_datafile.close()

    print "\nFinished writing processed Stress/Strain data to:", bold_text(os.path.join(post_processing, 'history.csv'))
    Delete(pv_grid_reader)


    # /// RESULTS VERIFICATION /// #
    # First checks to see if there are any broken joints
    # Terminates if none.
    # Else moves on to clustering
    # /// CLUSTER THE RESULTS /// #
    # return results to counter (creates bins for histogram plotting)
    # displays as a separate window and terminate parent
    # Flag set for the enable/disable display
    # Inserted here to avoid any early errors.

    if not crack_types:
        print("\n\033[1;31;0mNO BROKEN JOINTS TO PROCESS\nTERMINATING\033[0m ")
    else:
        if list_of_files_seismic:
            print "Processing Seismic Clustering"
            seismicclustering(post_processing)
        # if answer2 in ('Y', 'y'):
        import matplotlib.pyplot as plt
        # print Counter(cluster) # Uncomment if you wish to see the histogram results
        if node_skip == 4:
            plot_bar_from_counter(Counter(cluster), breakdown, fail_count, crack_count, rose_mode, rose_angle, "main_graph")
        else:
            plot_bar_from_counter(Counter(cluster), breakdown, fail_count, crack_count, initial_azimuth, initial_dip_angle)
    # plt.show()


    single_graph(post_processing)
    global_removal(output, cells, points, data_object)
    print ("\nProcessing Complete: %s" % bold_text(calc_timer_values(time.time() - start_initial_process)))
    print (red_text("-------------------------------"))

def read_history_file(filename):
    global axial_strain, stress
    axial_strain, stress = [], []
    with open(filename) as history_file:
        reader = csv.DictReader(history_file, delimiter=',')
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        for row in reader:
            axial_strain.append(float(row['strain y from platen (%)']))
            stress.append(float(row['axial stress (mpa)']))
    history_file.close()

def single_graph(sub_dirs):
    global compare_stress

    compare_stress, compare_strain = {}, {}
    fig3 = plt.figure(figsize=[12.8, 7.2])
    names = ["history.csv"]  # Name of file to look for

    filename = os.path.join(sub_dirs, names[0])
    read_history_file(filename)
    # Configure name of file (Name of sub_folder)
    # os.chdir(sub_dirs)
    parent_dir_of_file = os.path.dirname(sub_dirs)
    pdf_name = os.path.basename(parent_dir_of_file)
    # Check is Simulation results OR EXP results
    if 'EXP' in pdf_name:
        plt.plot(axial_strain, stress, linestyle="-.", label=pdf_name)
    else:
        plt.plot(axial_strain, stress, linestyle="-", label=pdf_name)
    # Manipulate Figure
    plt.title("Stress Strain")
    plt.xlabel("Strain (%)")
    plt.xlim(xmin=0)
    plt.ylabel("Stress (MPa)")
    plt.ylim(ymin=0)
    plt.legend(loc=0)

    try:
    # Obtain the Max_Stress in each file and its realtive index in the list
        max_value_index, max_values_stress =  max(enumerate(stress), key=operator.itemgetter(1)) # Enumerate to get index of max value in list
        elastic_modulus = ((stress[int(max_value_index) / 2]) / (axial_strain[int(max_value_index) / 2]) / 10 )
        print("%s - \t Max Stress %10.2f & \t Elastic Modulus %10.2f" % (pdf_name, max_values_stress, elastic_modulus))
    except ZeroDivisionError:
        print red_text("There is a problem in the history.csv file\nThe stress/strain values are extremely low.")

    # compare_stress[pdf_name] = (max_values_stress, elastic_modulus)
    # compare_strain[pdf_name] = axial_strain[int(max_value_index) / 2]

    # Save in the relevant sub-directory directory
    plt.savefig(os.path.join(sub_dirs, "stress_strain.pdf")) # Saves in the location of the sub directory
    plt.savefig(os.path.join(sub_dirs, "stress_strain.svg")) # Saves in the location of the sub directory
    plt.savefig(os.path.join(sub_dirs, "stress_strain.png"), dpi=100) # Saves in the location of the sub directory



# /// GRAPHICAL REPRESENTATION /// #
# ax1 = Histogram chart of the various brokenjoints based on the material contact
# ax2 = Stacked Bar Chart of the brokenjoint events in that time window.
# ax3 = Stacked Bar Chart for Fracture Type
# ax4 = Stacked Bar Chart of Crack Type

# def plot_bar_from_counter(counter, breaks, failures, types, rose_angle, rose_mode, list_azimuth, list_dip, ax1=None):
# Based on failure mode!
def plot_bar_from_counter(counter, breaks, failures, types, list_azimuth, list_dip, graph_name, ax1=None):

    import matplotlib.pyplot as plt
    try:
        from cycler import cycler
    except AttributeError:
        print "pvpython and matplotlib may not be linked properly. Try running \n$ sudo apt-get install paraview-python \nElse rebuild matplotlib and cycler, to try to resolve this."
    fig = plt.figure(figsize=[12.8, 7.2])
    if graph_name == "main_graph":
        if ax1 is None:
            fig = plt.figure(figsize=[12.8, 7.2])
            ax1 = fig.add_subplot(221) # General Plot
            ax2 = fig.add_subplot(222) # Fracture Group No.
            ax3 = fig.add_subplot(224) # Fracture Mode Type
            ax4 = fig.add_subplot(223) # Crack Type

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

        N, break_series = len(breaks), {}

        for key in {key for keys in breaks for key in keys}:
            break_series[key] = [(0 if key not in item else item[key]) for item in breaks]

        if len(break_series.keys()) == 1:
            colors = 2
            print("\n\033[0;32;0mTHERE IS 1 MINERAL PHASE IN THE MODEL\n\033[0m ")
        else:
            colors = len(break_series.keys())
        cmap = plt.get_cmap('viridis') # Refer color_maps Python documentation for more maps!
        ax2.set_prop_cycle(cycler('color', [cmap(1. * (i) / (colors - 1)) for i in range(colors)]))

        bottom, x = numpy.zeros(N), range(N)

        for key in break_series:
            ax2.bar(x, break_series[key], alpha= 0.5, label=key, bottom=bottom)
            bottom += numpy.array(break_series[key])

        ax2.set_xticks(numpy.arange(N))
        ax2.set_xticklabels(numpy.array(alist), rotation=45)
        ax2.set_ylabel("Frequency")
        ax2.set_xlabel("Time Step")
        ax2.legend()

        ''' AX3 '''

        colors_ax3 = len(failure_mode.keys()) + 1
        ax3.set_prop_cycle(cycler('color', [cmap(1. * (i) / (colors_ax3 - 1)) for i in range(colors_ax3)]))
        M, fail_series = len(failures), {}
        for key in {key for keys in failures for key in keys}:
            fail_series[key] = [(0 if key not in item else item[key]) for item in failures]
        bottom, x = numpy.zeros(M), range(M)

        for key in fail_series:
            ax3.bar(x, fail_series[key], alpha= 0.5, label=key, bottom=bottom)
            bottom += numpy.array(fail_series[key])

        ax3.set_xticks(numpy.arange(M))
        ax3.set_xticklabels(numpy.array(alist), rotation=45)
        ax3.set_ylabel("Frequency")
        ax3.set_xlabel("Time Step")
        ax3.legend()

        ''' AX4 '''

        colors_ax4 = len(crack_dir.keys())
        ax4.set_prop_cycle(cycler('color', [cmap(1. * (i) / (colors_ax4 - 1)) for i in range(colors_ax4)]))
        P, fail_type_series = len(types), {}
        for key in {key for keys in types for key in keys}:
            fail_type_series[key] = [(0 if key not in item else item[key]) for item in types]
        bottom, x = numpy.zeros(P), range(P)

        for key in fail_type_series:
            ax4.bar(x, fail_type_series[key], alpha= 0.5, label=key, bottom=bottom)
            bottom += numpy.array(fail_type_series[key])

        ax4.set_xticks(numpy.arange(P))
        ax4.set_xticklabels(numpy.array(alist), rotation=45)
        ax4.set_ylabel("Frequency")
        ax4.set_xlabel("Time Step")
        ax4.legend()

        plt.suptitle("Time Window %s to %s" % (begin, end))
        plt.savefig(os.path.join(post_processing, "histogram.png"), dpi = 600)

    ''' FIGURE 2'''
    if skip == 4:
        ax = WindroseAxes.from_ax()
        viridis = plt.get_cmap('viridis')
        ax.bar(list_dip, list_azimuth, normed=False, opening=0.8, edgecolor='white', nsector=36, cmap=viridis, bins=3)
        ax.set_legend()
        if graph_name == "main_graph":
            plt.savefig(os.path.join(post_processing, "2D_final_damage_rosette.png"), dpi=600)
        else:
            sub_post_processing = os.path.join(post_processing, "additional")
            if not os.path.exists(sub_post_processing):  # Check to see if the folder exists
                os.makedirs(sub_post_processing)  # if not then makes the folder
            plt.savefig(os.path.join(sub_post_processing, (str(graph_name) + "2D_final_damage_rosette.png")), dpi=600)

    else:
        import matplotlib.pyplot as plt
        import mplstereonet
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='stereonet')
        # ax.plane(azimuth, angle_rad, c='b', label='Bedding %03d/%02d' % (azimuth, angle_rad))
        # print list_azimuth
        # print list_dip
        # print len(list_azimuth), len(list_dip)
        list_azimuth = [x + 180 for x in list_azimuth]
        strikes, dips = (list_azimuth, list_dip)
        strike, dip = mplstereonet.fit_girdle(strikes, dips)
        cax = ax.density_contourf(strikes, dips, measurement='poles')
        ax.plane(list_azimuth, list_dip, c='k', label='Planes')
        ax.pole(strikes, dips, c='b')
        ax.grid()
        ax.legend()
        # plt.show()

    # return ax1, ax2, ax3, ax4, ax
# /// SEISMIC CLUSTERING /// #
# Identifies and sorts similar events by the values of the event energy generated.
# Returns the time step start and end of all the events identified during that time window.

def seismicclustering(post_processing):
    skip = 4
    global rose_angle, rose_mode, initial_azimuth, initial_dip_angle
    reader = csv.reader(open(os.path.join(post_processing, "brokenjoints.csv")), delimiter=",")
    if skip == 4:
        for TimeStep, BrokenJointID, Cell1, Cell2, Validation_cell_seismic, MaterialID_1, MaterialID_2, GroupType, CrackType, FailureMode, BrokenJoint, eventenergy, kineticenergyatfailure, kineticenergyatyielding, Angle in reader:
            seismicsortedlist = sorted(reader, key=operator.itemgetter(11), reverse=True)
    else:
        for TimeStep, BrokenJointID, Cell1, Cell2, Validation_cell_seismic, MaterialID_1, MaterialID_2, GroupType, CrackType, FailureMode, BrokenJoint, eventenergy, kineticenergyatfailure, kineticenergyatyielding, dip_angle, azimuth in reader:
            seismicsortedlist = sorted(reader, key=operator.itemgetter(11), reverse=True)
    s_clust, break_angle = [], []
    rose_mode, rose_angle = [], []
    azim, initial_dip_angle, initial_azimuth = [], [], []
    mx, mn, b_angle, azim_azim = 0 , 0, 0, 0
    if not seismicsortedlist:
        print("\n\033[1;31;0mNO BROKEN JOINTS TO PROCESS\n\033[0m ")
    else:
        if skip == 4:
            with open(os.path.join(post_processing, "seismicclustering.csv"), 'w') as ostr:
                row = csv.writer(ostr, delimiter=',')
                row_titles = ["End Step", "Start Step", "Cell1", "Cell2", "Group Type", "CrackType", "Failure Mode",
                              "Event Energy", "Angle"]
                row.writerow(row_titles)
                for i in range(len(seismicsortedlist) - 1):
                    if seismicsortedlist[i + 1][11] == seismicsortedlist[i][11]:
                        s_clust.append(float(seismicsortedlist[i][0]))
                        s_clust.append(float(seismicsortedlist[i + 1][0]))
                        break_angle.append((seismicsortedlist[i][14]))
                    else:
                        if len(s_clust) < 1:
                            s_clust.append(float(seismicsortedlist[i][0]))
                            break_angle.append(float(seismicsortedlist[i][14]))
                            mx, mn = s_clust[0], s_clust[0]
                            b_angle = break_angle[0]
                        else:
                            mx, mn = max(s_clust), min(s_clust)
                            b_angle = break_angle[0]
                        row.writerow(
                            [mx, mn, seismicsortedlist[i][2], seismicsortedlist[i][3], seismicsortedlist[i][7],
                             seismicsortedlist[i][8], seismicsortedlist[i][9], seismicsortedlist[i][11], b_angle])
                        s_clust, break_angle = [], []
                        rose_mode.append(float(seismicsortedlist[i][9]))
                        rose_angle.append(float(b_angle))
                row.writerow(
                    [mx, mn, seismicsortedlist[i][2], seismicsortedlist[i][3], seismicsortedlist[i][7],
                     seismicsortedlist[i][8], seismicsortedlist[i][9], seismicsortedlist[i][11], seismicsortedlist[i][14]])
                rose_mode.append(float(seismicsortedlist[i][9]))
                rose_angle.append(float(seismicsortedlist[i][14]))
                ostr.close()
        else:
            # if skip == 6:
            with open(os.path.join(post_processing, "seismicclustering.csv"), 'w') as ostr:
                row = csv.writer(ostr, delimiter=',')
                row_titles = ["End Step", "Start Step", "Cell1", "Cell2", "Group Type", "CrackType", "Failure Mode",
                              "Event Energy", "Dip Angle", "Azimuth"]
                row.writerow(row_titles)
                for i in range(len(seismicsortedlist) - 1):
                    if seismicsortedlist[i + 1][11] == seismicsortedlist[i][11]:
                        s_clust.append(float(seismicsortedlist[i][0]))
                        s_clust.append(float(seismicsortedlist[i + 1][0]))
                        break_angle.append((seismicsortedlist[i][14]))
                        azim.append(seismicsortedlist[i][15])
                    else:
                        if len(s_clust) < 1:
                            s_clust.append(float(seismicsortedlist[i][0]))
                            break_angle.append(float(seismicsortedlist[i][14]))
                            azim.append(float(seismicsortedlist[i][15]))
                            mx, mn = s_clust[0], s_clust[0]
                            b_angle = break_angle[0]
                            azim_azim = azim[0]
                        else:
                            mx, mn = max(s_clust), min(s_clust)
                            b_angle = break_angle[0]
                            azim_azim = azim[0]
                        row.writerow(
                            [mx, mn, seismicsortedlist[i][2], seismicsortedlist[i][3], seismicsortedlist[i][7],
                             seismicsortedlist[i][8], seismicsortedlist[i][9], seismicsortedlist[i][11], b_angle, azim_azim])
                        s_clust, break_angle, azim = [], [], []
                        rose_mode.append(float(seismicsortedlist[i][9]))
                        initial_dip_angle.append(float(b_angle))
                        initial_azimuth.append(float(azim_azim))
                row.writerow(
                    [mx, mn, seismicsortedlist[i][2], seismicsortedlist[i][3], seismicsortedlist[i][7],
                     seismicsortedlist[i][8], seismicsortedlist[i][9], seismicsortedlist[i][11],
                     seismicsortedlist[i][14], seismicsortedlist[i][15]])
                rose_mode.append(float(seismicsortedlist[i][9]))
                initial_dip_angle.append(float(seismicsortedlist[i][14]))
                initial_azimuth.append(float(seismicsortedlist[i][15]))
                ostr.close()

    source_reader = csv.reader(open(os.path.join(post_processing, "sourcejoints.csv")), delimiter=",")
    if skip == 4:
        for TimeStep, node0_x, node0_y, node0_z, node1_x, node1_y, node1_z, node2_x, node2_y, node2_z, node3_x, node3_y, node3_z, event_energy, event_time_step, failure_mode, failure_time_step, kinetic_energy_at_failure, kinetic_energy_at_yielding, yielding_time_step, x_avg, y_avg, z_avg in source_reader:
            source_seismicsortedlist = sorted(source_reader, key=operator.itemgetter(13), reverse=True)
        s_clust = []
        with open(os.path.join(post_processing, "source_mat.csv"), 'w') as ostr:
            source_row = csv.writer(ostr, delimiter=',')
            row_titles = ["Time Step", "X", "Y", "Z", "X", "Y", "Z", "X", "Y", "Z", "X", "Y", "Z","event energy", "event time step", "Failure Mode", "Failure Time Step", "Kinetic Energy - Release", "Kinetic Energy - Yielding", "Yielding Time Step", "X", "Y", "Z"]
            source_row.writerow(row_titles)
            for i in range(len(source_seismicsortedlist) - 1):
                if source_seismicsortedlist[i + 1][13] == source_seismicsortedlist[i][13]:
                    s_clust.append(float(source_seismicsortedlist[i][0]))
                    s_clust.append(float(source_seismicsortedlist[i + 1][0]))
                else:
                    if len(s_clust) < 1:
                        s_clust.append(float(source_seismicsortedlist[i][0]))
                    source_row.writerow([
                        source_seismicsortedlist[i][0], source_seismicsortedlist[i][1], source_seismicsortedlist[i][2],
                        source_seismicsortedlist[i][3], source_seismicsortedlist[i][4], source_seismicsortedlist[i][5],
                        source_seismicsortedlist[i][6], source_seismicsortedlist[i][7], source_seismicsortedlist[i][8],
                        source_seismicsortedlist[i][9], source_seismicsortedlist[i][10],
                        source_seismicsortedlist[i][11], source_seismicsortedlist[i][12],
                        source_seismicsortedlist[i][13], source_seismicsortedlist[i][14],
                        source_seismicsortedlist[i][15], source_seismicsortedlist[i][16],
                        source_seismicsortedlist[i][17], source_seismicsortedlist[i][18],
                        source_seismicsortedlist[i][19], source_seismicsortedlist[i][20],
                        source_seismicsortedlist[i][21], source_seismicsortedlist[i][22]])
                    s_clust = []
            source_row.writerow([
                source_seismicsortedlist[i][0], source_seismicsortedlist[i][1], source_seismicsortedlist[i][2],
                source_seismicsortedlist[i][3], source_seismicsortedlist[i][4], source_seismicsortedlist[i][5],
                source_seismicsortedlist[i][6], source_seismicsortedlist[i][7], source_seismicsortedlist[i][8],
                source_seismicsortedlist[i][9], source_seismicsortedlist[i][10],
                source_seismicsortedlist[i][11], source_seismicsortedlist[i][12],
                source_seismicsortedlist[i][13], source_seismicsortedlist[i][14],
                source_seismicsortedlist[i][15], source_seismicsortedlist[i][16],
                source_seismicsortedlist[i][17], source_seismicsortedlist[i][18],
                source_seismicsortedlist[i][19], source_seismicsortedlist[i][20],
                source_seismicsortedlist[i][21], source_seismicsortedlist[i][22]])
            ostr.close()
    else:
        for TimeStep, node0_x, node0_y, node0_z, node1_x, node1_y, node1_z, node2_x, node2_y, node2_z, node3_x, node3_y, node3_z, node4_x, node4_y, node4_z, node5_x, node5_y, node5_z,event_energy, event_time_step, failure_mode, failure_time_step, kinetic_energy_at_failure, kinetic_energy_at_yielding, yielding_time_step, x_avg, y_avg, z_avg in source_reader:
            source_seismicsortedlist = sorted(source_reader, key=operator.itemgetter(19), reverse=True)
        s_clust = []
        with open(os.path.join(post_processing, "source_mat.csv"), 'w') as ostr:
            source_row = csv.writer(ostr, delimiter=',')
            row_titles = ["Time Step", "X", "Y", "Z", "X", "Y", "Z", "X", "Y", "Z", "X", "Y", "Z", "X", "Y", "Z", "X", "Y", "Z", "event energy", "event time step", "Failure Mode", "Failure Time Step", "Kinetic Energy - Release", "Kinetic Energy - Yielding", "Yielding Time Step", "X", "Y", "Z"]
            source_row.writerow(row_titles)
            for i in range(len(source_seismicsortedlist) - 1):
                if source_seismicsortedlist[i + 1][13] == source_seismicsortedlist[i][13]:
                    s_clust.append(float(source_seismicsortedlist[i][0]))
                    s_clust.append(float(source_seismicsortedlist[i + 1][0]))
                else:
                    if len(s_clust) < 1:
                        s_clust.append(int(source_seismicsortedlist[i][0]))
                    source_row.writerow([
                        source_seismicsortedlist[i][0], source_seismicsortedlist[i][1], source_seismicsortedlist[i][2],
                        source_seismicsortedlist[i][3], source_seismicsortedlist[i][4], source_seismicsortedlist[i][5],
                        source_seismicsortedlist[i][6], source_seismicsortedlist[i][7], source_seismicsortedlist[i][8],
                        source_seismicsortedlist[i][9], source_seismicsortedlist[i][10],
                        source_seismicsortedlist[i][11], source_seismicsortedlist[i][12],
                        source_seismicsortedlist[i][13], source_seismicsortedlist[i][14],
                        source_seismicsortedlist[i][15], source_seismicsortedlist[i][16],
                        source_seismicsortedlist[i][17], source_seismicsortedlist[i][18],
                        source_seismicsortedlist[i][19], source_seismicsortedlist[i][20],
                        source_seismicsortedlist[i][21], source_seismicsortedlist[i][22],source_seismicsortedlist[i][23], source_seismicsortedlist[i][24],source_seismicsortedlist[i][25], source_seismicsortedlist[i][26],source_seismicsortedlist[i][27], source_seismicsortedlist[i][28]])
                    s_clust = []
            source_row.writerow([
                source_seismicsortedlist[i][0], source_seismicsortedlist[i][1], source_seismicsortedlist[i][2],
                source_seismicsortedlist[i][3], source_seismicsortedlist[i][4], source_seismicsortedlist[i][5],
                source_seismicsortedlist[i][6], source_seismicsortedlist[i][7], source_seismicsortedlist[i][8],
                source_seismicsortedlist[i][9], source_seismicsortedlist[i][10],
                source_seismicsortedlist[i][11], source_seismicsortedlist[i][12],
                source_seismicsortedlist[i][13], source_seismicsortedlist[i][14],
                source_seismicsortedlist[i][15], source_seismicsortedlist[i][16],
                source_seismicsortedlist[i][17], source_seismicsortedlist[i][18],
                source_seismicsortedlist[i][19], source_seismicsortedlist[i][20],
                source_seismicsortedlist[i][21], source_seismicsortedlist[i][22],source_seismicsortedlist[i][23], source_seismicsortedlist[i][24],source_seismicsortedlist[i][25], source_seismicsortedlist[i][26],source_seismicsortedlist[i][27], source_seismicsortedlist[i][28]])
            ostr.close()

def global_removal(output, cells, points, data_object):
    del output, cells, points, data_object
    gc.collect()
    Delete(pv_grid_reader)
    print bold_text("Memory Cleared")


'''
Loop over platen cells and calculate the average displacements and sum of forces 
    Inputs:
        - Data object containing the output (main point/cell data)
        - Point data 
        - Cell index (IDs) of platen cells
    Returns:
        - Sum of platen forces (list of 3 floats, [Fx, Fy, Fz])
        - Average of platen displacements (list of 3 floats, [dx, dy, dz])
'''


def calculatePlatenForceAndDisplacement(data_object_output, data_points, platen_cell_ids):
    # Get the number of points per cell (2D: 3; 3D: 4)
    cell_0 = data_object_output.GetCell(0)
    number_of_points_per_cell = float(cell_0.GetNumberOfPoints())

    platen_force = [0.0, 0.0, 0.0]
    platen_disp = [0.0, 0.0, 0.0]
    # Loop over platen cell ids
    for idx in platen_cell_ids:
        # Get the cell for the give idx
        cell = data_object_output.GetCell(idx)
        # Loop over the points of each cell
        for j in xrange(cell.GetNumberOfPoints()):
            # Get the force & displacement tuples for each point
            platen_point_id = cell.GetPointId(j)
            force = data_points.GetArray('force').GetTuple3(platen_point_id)
            displacement = data_points.GetArray('displacement').GetTuple3(platen_point_id)
            for k in range(0, 3):
                platen_force[k] = platen_force[k] + force[k]
                platen_disp[k] = platen_disp[k] + displacement[k]
    # Get the average of displacements
    for k in range(0, 3):
        # divide by the number of points per cell  (3 in 2D and 4 in 3D)
        platen_disp[k] = platen_disp[k] / (number_of_points_per_cell * len(platen_cell_ids))

    return platen_force, platen_disp


# /// Calculate Area & Volume /// #

def area_calculation(output, modelextend, specimen_center):
    # Get the X/Y/Z of the edge co-ordinates
    coordinates = []
    global init_area
    ## CALCULATIONS ##
    for j in modelextend:
        xs, ys, zs = output.GetPoint(j)
        coordinates.append([xs, ys, zs])
    # If 2D Sort points in clockwise direction using Origin as sample center
    # Sort on Angle tan(dy/dx); then on radial distance sqrt(dx^2+dy^2)
    # Calculate Area
    if node_skip == 4:
        coordinates.sort(key=lambda c: (math.atan2(c[1] - specimen_center[1], c[0] - specimen_center[0]),
                                       math.sqrt(
                                           (c[1] - specimen_center[1]) ** 2 + (c[0] - specimen_center[0]) ** 2)) )
        polygon = Polygon(coordinates)  # create polygon based on Point ID on boundary
        init_area = polygon.area  # calculate area of the polygon
    else:
        init_area = ss.ConvexHull(coordinates).volume
    return init_area

# /// Function to enable data collection for Rosette & Stereonets /// #
# /// Can be four conditions
# CONDITION A) 2D with only one type of DFN property (DATA FROM DFN FILE)
# CONDITION B) 2D with multi type of DFN property (DATA FROM BROKEN JOINT FILE AND CHECKED FROM BASIC FILE)
# CONDITION C) 3D with only one type of DFN property (DATA FROM DFN FILE)
# CONDITION D) 3D with multi type of DFN property

def illustration(output, cells, output_broken, broken_joint_info, output_dfn, dfn_joint_info, platen_cells_prop_id, dfn_output_types, out_name):
    print green_text("DFN Found - Processing Initial DFN")
    rose_angle, damage = [], []
    list_dip, list_azimuth = [], []
    if node_skip == 4 and len(dfn_output_types) == 1: # CONDITION A
        print("Condition A - 2D Single DFN")
        for i in range(0, output.GetNumberOfCells(), 2):
            ax, ay, az = point_a = output.GetPoint(i * node_skip) # Get Point A of the line
            bx, by, bz = point_b = output.GetPoint((i * node_skip) + node_skip) # Get Point B of the line
            angle_deg = (math.degrees(math.atan2(by - ay, bx - ax)))  # Calculate the slope of the line in degrees
            cal_len = math.sqrt((by - ay) ** 2 + (bx - ax) ** 2)  # Calculate the length of the line
            rose_angle.append(angle_deg)
            damage.append(cal_len)
        print("Condition A - 2D Single DFN - %s" % output.GetNumberOfCells())
        print("Length of DFN: Max Length %.2f \tMin Length %.2f\tAverage %.2f" % (
        max(damage), min(damage), sum(damage) / len(damage)))
    elif node_skip == 4 and len(dfn_output_types) != 1: # CONDITION B
        print("Condition B")
        for i in range(0, output_broken.GetNumberOfCells(), 2):
            # if broken_joint_info.GetArray('failure mode').GetTuple1(i) == 0:
            ax, ay, az = point_a = output_broken.GetPoint(i * node_skip)  # Get Point A of the line
            bx, by, bz = point_b = output_broken.GetPoint((i * node_skip) + node_skip)  # Get Point B of the line
            if ax not in [sample_x_min, sample_x_max] and bx not in [sample_x_min, sample_x_max] and ay not in [sample_y_min, sample_y_max] and by not in [sample_y_min, sample_y_max]:
                cell_1 = output.FindCell(point_a, None, 0, 1e-4, subId, pcoords, w) # Lookup the cell of the point
                cell_2 = output.FindCell(point_b, None, 0, 1e-4, subId, pcoords, w) # Lookup the cell of the point
                material_id_cell_1 = cells.GetArray('material property ID').GetTuple1(cell_1) # Lookup the material ID of the cell
                material_id_cell_2 = cells.GetArray('material property ID').GetTuple1(cell_2) # Lookup the material ID of the cell
                if material_id_cell_1 != platen_cells_prop_id or material_id_cell_2 != platen_cells_prop_id:
                    angle_deg = (math.degrees(math.atan2(by-ay, bx-ax))) # Calculate the slope of the line in degrees
                    cal_len = math.sqrt((by-ay)**2 + (bx-ax)**2) # Calculate the length of the line
                    rose_angle.append(angle_deg)
                    damage.append(cal_len)
        print("2D Single Multi-DFN - Total No.: %s" % len(rose_angle))
        print("Length of DFN elements:\tMax %.2f\tMin %.2f\tAverage %.2f" % (
        max(damage), min(damage), sum(damage) / len(damage)))
    elif node_skip == 6 and len(dfn_output_types) == 1: # CONDITION C
        print("3D Single DFN - %s" % output_dfn.GetNumberOfCells())
        reference_vector = numpy.array([0, 1, 0])  # reference vector pointing upwards
        north_direction = [1, 0] # North direction
        for i in range(0, output_dfn.GetNumberOfCells(), 2):
            point_a = output_dfn.GetPoint(i * node_skip)  # Get Point A of the Plane
            point_b = output_dfn.GetPoint((i * node_skip) + 1)  # Get Point B of the Plane
            point_c = output_dfn.GetPoint((i * node_skip) + 2)  # Get Point C of the Plane
            sorted_points = [point_a, point_b, point_c]
            point_a, point_b, point_c = sorted_points = sorted(sorted_points, key=itemgetter(2, 1), reverse=True) # Sort points based on Axis
            # Convert points into array for further processing
            p1 = numpy.array(point_a)
            p2 = numpy.array(point_b)
            p3 = numpy.array(point_c)
            # These two vectors are in the plane
            v_1 = p3 - p1
            v_2 = p2 - p1

            # The cross product is a vector normal to the plane
            # vector equation of the plane
            cp_x, cp_y, cp_z = cp = numpy.cross(v_1, v_2)

            # Returns angle between the cross-product and reference vector.
            if angle_between(cp, reference_vector) > 90:
                dip_angle = (180 - angle_between(cp, reference_vector))
            else:
                dip_angle = angle_between(cp, reference_vector)
            list_dip.append(dip_angle)
            unit_a, unit_b, unit_c = unit_vector(cp)  # returns the values of the unit vectors
            if unit_b > 0:
                unit_a, unit_b, unit_c = cp = [x * -1 for x in cp]

            # Azimuth is the angle between the plane and the unit vector.
            azimuth = angle_between([unit_a, unit_c], north_direction) - 90
            list_azimuth.append(azimuth)

    elif node_skip == 6 and len(dfn_output_types) != 1:  # CONDITION D
        print("3D Single Multi-DFN - %s" % output_broken.GetNumberOfCells())
        exit("UNKNOWN CONDITION")

    if node_skip == 4:
        count_initial_dfn = len(rose_angle)
        rose_illustration(rose_angle, damage, "2D_initial_damage_intensity_rosette")

    # else:
    #     count_initial_dfn = len(list_azimuth)
    #     azimu_illustration(list_azimuth, list_dip)
    #     plt.suptitle("Density coutour of the Poles\nNorth Vector %s - Elevation %s" % ((reference_vector), (north_direction)),
    #                  fontsize=16)
    #     plt.savefig(os.path.join(post_processing, out_name), dpi=600)

# /// Function to return normal unit vectors /// #

def unit_vector(vector):
    return vector / numpy.linalg.norm(vector)

# /// Function to return angles between vectors /// #

def angle_between(v1, v2):
    # Calculate the unit vector of each vector
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    # Returns the cos angle in degrees
    return math.degrees(numpy.arccos(numpy.clip(numpy.dot(v1_u, v2_u), -1.0, 1.0)))

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
    return min_range, max_range,inter

# /// Function for 2D Rose Diagrams /// #

def rose_illustration(rose_angle, damage, fname, ax=None):
    illus_var(damage)
    ax = WindroseAxes.from_ax()
    viridis = plt.get_cmap('viridis')
    if fname in ("2D_mesh_rosette", "2D_initial_damage_intensity_rosette"):
        ax.bar(rose_angle, damage, normed=True, opening=0.8, edgecolor='white', nsector=36, cmap=viridis,
               bins=numpy.arange(min_range, max_range, inter))
        ax.set_xticklabels(['90', '', '0', '', '', '', '0', ''])
        ax.set_theta_zero_location("N")
        ax.set_legend(title="Length", loc="upper left")
    else:
        ax.bar(rose_angle, damage, normed=True, opening=0.8, edgecolor='white', nsector=36, cmap=viridis,
               bins=numpy.arange(1.0, 2.5, 0.5))
        ax.set_xticklabels(['0', '', '90', '', '0', '', '', ''])
        ax.set_theta_zero_location("E")
        ax.set_legend(title="Failure Mode", loc="upper left")
    # ax.bar(rose_angle, damage, normed=True, opening=0.8, edgecolor='white', nsector=36, cmap=viridis,
    #        bins=numpy.arange(min_range, max_range, inter))
    # ax.set_xticklabels(['90', '', '0', '', '', '', '0', ''])
    # ax.set_theta_zero_location("N")
    # ax.set_legend(title="Length", loc="upper left")
    if fname == "2D_initial_damage_intensity_rosette":
        plt.suptitle("Fracture Intensity $P_{21}$ %.2f $mm^{-1}$" % (sum(damage) / (width * height)),
                     fontsize=16)
        plt.savefig(os.path.join(post_processing, fname), dpi=600)
    elif fname ==  "2D_mesh_rosette":
        plt.suptitle("Mesh Plot - %s lines" % len(new_overall_points), fontsize=16)
        plt.savefig(os.path.join(post_processing, fname), dpi=600)
    else:
        sub_post_processing = os.path.join(post_processing, "additional")
        plt.suptitle("Time Step %s - %s cracks" % (int(fname), len(damage)), fontsize=16)
        if not os.path.exists(sub_post_processing):  # Check to see if the folder exists
            os.makedirs(sub_post_processing)  # if not then makes the folder
        plt.savefig(os.path.join(sub_post_processing, ("Time_Step_" + str("%0.4d" % int(fname)) + "_2D_final_damage_rosette.png")), dpi=600)
    plt.close('all')


# /// Function to validate User Input /// #

def validate(int_val):
    try:
        int(int_val)
    except ValueError:
        print("ERROR. Please input Integer Value")
        return False

# /// Function to capture User Input /// #

def user_data():
    global var_matplotlib, begin, end, answer1, answer2, alist, bin_freq, list_of_files_seismic
    while True:
        answer = raw_input("Do you wish to enter the time window manually? \033[1m[Y/N]\033[0m ")
        if str(answer) not in ('N', 'n', 'y', 'Y'):
            print("ERROR. Please Type \033[1m\'Y\' or \'N\'\033[0m")
        else:
            break

    if answer in ('Y', 'y'):
        while True:
            begin = (raw_input("Start Frame: "))
            while not validate(begin):
                begin = (raw_input("Start Frame: "))
            if int(begin) > int(max(time_steps)):
                print("ERROR. Value should be less than \033[1m%d\033[0m" % int(max(time_steps)))
            else:
                break
        while True:
            end = (raw_input("End Frame: "))
            while not validate(end):
                end = (raw_input("End Frame: "))
            if int(end) > int(max(time_steps)):
                print("ERROR. Value should be less than \033[1m%d\033[0m" % int(max(time_steps)))
            elif int(end) <= int(begin):
                print("ERROR. Value should be greater than \033[1m%d\033[0m" % int(begin))
            else:
                break
    elif answer in ('N', 'n'):
        begin = int(0)
        end = int(max(time_steps))

    while True:
        bin_freq = (raw_input("Enter Bin Frequency: "))
        while not validate(bin_freq):
                bin_freq = (raw_input("Enter Bin Frequency: "))
        if int(bin_freq) > (int(end) - int(begin)):
            print("Value of frequency should be less than %d" % (int(end) - int(begin)))
        else:
            break

    while True:
        if list_of_files_seismic:
            answer1 = (raw_input("Do you wish run seismic clustering? \033[1m[Y/N]\033[0m "))
            if answer1 not in ('N', 'n', 'y', 'Y'):
                print("ERROR. Please Type \033[1m\'Y\' or \'N\'\033[0m")
            else:
                break
        else:
            print("\n\033[1;31;0mNO SEISMIC FILES TO PROCESS\n\033[0m ")
            break

    while True:
        if not var_matplotlib:
            answer2 = raw_input("Do you wish to display graphical results? \033[1m[Y/N]\033[0m ")
            if str(answer2) not in ('N', 'n', 'y', 'Y'):
                print("ERROR. Please Type \033[1m\'Y\' or \'N\'\033[0m")
            else:
                break
        else:
            answer2 = "N"
            print("VISUALIZATION DISABLED - %s is missing modules/outdated" % var_matplotlib)
            break

    alist = range(int(begin), int(end), int(bin_freq))
    alist.pop(0) # Remove first element of the list
    alist.append(int(end) - 1) # Add the value "end" to the list

# /// Calculate Irregular Polygon Area /// #

def PolyArea(x,y):
    return 0.5*numpy.abs(numpy.dot(x,numpy.roll(y,1))-numpy.dot(y,numpy.roll(x,1)))

# http://en.wikipedia.org/wiki/Heron's_formula

def area(a, b, c):
    def distance(p1, p2):
        return math.hypot(p1[0]-p2[0], p1[1]-p2[1])

    side_a = distance(a, b)
    side_b = distance(b, c)
    side_c = distance(c, a)
    s = 0.5 * ( side_a + side_b + side_c)
    return math.sqrt(s * (s - side_a) * (s - side_b) * (s - side_c))


'''
Find the perimeter of the polygon defined by a list of sorted coordinates
Perimeter = Sum of the polygon edge lengths (projected on x-z plane, i.e, dy=0)
The perimeters parallel to two planes are also reported individually.
This will help in calculating the Poisson's ratio in two orthogonal directions
    Inputs:
        - Points coordinates: list of [x, y, z]
        - Coordinates of the center of the specimen [x, y, z]
        - Axis of loading: x=0; y=1; z=2
        - An epsilon value 
    Returns:
        - Perimeters along two axes and the overall perimeter
'''


def get_perimeter(coordinates, centre, axis_of_loading, epsilon):
    perimeter = 0.0
    dist_x = 0.0
    dist_z = 0.0

    if coordinates is None or not coordinates:
        print "Coordinate array is empty!"
    else:
        for i in range(len(coordinates)):
            p0 = coordinates[i]
            if (i < len(coordinates) - 1):
                p1 = coordinates[i + 1]
            else:
                p1 = coordinates[0]
            # Project points to the plane perpendicular to axis of loading
            p0[axis_of_loading] = centre[axis_of_loading]
            p1[axis_of_loading] = centre[axis_of_loading]

            distSquared = vtk.vtkMath.Distance2BetweenPoints(p0, p1)
            dist = math.sqrt(distSquared)
            perimeter = perimeter + dist

            if (math.fabs(p1[(axis_of_loading + 1) % 3] -
                          p0[(axis_of_loading + 1) % 3]) < epsilon):  # Along z
                dist_z = dist_z + dist
            elif (math.fabs(p1[(axis_of_loading) % 3] -
                            p0[(axis_of_loading) % 3]) < epsilon):  # Along x
                dist_x = dist_x + dist

    return dist_x, dist_z, perimeter

# /// Main function with command line option handling /// #

def main(argv):
    # Parse arguments from user
    parser = argparse.ArgumentParser()
    #3D - BD - No Flaws
    # parser.add_argument('-d', '--dir', type=str, default="/home/andrea/Paraview_output_files/BD_Flow_3D-R_11",
    #                     help="Base directory containing ParaView output files")
    # 3D - UCS - With Flaws
    # parser.add_argument('-d', '--dir', type=str, default="/hdd/OUTPUTS/ucs_with_3_flaws",
    #                     help="Base directory containing ParaView output files")
    # 2D - UCS - With Flaws
    parser.add_argument('-d', '--dir', type=str, default="/hdd/OUTPUTS/damaged_UCS/90Dam_25Red",
                        help="Base directory containing ParaView output files")
    # 2D - UCS - No Flaws
    # parser.add_argument('-d', '--dir', type=str, default="/hdd/OUTPUTS/Verification_Model")
    parser.add_argument('-b', '--batch', action="store_true",
                        help="Batch process of output files in the subdirectories of the given directory")
    parser.add_argument('-c', '--center', nargs=3, type=float, dest='specimen_center',
                        help="Coordinates of the center of the specimen: x y z")
    parser.add_argument('-e', '--element', type=int, dest='platen_cells_prop_id', default=3,
                        help="Element material property id of the platen cells")
    parser.add_argument('-p', '--point', nargs=2, type=int, dest='platen_points_prop_id',
                        help="Nodal boundary condition ids of the platen points [top platen, bottom platen]")
    parser.add_argument('--gl', type=float, default=10.0,
                        help="Strain gauge length in mm")
    parser.add_argument('--gw', type=float, default=3.0,
                        help="Strain gauge width in mm")
    parser.add_argument('-o', '--output', dest='output_file_name', default="history.csv",
                        help="Name of the output csv file")
    args = parser.parse_args()

    # If '.' is specified, the current directory is used as the directory path
    if args.dir == '.':
        args.dir = os.getcwd()

    # Placeholder for the list of directory(ies). Change relative path to absolute path.
    directories = [os.path.abspath(args.dir)]

    # If batch processing, we'll need to do the analysis on all subdirectories of the given directory
    if args.batch:
        print('\nBatch processing...\n')

        directories = findSubdirectories(os.path.abspath(args.dir))

    # Loop over directory(ies)
    if directories is not None and directories:
        for sub_dir in directories:
            processOutputPropID(sub_dir, "vtu", args.specimen_center, args.platen_cells_prop_id, args.platen_points_prop_id, args.output_file_name)


if __name__ == "__main__":
    try:
        # from multiprocessing import Pool
        # pool = Pool(multiprocessing.cpu_count())  # Create a multiprocessing Pool
        # q = sys.argv[1:]
        #
        # data_outputs = pool.map(main, q)

        # for process_idx in range(multiprocessing.cpu_count() / 2):

        q = sys.argv[1:]
        p = multiprocessing.Process(target=main, args=(q,))

        # print(process_idx % multiprocessing.cpu_count())
        p.start()
        # p.join()
        # p = psutil.Process(sys.settrace(main(sys.argv[1:])))
        # p = (target=sys.settrace(main), args=(sys.argv[1:]))

        # reset affinity against all CPUs
        # all_cpus = list(range(psutil.cpu_count()))
        # p.cpu_affinity(all_cpus)
        # q = sys.argv[1:]
        # pool = multiprocessing.Pool()
        # p = multiprocessing.Process(target=main, args=(q,))
        # all_cpus = list(range(psutil.cpu_count()))
        # # p.cpu_affinity(all_cpus)
        # p.start()
        # p.join()
        # print '\n\nfinished!\n\n'
    except KeyboardInterrupt:
        # print("\n\033[1;31;0mTERMINATED BY USER\n")
        exit("TERMINATED BY USER")