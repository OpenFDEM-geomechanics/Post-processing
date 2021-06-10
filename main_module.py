# /////////////////////////////////////////////////////////////// #
# !python2.7
# -*- coding: utf-8 -*-
# Python Script initially created on 15/03/17
# Compiled by Aly @ Grasselli's Geomechanics Group, UofT, 2017
# Base Script by Omid Mahabadi (Courtesy of Geomechanica)
# B value code courtesy of Qi Zhao
# Created using PyCharm
# /////////////////////////////////////////////////////////////// #

global reset_scale
reset_scale = 0

try:
    paraview.simple
    print('Just paraview.simple')
except:
    from paraview.simple import *
    print('paraview.simple import *')
import itertools, operator, csv, argparse, subprocess
import numpy, os, re, sys, time, math, gc
import platform
import paraview.vtk as vtk
from collections import Counter
import json

import formatting_codes
import rose_Illustrations
import geometry_calculator
import graphical_illustrations
import seismic_illustrations
import process_b_value

# Displays the pid on the system
# Check paraview Version before proceeding
try:
    if platform.system()=="Linux":
        print("Computer Running %s, Distribution %s, Release %s" % (platform.system(), platform.linux_distribution()[0], platform.linux_distribution()[1]))
    else:
        print("Computer Running %s, Release %s, Version %s" % (platform.system(), platform.release(), platform.version()))
except AssertionError:
    exit("Incompatible Linux Distribution 16.0+")

try:
    assert sys.version_info >= (2, 6)
    print("Python Version: %s" % sys.version.split('\n')[0])
except AssertionError:
    exit("Incompatible Python Version 2.7.6+")

try:
    if GetParaViewVersion() > 3:
        print("Using Paraview Source Version %s and Paraview Version %s." % (GetParaViewSourceVersion(), GetParaViewVersion()))
    else:
        print("Script may not be compatible")
        exit("Requires Paraview 3+")
except NameError:
    exit("Requires Paraview 4.2+")


''' DICTIONARIES '''

crack_dir = {"intergranular": 1, "intragranular": 2, "intergranular / intragranular": 3}
failure_mode = {1: "Mode 1", 2: "Mode 2", 3: "Simultaneous"}


''' 
Sort filenames based on natural order (e.g. 1, 2,..., 10, 11, ..., instead of 1, 10, 11, 2, ...) 
    Return:
        - Sorted list
'''

numbers = re.compile(r'(\d+)')

# Strip the numerical portion of the file
def numericalSort(value):
    parts = numbers.split(value)  # Split the numerical part of the file
    parts[1::2] = map(int, parts[1::2])  # Return the numerical portion of the file
    return parts


'''
Find ParaView output files in a given directory
    Inputs: 
        - Directory where output files are located (full path is required)
        - File extension of output files (usually '.vtu')
        - Type of output files, for instance use '_basic_' for the basic/main output files
    Returns:
        - Sort them based on their numerical value
        - List of filenames (complete path)
'''


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
Find all the subdirectories of a given path
    Inputs:
        - Starting directory (full path is required)
    Returns:
        - A list of subdirectories
'''


def findSubdirectories(dir_path):
    sub_dirs = []
    listtoremove = []
    for root, dirs, files in os.walk(dir_path):
        for dir_name in dirs:

            sub_dirs.append(os.path.join(root,dir_name))
    sub_dirs = sorted(list(set(sub_dirs))) # Sort directories alphabetically in ascending order
    for i in sub_dirs:
        if 'post_processing' in i:
            listtoremove.append(i)
            listtoremove.append(i.rsplit('/',1)[0])
            # print(i)
            # print(i.rsplit('/',1)[0])
    print("Found \033[1m%s\033[0m sub-directories" % len(sub_dirs))
    sub_dirs = [x for x in sub_dirs if x not in listtoremove]
    print("After remove already processed file\nFound \033[1m%s\033[0m sub-directories" % len(sub_dirs))
    return sub_dirs


'''
Create post_processing folder
Check to see if the folder exists, else creates it
    Inputs:
        - Main Folder
        - Name of postprocessing folder for outputs
    Returns:
        - Name of postprocessing directory (complete path)
'''


def create_post_processing(dir_path, sub_dir_name):
    global post_processing
    post_processing = os.path.join((dir_path), sub_dir_name)
    if not os.path.exists(post_processing):  # Check to see if the folder exists
        os.makedirs(post_processing)  # if not then makes the folder
    return post_processing


''' 
Print iterations progress
Source: https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
'''


def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=50):
    # Adjusted bar length to 50, to display on small screen
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


'''
Process ParaView output of a 2D/3D UCS simulation. Writes the data to a CSV file.
    Inputs:
        - Directory where output files are located
        - File extension of output files (usually '.vtu')
        - Type of output files, for instance use '_basic_' for the basic/main output files
        - Element (cell) property_id of the platens
        - Point (node or boundary condition) property_id of platen nodes as a list: top platen, bottom platen
        - Diameter (width) of the sample in mm
        - Height of the sample in mm
        - Coordinates of the center of the specimen [x, y, z]
        - Strain gauge length (mm)
        - Strain gauge width (mm)
        - Flag to indicate if strain gauge analysis should be done, given relevant data is provided
          (the flag could become False if this data is not available)
        - Axis of loading: x=0; y=1; z=2 
        - Shape of the specimen cross-section (circle or square)
'''


# /// PROCESS PARAVIEW OUTPUTS /// #
    # Inputs:
        # Directory where output files are located
        # Looks for the various file extensions to decide on processing

def processOutputPropID(dir_path, file_extension, platen_cells_prop_id, platen_points_prop_id, output_file_name, bin_frequency, gauge_length=0, gauge_width=0):
    start_prop = time.time() # Processing Start Time

    print("Processing data files in:" + formatting_codes.red_text(dir_path))

    ## LOAD FILES ##
    #  List BASIC Files => Fatal Error (Move to next folder)
    global list_of_files
    list_of_files = findOutputFiles(dir_path, file_extension, "_basic_")
    #  If no output files were found, warn the user and abort
    if list_of_files is None or not list_of_files:
        print(formatting_codes.bold_text("\nERROR: ") +
            "No output files matching extension '%s' containing '%s' were found in %s\n%s\n" % (file_extension, "_basic_", dir_path, formatting_codes.red_text("ABORTING")))
        return

    # List BROKEN JOINTS Files => Fatal Error (Move to next folder)
    global list_of_files_broken
    list_of_files_broken = findOutputFiles(dir_path, file_extension, "_broken_joint_")
    #  If no output files were found, warn the user and abort
    if list_of_files_broken is None or not list_of_files_broken:
        print(formatting_codes.bold_text("\nERROR: ") +
            "No output files matching extension '%s' containing '%s' were found in %s\n%s\n" % (file_extension, "_broken_joint_", dir_path, formatting_codes.red_text("ABORTING")))
        return

    # List SEISMIC Files => Warning Error
    global list_of_files_seismic
    list_of_files_seismic = findOutputFiles(dir_path, file_extension, "_seismic_event_")
    #  If no output files were found, warn the user and continue
    if list_of_files_seismic is None or not list_of_files_seismic:
        print(formatting_codes.bold_text("\nERROR: ") +
            "No output files matching extension '%s' containing '%s' were found in %s\nProceeding\n"
            % (file_extension, "_seismic_event_", dir_path))
        list_of_files_seismic = None

    # List DFN Files => Warning Error
    global list_of_files_dfn
    list_of_files_dfn = findOutputFiles(dir_path, file_extension, "_dfn_property_")
    #  If no output files were found, warn the user and continue
    if list_of_files_dfn is None or not list_of_files_dfn:
        print(formatting_codes.bold_text("\nERROR: ") +
            "No output files matching extension '%s' containing '%s' were found in %s\n\033[1;31;0mWill not identify DFN and/or intragranular cracks\033[0m\nProceeding\n"
            % (file_extension, "_dfn_property_", dir_path))

    print("Found %d output files (=time steps)" % len(list_of_files))
    ## LOAD FILES COMPLETE ##

    # Set global variables
    global subId, pcoords, w, node_skip, time_steps, begin, end, mineral_bound_2, angle_deg, list_azimuth, list_dip, output, output_broken, output_dfn, output_seismic, cells, points, broken_joint_info, data_object, pv_grid_reader, pv_grid_reader_broken, pv_grid_reader_seismic, pv_grid_reader_dfn, seismic_point_info, alist

    # Read the vtu files using an XMLUnstructuredGridReader
    pv_grid_reader = XMLUnstructuredGridReader(FileName=list_of_files)  # Basic File
    pv_grid_reader_broken = XMLUnstructuredGridReader(FileName=list_of_files_broken)  # Broken Joint File
    if list_of_files_seismic:  # Seismic File
        pv_grid_reader_seismic = XMLUnstructuredGridReader(FileName=list_of_files_seismic)
    if list_of_files_dfn:  # DFN File
        pv_grid_reader_dfn = XMLUnstructuredGridReader(FileName=list_of_files_dfn)

    # Show the grid
    Show(pv_grid_reader)  # Basic File
    Show(pv_grid_reader_broken)  # Broken Joint File
    if list_of_files_seismic:  # Seismic File
        Show(pv_grid_reader_seismic)
    if list_of_files_dfn:  # DFN File
        Show(pv_grid_reader_dfn)

    # Define a new cell array for the output
    cellArray = vtk.vtkCellArray()

    # Get the Client Side Data Object of the grid reader
    data_object = pv_grid_reader.GetClientSideObject()  # Basic File
    data_object_broken = pv_grid_reader_broken.GetClientSideObject()  # Broken Joint File
    if list_of_files_seismic:
        data_object_seismic = pv_grid_reader_seismic.GetClientSideObject()
    if list_of_files_dfn:  # DFN File
        data_object_dfn = pv_grid_reader_dfn.GetClientSideObject()

    # Get the output
    output = data_object.GetOutput()  # Basic File
    output_broken = data_object_broken.GetOutput()  # Broken Joint File
    if list_of_files_seismic:  # Seismic File
     output_seismic = data_object_seismic.GetOutput()
    if list_of_files_dfn:  # DFN File
        output_dfn = data_object_dfn.GetOutput()

    # Get the required information of points and cells
    # Basic File => Cell // Point
    cells = output.GetCellData()
    points = output.GetPointData()
    # Broken Joint File (Cell)
    broken_joint_info = output_broken.GetCellData()
    # Seismic Points (Point)
    if list_of_files_seismic:
        seismic_point_info = output_seismic.GetPointData()
    # DFN => Joint (Cell) // Point
    if list_of_files_dfn:
        dfn_point_info = output_dfn.GetPointData()
        dfn_joint_info = output_dfn.GetCellData()

    # Get the available time steps
    time_steps = pv_grid_reader.TimestepValues
    # print 'time steps', time_steps

    # Older vtu files used a generic 'property_id' while newer files user more
    # descriptive names. We first try to parse the old name for old files.
    # If this fails, which means we're using a more recent vtu file, then we get the new name/field.
    cell_pr_txt = 'property_id'
    point_pr_txt = 'property_id'
    if cells.GetArray('property_id') is None:
        cell_pr_txt = 'material property ID'
    if points.GetArray('property_id') is None:
        point_pr_txt = 'boundary condition ID'

    # Check the current view time
    view = GetActiveViewOrCreate('RenderView')
    view.ResetCamera()  # Reset the camera to display the entire model.

    # Check 2D (3 Points - Triangle) or 3D (4 Points - Tetrahedral) Simulations from the cell vertex.
    global node_skip, node_skip
    check_dim_cell = output.GetCell(0)
    number_of_points_per_cell = float(check_dim_cell.GetNumberOfPoints())
    if number_of_points_per_cell == 3: #2D (Triangles)
        print (formatting_codes.green_text("2D Simulation"))
        node_skip = 4
    else: #3D (Tetrahedral)
        print (formatting_codes.green_text("3D Simulation"))
        node_skip = 6
        exit("3D Simulation not supported")

    # Variables
    top_platen_cellArray = vtk.vtkCellArray()
    bottom_platen_cellArray = vtk.vtkCellArray()
    subId = vtk.mutable(0)  # Dummy for FindCell
    pcoords = [0.0, 0.0, 0.0]  # Dummy for FindCell
    w = [0.0, 0.0, 0.0]  # Dummy for FindCell
    idlist = vtk.vtkIdList() # Dummy for idList
    error_array, cluster, failure_modes = [], [], []  # For QAQC Error & Cluster CrackType
    mineral_bound, combos, breakdown, fail_count, crack_types, crack_count = [], [], [], [], [], []
    global material_combinations
    material_combinations ={}
    list_azimuth, list_dip, list_points = [], [], []
    dfn_line_list = []  # reset every time step to match dfn lines against respective broken joints
    cluster, failure_modes, crack_types = [], [], []  # reset for noncumulative graphs
    modelextend, cellextent = [], []
    global top_platen_cell_ids, bottom_platen_cell_ids
    top_platen_cells, top_platen_cell_ids  = [], []
    bottom_platen_cells, bottom_platen_cell_ids = [], []

    # Prompt for user input first
    global begin, end, bin_freq
    if not bin_frequency:
        if int(max(time_steps)) < 11:
            bin_freq = 1
        else:
            bin_freq = int(max(time_steps)) / 10
    else:
        bin_freq = bin_frequency
    begin, end = int(0), int(max(time_steps))
    alist = list(range(int(begin), int(end), int(bin_freq)))
    alist.pop(0)  # Remove first element of the list
    alist.append(int(end))  # Add the value "end" to the list


    ''' Initialization Data '''
    # Reads ONLY 1 time step and returns the dictionary of possible combinations
    t = time_steps[0]
    t_end = time_steps[-1]
    view.ViewTime = t
    Render()
    sample_bound_x, sample_bound_y, sample_bound_z = [], [], []
    start_initial = time.time()
    print(formatting_codes.green_text("Starting Initialization"))
    print(formatting_codes.green_text("Bin Frequency for Post-processing Automatically set at %d" % bin_freq))

    ''' 
    MODEL PROPERTIES // BOUNDARIES
    '''

    ''' Combinations of the the various material boundaries in the simulation results '''
    # Also returns % composition within the model
    # Create a list of the material types in the model
    for i in range(0, output.GetNumberOfCells()):
        mineral_bound.append(int(cells.GetArray('material property ID').GetTuple1(i)))

    # Return a list of unique material types
    mineral_bound_2 = list(set(mineral_bound))
    const = Counter(mineral_bound)  # counts the frequency of each mineral type
    list_combo = [x for x in itertools.combinations_with_replacement(mineral_bound_2, 2)]  # Creates a list of all possible combinations

    # Combines the Groups, creating a dictionary
    # global material_combinations
    for key, val in enumerate(list_combo):  # Enumerates the listcombo making it "Group #")
        material_combinations["Group %d" % key] = val

    # Display the various elements in the models
    print(formatting_codes.green_text("Model Element Properties:"))
    for key, val in const.items():
        per = "{0:.2f}".format(float(val) / float(output.GetNumberOfCells()) * 100)
        print ("\tElement ID %s - Percentage Composition %s %%" % (formatting_codes.bold_text(key), formatting_codes.bold_text(per)))
    print ("\tPlease be aware of the presence of the Platens in the %'s")

    # If specimen center is not defined as input
    # Obtains and Calculates the model bounds and specimen center.
    # Assumes symmetry.

    model_bounds = output.GetBounds()  # xmin, xmax, ymin, ymax, zmin, zmax
    x_center = (model_bounds[0] + model_bounds[1]) / 2
    y_center = (model_bounds[2] + model_bounds[3]) / 2
    z_center = (model_bounds[4] + model_bounds[5]) / 2
    specimen_center = [x_center, y_center, z_center]

    # If Platen Property ID not defined.
    # Lookup cell element ID on the top center and then trace points
    # Using this information, we obtain the platen prop ID
    if platen_cells_prop_id == -1:
        top_center_point = [(model_bounds[0] + model_bounds[1]) / 2, model_bounds[3], model_bounds[5]]
        top_center_cell = output.FindCell(top_center_point, None, 0, 1e-4, subId, pcoords, w)
        if top_center_cell == -1:
            print(formatting_codes.red_text("Unable to identify Platen ID Correctly. Proceed with caution or use -e flag."))
        platen_cells_prop_id = int(cells.GetArray('material property ID').GetTuple1(top_center_cell))
        print("\tPlaten Material ID found as %s" % formatting_codes.bold_text(platen_cells_prop_id))
    else:
        print("\tPlaten Material ID user defined as %s" % formatting_codes.bold_text(platen_cells_prop_id))

    if platen_points_prop_id is None:
        platen_status = "undefined"
        platen_points_prop_id = []
    # else:
    #     platen_status = "defined"


    ''' HAVE TO LOOP ONCE TO GET EXTENTS THEN LOOP AGAIN TO GET POINTS ON BOUNDARY '''
    for i in xrange(output.GetNumberOfCells()):
        # Get the Cell Object & Point Data
        cell = output.GetCell(i)
        # point_data = output.GetPointData()
        # Get the property_id of the current Cell
        cell_pr = int(cells.GetArray(cell_pr_txt).GetTuple1(i))
        # Get the id list of cell having the prescribed element property_id
        # idlist = vtk.vtkIdList()
        id_list = output.GetCellPoints(i, idlist)
        cellidlist = ([idlist.GetId(k) for k in range(idlist.GetNumberOfIds())])  # Nodes of the CellID
        # Check if the property_id is equal to the prescribed element property_id
        if (cell_pr == platen_cells_prop_id):
            if platen_status == "undefined":
                # Lookup each point of the cell and return the boundary condition ID.
                for j in range(0, len(cellidlist)):
                    # print point_data.GetArray('boundary condition ID')
                    b_cond = points.GetArray('boundary condition ID').GetTuple1(cellidlist[j])
                    platen_points_prop_id.append(int(b_cond))
        else:
            # Make a list of all X, Y, and Z of the points not in the element property ID
            for j in range(0, len(cellidlist)):
                all_xs, all_ys, all_zs = output.GetPoint(cellidlist[j])
            sample_bound_x.append(all_xs)
            sample_bound_y.append(all_ys)
            sample_bound_z.append(all_zs)

    # Return unique elements (i.e., unique boundary id of platens)
    platen_points_prop_id = list(set(platen_points_prop_id))


    ''' 
    OBTAINING VITAL SAMPLE INFORMATION FOR FURTHER PROCESSING
    '''

    global sample_x_min, sample_x_max, sample_y_min, sample_y_max, width, height, thickness
    sample_x_min, sample_x_max  = min(sample_bound_x), max(sample_bound_x)
    sample_y_min, sample_y_max  = min(sample_bound_y), max(sample_bound_y)
    width = sample_x_max - sample_x_min
    height = sample_y_max - sample_y_min
    thickness = max(sample_bound_z) - min(sample_bound_z)

    # Lookup cell element ID on the top left of the SAMPLE to define test type
    # If no material is found, cell ID == -1 is returned
    # Meaning that probably it is a BD Sample.
    check_edge_point = [sample_x_max, sample_y_max, 0]
    check_edge_cell = output.FindCell(check_edge_point, None, 0, 1e-4, subId, pcoords, w)
    global sim_type
    if check_edge_cell == -1:
        sim_type = "BD Simulation"
    else:
        sim_type = "UCS Simulation"


    # If Strain gauges on, identify location using specimen center
    # List of points for the vertical (v) and horizontal (h) strain gauges
    global pv, ph
    pv, ph = [],[]

    if gauge_length != gauge_width or gauge_length != 0:
        print(formatting_codes.green_text("Strain gauge (SG) calculation Enabled"))
        do_strain_gauge_analysis = True
    else:
        do_strain_gauge_analysis = False

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
        print('\tDimensions of SG are %s x %s' % (gauge_length, gauge_width))
        cv, ch = [], []
        for ps in xrange(0, len(pv)):
            cv.append(int(output.FindCell(pv[ps], None, 0, 1e-3, subId, pcoords, w)))
            ch.append(int(output.FindCell(ph[ps], None, 0, 1e-3, subId, pcoords, w)))
            cellArray.InsertNextCell(output.GetCell(cv[ps]))
            cellArray.InsertNextCell(output.GetCell(ch[ps]))
            # SetCells of the output with only the selected analyzed cells
        # output.SetCells(vtk.VTK_TRIANGLE, cellArray)
        # print cv, ch
        print('\tVertical Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (pv, cv))
        print('\tHorizontal Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (ph, ch))

    '''
    PROCESS THE INITIAL MODEL MESH
        - Platens excluded
        - Calculate Size / Area / Length
        - Draw Rosette Diagram
    '''

    print(formatting_codes.green_text("\nMesh Element Properties:"))
    print(formatting_codes.red_text("\tExcluding Platens"))
    # dfn_output_types, initial_dfn_line_list = [], []
    create_post_processing(dir_path, "post_processing") # Create post-processing folder
    areas, overall_lines, rose_angle, damage = [], [], [], []
    for i in range(0, output.GetNumberOfCells()):
        # print(cells.GetArray('material property ID').GetTuple1(i), platen_cells_prop_id)
        if cells.GetArray('material property ID').GetTuple1(i) != float(platen_cells_prop_id):
            # dfn_output_types.append(cells.GetArray('material property ID').GetTuple1(i))
            # Triangle points rotate anti-clockwise
            # Get Points
            point_a = output.GetPoint((i * 3) + 0)  # Get Point A of the Triangle
            point_b = output.GetPoint((i * 3) + 1)  # Get Point B of the Triangle
            point_c = output.GetPoint((i * 3) + 2)  # Get Point C of the Triangle

            # Use points to construct Triangles
            overall_lines.append([point_a, point_b])  # Construct line #1 of Triangle
            overall_lines.append([point_b, point_c])  # Construct line #2 of Triangle
            overall_lines.append([point_c, point_a])  # Construct line #3 of Triangle

            # Use the points to calculate the area of the triangle
            tri_area = area(point_a, point_b, point_c)
            areas.append(tri_area)

    # Cleanup for lines that coincide (triangle edges)
    # remove duplicate items from nested list
    fset = set(frozenset(x) for x in overall_lines)
    global new_overall_lines
    # return list of unique items
    new_overall_lines = [list(x) for x in fset]

    # Calculate edge length and angle for each mesh edge (line)
    for i in new_overall_lines:
        ax, ay, az = list(i[0])
        bx, by, bz = list(i[1])

        # Calculate the slope of the line in degrees
        angle_deg = (math.degrees(math.atan2(by - ay, bx - ax)))

        # Calculate the length of the line
        cal_len = math.sqrt((by - ay) ** 2 + (bx - ax) ** 2)

        rose_angle.append(angle_deg)
        damage.append(cal_len)
        # print ax, ay, az, bx, by, bz, angle_deg, cal_len  # Uncomment to see line X Y Z Angle Length

    # Display Mesh Information
    print("\tMesh Element Length:\tMax %.2f\tMin %.2f\tAverage %.2f" % (max(damage), min(damage), sum(damage) / len(damage)))
    print("\tMesh Element Area:\t\tMax %.2f\tMin %.2f\tAverage %.2f" % (max(areas), min(areas), sum(areas) / len(areas)))
    print(formatting_codes.red_text("\tCreating Rosette"))
    # Draw Rosette Diagram
    # rose_illustration(post_processing, width, height, len(new_overall_lines), reset_scale, rose_angle, damage, "2D_mesh_rosette")
    rose_illustration(rose_angle, damage, "2D_mesh_rosette")

    global mesh_avg_length, mesh_avg_area
    mesh_avg_area = sum(areas) / len(areas)
    mesh_avg_length = sum(damage) / len(damage)



    ''' 
    GET MORE MODEL INFORMATION 
        - TOP//BOTTOM PLATEN
        - POINTS ON MODEL EXTERIOR
    '''

    global circum_coord
    circum_coord = []
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
                # Read the co-ordinates of all the cells
                xs, ys, zs = circum_point = output.GetPoint(cellidlist[j])
                if sim_type == "UCS Simulation":
                    # Identifies if cell co-ordinate is on the edge.
                    if xs in [sample_x_min, sample_x_max] or ys in [sample_y_min, sample_y_max]:
                        modelextend.append(cellidlist[j])  # "Point ID on boundary"
                        cellextent.append(i)  # "cells on boundary"
                        circum_coord.append(circum_point)  # X/Y/Z of coordinates on edge
                else:
                    # Calculate Radial Distance
                    # (x - center_x) ^ 2 + (y - center_y) ^ 2 < radius ^ 2
                    r = ((xs - specimen_center[0])**2) + ((ys - specimen_center[1])**2)
                    # Identifies if co-ordinate is on the edge.
                    if ((sample_x_min - specimen_center[0])**2 - mesh_avg_length) <= r <=  ((sample_x_min - specimen_center[0])**2 + mesh_avg_length):
                        modelextend.append(cellidlist[j])  # "Point ID on boundary"
                        cellextent.append(i)  # "cells on boundary"
                        circum_coord.append(list(circum_point))  # X/Y/Z of coordinates on edge

    # Unique list of cells, as cell may repeat if two points are on the edge.
    cellextent = list(set(cellextent))
    circum_coord = [list(y) for y in set([tuple(x) for x in circum_coord])]

    # Calculate initial area using points on boundary and center of model

    global init_area
    init_area = area_calculation(output, modelextend, specimen_center)

    '''
    IDENTIFY IF THE MODEL HAS PRE-EXISTING FRACTURES 
        - Build initial DFN list based on number of points
        - Plot the Rosette diagram for the existing DFN
        - Platens excluded
        - Model Extents excluded
    '''

    dfn_output_types, initial_dfn_line_list = [], []

    if list_of_files_dfn is None or not list_of_files_dfn:
        print(formatting_codes.red_text("\nNo initial DFN found"))
    else:
        # Read and obtain the different types of ID of the DFN
        for i in range(0, output_dfn.GetNumberOfCells(), 2):
            dfn_output_types.append(dfn_joint_info.GetArray('property ID').GetTuple1(i))
        dfn_output_types = list(set(dfn_output_types))

        # If there is only one type of DFN, use that information to build the DFN Rosette
        if len(dfn_output_types) == 1:
            for i in range(0, output_dfn.GetNumberOfPoints()):
                # Build initial dfn list based on number of points.
                initial_dfn_line_list.append(output_dfn.GetPoint(i))
            initial_dfn_line_list = list(set(initial_dfn_line_list))
            print("Creating DFN Rosette")
            illustration(output, cells, output_broken, broken_joint_info, output_dfn, dfn_joint_info, platen_cells_prop_id, dfn_output_types, "2D_initial_fracture_intensity_rosette.pdf")
        else:
            # Since broken joints is only Points.
            # Use point to find cell and check if cell is on the model extent.
            for i in range(0, output_broken.GetNumberOfCells(), 2):
                if broken_joint_info.GetArray('failure mode').GetTuple1(i) == 0:
                    ax, ay, az = point_a = output_broken.GetPoint(i * node_skip)  # Get Point A of the line
                    bx, by, bz = point_b = output_broken.GetPoint(
                        (i * node_skip) + node_skip)  # Get Point B of the line
                    cell_1 = output.FindCell(point_a, None, 0, 1e-4, subId, pcoords,
                                             w)  # Lookup the cell of the point
                    cell_2 = output.FindCell(point_b, None, 0, 1e-4, subId, pcoords,
                                             w)  # Lookup the cell of the point
                    if cell_1 not in cellextent and cell_2 not in cellextent:
                        material_id_cell_1 = cells.GetArray('material property ID').GetTuple1(
                            cell_1)  # Lookup the material ID of the cell
                        material_id_cell_2 = cells.GetArray('material property ID').GetTuple1(
                            cell_2)  # Lookup the material ID of the cell
                        # Then confirm that cell is not in platen.
                        if material_id_cell_1 != platen_cells_prop_id or material_id_cell_2 != platen_cells_prop_id:
                            initial_dfn_line_list.append(
                                output_dfn.GetPoint(i))  # Build initial dfn list based on number of points.
                initial_dfn_line_list = list(set(initial_dfn_line_list))
        illustration(output, cells, output_broken, broken_joint_info, output_dfn, dfn_joint_info, platen_cells_prop_id, dfn_output_types, "2D_initial_fracture_intensity_rosette.pdf")


    '''
    OUTPUT Information about model to Terminal
    '''
    # Display some important information on Terminal
    print(formatting_codes.green_text("\nModel Properties:"))
    print (formatting_codes.green_text("\t" + sim_type))
    print ("\tNo. of Nodes in Model: %s \n\tNo. of Mesh Elements: %s" % (formatting_codes.bold_text(output.GetNumberOfPoints()), formatting_codes.bold_text(output.GetNumberOfCells())))
    print ("\tModel Bounds: \txmin = %.2f,\txmax = %.2f\n\t\t\t\t\tymin = %.2f,\tymax = %.2f\n\t\t\t\t\tzmin = %.2f,\tzmax = %.2f" % tuple(
        model_bounds))
    print ("\tSpecimen Dimensions: \tX = %.2f, Y = %.2f, Z = %.2f" % (width, height, thickness))
    print ("\tSpecimen Center: \tX = %.2f, Y = %.2f, Z = %.2f" % tuple(specimen_center))
    if not platen_points_prop_id:
        print(formatting_codes.red_text("Unknown platen boundary conditions. Should be one of these values %s" % const.keys()))
    elif len(platen_points_prop_id) == 1 or -1 in platen_points_prop_id:
        print(formatting_codes.red_text("Platen Boundary Conditions: %s" % formatting_codes.bold_text(platen_points_prop_id)[0:]))
        print(formatting_codes.red_text("Please recheck platen element property."))
    else:
        print ("\tPlaten Boundary Conditions: %s" % formatting_codes.bold_text(platen_points_prop_id)[0:])
    print ("\tNo. of Nodes on Specimen Boundary: %s " % formatting_codes.bold_text(len(modelextend)))
    print ("\tInitial Area: %s" % formatting_codes.bold_text(round(init_area,2)))
    print ("\nInitialization Complete: %s" % formatting_codes.bold_text(formatting_codes.calc_timer_values(time.time() - start_initial)))

    # exit("Initialization Complete")

    '''
    QUICK Crack Analysis
        - Goes through the time steps (specific to the frequency) to identify the maximum number of cracks in a specific direction
        - Uses this information while creating the rosette diagram as the maximum scale
        - This avoids the fact that the scale is dynamic, its static to the maximum mode of cracks in the entire simulation
    '''

    global reset_scale
    reset_scale = 0
    for t in time_steps:
        # print t
        if t in alist:
            view.ViewTime = t
            view.StillRender()
            UpdatePipeline(t)
            global max_scale_angle
            max_scale_angle = []
            for i in range(0, output_broken.GetNumberOfCells(), 2):
                if broken_joint_info.GetArray('failure mode').GetTuple1(i) not in [0, 4]:
                # Co-relating Broken Joints TO Basic Cell #
                # // THIS IS 2D CALCULATIONS! // #
                    ax, ay, az = point_a = output_broken.GetPoint(i * node_skip)
                    ax1, ay1, az1 = point_a1 = output_broken.GetPoint((i * node_skip) + 1)
                    bx, by, bz = point_b = output_broken.GetPoint((i * node_skip) + node_skip)
                    bx1, by1, bz1 = point_b1 = output_broken.GetPoint((i * node_skip) + node_skip + 1)
                    cen_ax, cen_ay, cen_az = (ax + bx1) / 2, (ay + by1) / 2, (az + bz1) / 2
                    cen_bx, cen_by, cen_bz = (bx + ax1) / 2, (by + ay1) / 2, (bz + az1) / 2

                    ## ROSE DIAGRAM CALCULATIONS
                    #  Calculate the length of the line
                    cal_len = math.sqrt((cen_by - cen_ay) ** 2 + (cen_bx - cen_ax) ** 2)

                    #  Calculate the slope of the line in degrees
                    angle_deg = (math.degrees(math.atan2(cen_by - cen_ay, cen_bx - cen_ax)))
                    # print (angle_deg)
                    if angle_deg <= 0:
                        angle_deg = (angle_deg + 180)  # All to Positive

                    max_scale_angle.append(angle_deg)
            counts, bins = numpy.histogram(max_scale_angle, bins=18, range=(0, 180))

            # Updates the scale for the Rosette of the crack elements
            if reset_scale < max(counts):
                reset_scale = max(counts)

    # exit("Reset Scale of Rose Diagram")

    UpdatePipeline(time_steps[0])

    processOutput(dir_path, list_of_files, width, height, specimen_center, platen_cells_prop_id, platen_points_prop_id, gauge_length, gauge_width, do_strain_gauge_analysis, modelextend, output_file_name)


'''
Process ParaView output of a 2D UCS simulation. Writes the data to a CSV file.
    Inputs:
        - Directory where output files are located
        - File extension of output files (usually '.vtu')
        - Type of output files, for instance use '_basic_' for the basic/main output files
        - Element (cell) property_id of the platens
        - Point (node or boundary condition) property_id of platen nodes as a list: top platen, bottom platen
        - Diameter (width) of the sample in mm
        - Height of the sample in mm
        - Coordinates of the center of the specimen [x, y, z]
        - Strain gauge length (mm)
        - Strain gauge width (mm)
        - Flag to indicate if strain gauge analysis should be done, given relevant data is provided
          (the flag could become False if this data is not available)
        - Axis of loading: x=0; y=1; z=2 
        - Shape of the specimen cross-section (circle or square)
'''


def processOutput(dir_path, list_of_files, diameter, height, specimen_center, platen_cells_prop_id, platen_points_prop_id, gauge_length, gauge_width, do_strain_gauge_analysis, modelextend, output_file_name):
    # By default this is a 2D model and analysis
    print(formatting_codes.green_text("\nProcessing Started"))
    start_initial_process = time.time()
    broken_time = 0
    seis_time = 0

    num_dimensions = 2

    datafile = open(os.path.join(post_processing, output_file_name), 'w',newline='')
    row = csv.writer(datafile, delimiter=',')

    brokenjoint_file = 'brokenjoints.csv'
    broken_datafile = open(os.path.join(post_processing, brokenjoint_file), 'w',newline='')
    brokenjoint_row = csv.writer(broken_datafile, delimiter=',')

    if list_of_files_seismic:
        source_file = 'sourcejoints.csv'
        source_datafile = open(os.path.join(post_processing, source_file), 'w',newline='')
        source_row = csv.writer(source_datafile, delimiter=',')

    # Data titles for crack type and number of broken elements
    br_row_titles = ['Time Step', "Broken Joint ID", "Cell_1", "Cell_2",
                     "Validation_cell_seismic", "Material ID_1", "Material ID_2",
                     "Group Type", "Crack Type", "Failure Mode", "Broken Joint Area",
                     "event energy", "kinetic energy at failure", "kinetic energy at yielding", "Angle"]

    brokenjoint_row.writerow(br_row_titles)

    if list_of_files_seismic:
        # Data titles for crack type and number of seismic post-processing
        source_row_titles = ["Time Step", "Node_0_X", "Node_0_Y", "Node_0_Z",
                             "Node_1_X", "Node_1_Y", "Node_1_Z",
                             "Node_2_X", "Node_2_Y", "Node_2_Z",
                             "Node_3_X", "Node_3_Y", "Node_3_Z",
                             "Event Energy", "Event Time Step", "Failure Mode", "Failure Time Step",
                             "Kinetic Energy at Failure", "Kinetic Energy at Yielding",
                             "Yielding Time Step", "Center X", "Center Y", "Center Z"]

        source_row.writerow(source_row_titles)

    if sim_type == "UCS Simulation":
        # Data titles and data of mechanical properties
        row_titles = ["Time Step (-)", "Platen displacement y (mm)", "Strain y from platen (%)",
                      "Strain x (%)", "Strain y (%)", "Volume (mm3)", "Delta Volume (mm3)",
                      "Platen Force (kN)", "Axial stress (MPa)", "Event Count",
                      "Volumetric strain (%)",
                      "Material property ID of platens = " + str(platen_cells_prop_id),
                      "Boundary condition ID of platen = " + str(platen_points_prop_id),
                      "Specimen diameter = " + str(diameter),
                      "Specimen height = " + str(height),
                      "Specimen centre = " + str(specimen_center),
                      "Strain gauge length = " + str(gauge_length),
                      "Strain gauge width = " + str(gauge_width)]

        # Perform strain gauge analysis only if the relevant data was provided as inputs
        # if specimen_center is None or not specimen_center or not do_strain_gauge_analysis:
        if gauge_length == gauge_width and gauge_length == 0:
            do_strain_gauge_analysis = False
            # Remove row titles related to strain gauges
            row_titles.remove('Strain x (%)')
            row_titles.remove('Strain y (%)')
            row_titles.remove("Specimen centre = " + str(specimen_center))
            row_titles.remove("Strain gauge length = " + str(gauge_length))
            row_titles.remove("Strain gauge width = " + str(gauge_width))
    else:
        # Data titles and data of the BD simulation
        row_titles = ['Time Step (-)', 'Y displacement(mm)',
                      'Platen Force (kN)', 'Indirect tensile stress (MPa)',
                      "Material property ID of platens = " + str(platen_cells_prop_id),
                      "Boundary condition ID of platen = " + str(platen_points_prop_id),
                      "Disc diameter = " + str(diameter)]

    row.writerow(row_titles)

    # Show the grid
    Show(pv_grid_reader)

    # Define a new cell array for the output
    cellArray = vtk.vtkCellArray()
    cellTypes = []

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

    # Dummy variables for FindCell
    p = [0.0, 0.0, 0.0]
    subId = vtk.mutable(0)
    pcoords = [0.0, 0.0, 0.0]
    w = [0.0, 0.0, 0.0]

    # Get the available time steps
    time_steps = pv_grid_reader.TimestepValues
    # print 'time steps', time_steps

    # Check the current view time
    global view
    view = GetActiveView()
    # UpdatePipeline(time_steps[0])

    # Initial sample volume
    initial_volume = init_area
    area = init_area
    initial_perimeter = init_area

    # Show progress bar
    print_progress(0, len(list_of_files), prefix='Progress:', suffix='Complete')

    # Loop over time steps and perform calculations
    # This portion needs to go to pool...
    error_array, cluster, failure_modes = [], [], []  # For QAQC Error & Cluster CrackType
    mineral_bound, combos, breakdown, fail_count, crack_types, crack_count = [], [], [], [], [], []
    fail_modes, fail_angles = [], []  # reset every time step to match dfn lines against respective broken joints
    for t in time_steps:
        cluster, failure_modes, crack_types = [], [], []
        event_count = 0
        gc_counter = t

        view.ViewTime = t
        Render()

        ###----------------------------------
        ###      VOLUMETRIC STRAIN ANALYSIS
        ###----------------------------------

        # Loop over the nodes on the specimen boundary to obtain their current co-ordinate
        # Calculate area using points on boundary and center of model
        poly_area = area_calculation(output, modelextend, specimen_center)
        # Volumetric strain calculation for the entire sample
        d_volume = poly_area - init_area
        if d_volume is not None and initial_volume >= 0.0:
            volumetric_strain = ((d_volume) * 100.0 / initial_volume)

        # print t, volume, volumetric_strain

        ###----------------------------------
        ###      STRAIN GAUGE ANALYSIS
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
        ###        PLATEN ANALYSIS
        ###----------------------------------
        top_platen_cellArray = vtk.vtkCellArray()
        bottom_platen_cellArray = vtk.vtkCellArray()

        # Initialize the array of platen forces and displacements with zero (x,y,z)
        top_platen_force = [0.0, 0.0, 0.0]
        bottom_platen_force = [0.0, 0.0, 0.0]
        avg_top_platen_disp = [0.0, 0.0, 0.0]
        avg_bottom_platen_disp = [0.0, 0.0, 0.0]

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
        # The skip is made in increments of 4 (2D crack rectangle)
        # Match is made on the coordinates of the brokenjoint element.
        idlist = vtk.vtkIdList()  # Dummy for idList
        dfn_line_list = [] # reset every time step to match dfn lines against respective broken joints

        # Build DFN list based current position.
        if list_of_files_dfn:
            for f in range(0, output_dfn.GetNumberOfPoints()):
                dfn_line_list.append(output_dfn.GetPoint(f))
        dfn_line_list = [list(y) for y in set([tuple(x) for x in dfn_line_list])]  # reduces list to unique sets

        crack_cluster = time.time()
        for i in range(0, output_broken.GetNumberOfCells(), 2):

            # if 0 < broken_joint_info.GetArray('failure mode').GetTuple1(i) < 4:
            if broken_joint_info.GetArray('failure mode').GetTuple1(i) not in [0, 4]:
                # Co-relating Broken Joints TO Basic Cell #
                # // THIS IS 2D CALCULATIONS! // #
                ax, ay, az = point_a = output_broken.GetPoint(i * node_skip)
                ax1, ay1, az1 = point_a1 = output_broken.GetPoint((i * node_skip) + 1)
                bx, by, bz = point_b = output_broken.GetPoint((i * node_skip) + node_skip)
                bx1, by1, bz1 = point_b1 = output_broken.GetPoint((i * node_skip) + node_skip + 1)
                cen_ax, cen_ay, cen_az = (ax + bx1) / 2, (ay + by1) / 2, (az + bz1) / 2
                cen_bx, cen_by, cen_bz = (bx + ax1) / 2, (by + ay1) / 2, (bz + az1) / 2

                ## ROSE DIAGRAM CALCULATIONS
                #  Calculate the length of the line
                cal_len = math.sqrt((cen_by - cen_ay) ** 2 + (cen_bx - cen_ax) ** 2)

                #  Calculate the slope of the line in degrees
                angle_deg = (math.degrees(math.atan2(cen_by - cen_ay, cen_bx - cen_ax)))
                if angle_deg < 0:
                    angle_deg = (angle_deg + 180) # All to Positive

                ## FIND CELL ID USING POINTS (from basic)
                cell_1 = output.FindCell(point_a, None, 0, 1e-4, subId, pcoords, w)
                cell_2 = output.FindCell(point_b, None, 0, 1e-4, subId, pcoords, w)

                if t in alist[:]:
                    fail_angles.append(angle_deg)

                #  SEISMIC DATA MANIPULATION #
                try:
                    if list_of_files_seismic:
                        # print seismic_point_info
                        seismic_point = seismic_point_info.GetArray('node 0').GetTuple3(i // 2)
                        seismic_point_coordinates = [seismic_point[0], seismic_point[1], seismic_point[2]]
                        cell_seismic = int(output.FindCell(seismic_point_coordinates, None, 0, 1e-4, subId, pcoords, w))
                except AttributeError:
                    print("ERROR IN INPUT FILE")
                    cell_seismic = cell_1

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
                    if broken_joint_info.GetArray('failure mode').GetTuple1(i) < 5:
                        # print(broken_joint_info.GetArray('failure mode').GetTuple1(i))
                        fail_mode = broken_joint_info.GetArray('failure mode').GetTuple1(i)  # Get Array "FAILURE MODE" of SPECIFIC CELL
                    else:
                        # print("GREATER THAN 5",broken_joint_info.GetArray('failure mode').GetTuple1(i))
                        fail_mode = broken_joint_info.GetArray('failure mode').GetTuple1(i) - 5  # Get Array "FAILURE MODE" of SPECIFIC CELL
                    if t in alist[:]:
                        fail_modes.append(fail_mode)
                except AttributeError:
                    print("Minor Error")
                    material_id_cell_1 = 0
                    material_id_cell_2 = 1
                    # fail_mode = 0

                # Clustering Failure Mode (material ids) - sorted in ascending order to avoid repetition
                if t in alist[:]:
                    cluster.append(tuple(sorted([int(material_id_cell_1), int(material_id_cell_2)])))
                for key, values in material_combinations.items():
                    if values == tuple(sorted([int(material_id_cell_1), int(material_id_cell_2)])):
                        group_type = key

                # Clustering Failure Modes (Tensile Dominant / Shear Dominant / Mix Mode)
                if t in alist[:]:
                    if 1 <= fail_mode <= 1.5:
                        failure_modes.append("Tensile Dominant")
                    elif 1.5 < fail_mode <= 2.0:
                        failure_modes.append("Shear Dominant")
                    elif fail_mode == 3.0:
                        failure_modes.append("Mixed Mode")

                # Identify the type of Crack
                # A) NOT Same material => IntERgranular
                # B) Same material AND on DFN => IntERgranular
                # C) Same material NOT on DFN => IntRAgranular

                if material_id_cell_1 != material_id_cell_2:
                    crack_type = "intergranular_Material"
                else:
                    # Since DFN is not defined.
                    # There is no possibility of knowing the location of the crack with respect grain boundary.
                    # A) NOT Same material => IntERgranular
                    if list_of_files_dfn is None or not list_of_files_dfn:
                        crack_type = "intergranular / intragranular"
                    else:
                        # Point List of Broken Joint
                        output_broken.GetCellPoints(i, idlist)
                        broken_joint_idList_1 = ([idlist.GetId(k) for k in range(idlist.GetNumberOfIds())])
                        output_broken.GetCellPoints(i + 1, idlist)
                        broken_joint_idList_2 = ([idlist.GetId(k) for k in range(idlist.GetNumberOfIds())])
                        # Looking for X/Y/Z of Broken Joint in DFN
                        for cor_idList_1, cor_idList_2 in zip(range(0, len(broken_joint_idList_1)), range(0, len(broken_joint_idList_2))):
                            xyz_idList_1 = output_broken.GetPoint(broken_joint_idList_1[cor_idList_1])
                            xyz_idList_2 = output_broken.GetPoint(broken_joint_idList_2[cor_idList_2])
                            if list(xyz_idList_1) in dfn_line_list and list(xyz_idList_2) in dfn_line_list:
                                # B) Same material AND on DFN => IntERgranular
                                crack_type = "intergranular_DFN"
                            else:
                                # C) Same material NOT on DFN => IntRAgranular
                                crack_type = "intragranular"

                if t in alist[:]:
                    crack_types.append(crack_type)

                ###----------------------------------
                ###        FILE OUTPUT
                ###----------------------------------

                ###-----------------------------
                ### BROKEN JOINT FILE OUTPUT ###
                ###-----------------------------

                if list_of_files_seismic:
                    if cell_seismic == cell_1 or cell_seismic == cell_2:
                        br_row = [t]
                        br_row.extend([i, int(cell_1), int(cell_2), int(cell_seismic), material_id_cell_1, material_id_cell_2, group_type, crack_type, fail_mode, broken_joint_info.GetArray('area').GetTuple1(i), seismic_point_info.GetArray('event energy').GetTuple1(i // 2), seismic_point_info.GetArray('kinetic energy at failure').GetTuple1(i // 2), seismic_point_info.GetArray('kinetic energy at yielding').GetTuple1(i // 2), angle_deg])
                        broken_time = time.time()
                        brokenjoint_row.writerow(br_row)
                        broken_time -= time.time()
                        # print(seis_time)
                else:
                    # All Seismic info removed as no Seismic files found.
                    if cell_2 != cell_1:
                        rowData = [t]
                        rowData.extend([i, int(cell_1), int(cell_2), "NA", material_id_cell_1, material_id_cell_2, group_type, crack_type, fail_mode, broken_joint_info.GetArray('area').GetTuple1(i), "NA", "NA", "NA", angle_deg])
                        broken_time = time.time()
                        brokenjoint_row.writerow(rowData)
                        broken_time -= time.time()
                        # print(broken_time)
                    else:
                        print("ERROR in %d: Check Cell contact at # %s # %s" % (t, cell_1, cell_2))

                ###------------------------
                ### SEISMIC FILE OUTPUT ###
                ###------------------------

                if list_of_files_seismic:
                    if cell_seismic == cell_1 or cell_seismic == cell_2:
                        x_tot, y_tot, z_tot = 0, 0, 0
                        event_count += 1
                        for n in range(0, node_skip):
                            x_tot += seismic_point_info.GetArray(n).GetTuple3(i // (node_skip // 2))[0]
                            y_tot += seismic_point_info.GetArray(n).GetTuple3(i // (node_skip // 2))[1]
                            z_tot += seismic_point_info.GetArray(n).GetTuple3(i // (node_skip // 2))[2]
                        x_avg = x_tot / node_skip
                        y_avg = y_tot / node_skip
                        z_avg = z_tot / node_skip
                        source_rowData = [t]
                        source_rowData.extend(
                            [seismic_point_info.GetArray('node 0').GetTuple3(i // 2)[0],
                             seismic_point_info.GetArray('node 0').GetTuple3(i // 2)[1],
                             seismic_point_info.GetArray('node 0').GetTuple3(i // 2)[2],
                             seismic_point_info.GetArray('node 1').GetTuple3(i // 2)[0],
                             seismic_point_info.GetArray('node 1').GetTuple3(i // 2)[1],
                             seismic_point_info.GetArray('node 1').GetTuple3(i // 2)[2],
                             seismic_point_info.GetArray('node 2').GetTuple3(i // 2)[0],
                             seismic_point_info.GetArray('node 2').GetTuple3(i // 2)[1],
                             seismic_point_info.GetArray('node 2').GetTuple3(i // 2)[2],
                             seismic_point_info.GetArray('node 3').GetTuple3(i // 2)[0],
                             seismic_point_info.GetArray('node 3').GetTuple3(i // 2)[1],
                             seismic_point_info.GetArray('node 3').GetTuple3(i // 2)[2],
                             seismic_point_info.GetArray('event energy').GetTuple1(i // 2),
                             seismic_point_info.GetArray('event time step').GetTuple1(i // 2),
                             seismic_point_info.GetArray('failure mode').GetTuple1(i // 2),
                             seismic_point_info.GetArray('failure time step').GetTuple1(i // 2),
                             seismic_point_info.GetArray('kinetic energy at failure').GetTuple1(i // 2),
                             seismic_point_info.GetArray('kinetic energy at yielding').GetTuple1(i // 2),
                             seismic_point_info.GetArray('yielding time step').GetTuple1(i // 2),
                             x_avg, y_avg, z_avg])
                        seis_time = time.time()
                        source_row.writerow(source_rowData)
                        seis_time -= time.time()
                    else:
                        print("ERROR in %d: Check Cell contact at # %s # %s # %d" % (t, cell_1, cell_2, cell_seismic))

        # print(time.time() - crack_cluster)
        ###-------------------------
        ### STRESS/STRAIN OUTPUT ###
        ###-------------------------
        # Convert forces from microN to kN and get the average forces & displacements
        avg_platen_force = [0.0, 0.0, 0.0]
        avg_platen_disp = [0.0, 0.0, 0.0]
        axis_of_loading = 1
        # print("No. of events ", t, event_count, crack_count)
        rowData = [t]
        for i in range(0, 3):
            # microN to kN
            avg_platen_force[i] = 0.5 * (abs(top_platen_force[i]) + abs(bottom_platen_force[i])) / 1.0e9
            avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

        if sim_type == "UCS Simulation":
            # stress in MPa (force in kN & area in mm^2)
            stress = avg_platen_force[axis_of_loading] / width * 1.0e3
            # strain calculated from platen displacements
            strain_from_platen = avg_platen_disp[axis_of_loading] / height * 100.0

            if do_strain_gauge_analysis:
                rowData.extend([avg_platen_disp[axis_of_loading], strain_from_platen, strain_x, strain_y,
                                poly_area, d_volume,
                                avg_platen_force[axis_of_loading], stress, event_count])
            else:
                rowData.extend([avg_platen_disp[axis_of_loading], strain_from_platen,
                                poly_area, d_volume,
                                avg_platen_force[axis_of_loading], stress, event_count])

            if d_volume is not None and initial_volume >= 0.0:
                rowData.extend([volumetric_strain, event_count])
            # print "\nRow data:", rowData

            # Append the data to the csv file
            his_time = time.time()
            row.writerow(rowData)
            his_time -= time.time()

        else:
            disc_area = math.pi * width
            # stress in MPa (force in kN & lengths in mm)
            stressY = 2.0 * avg_platen_force[axis_of_loading] / disc_area * 1.0e3
            # strain calculated from platen displacements
            rowData.extend([avg_platen_disp[axis_of_loading], avg_platen_force[axis_of_loading], stressY])

            # print rowData
            # Append the data to the csv file
            row.writerow(rowData)

        ###----------------------------------
        ###        GRAPHICAL OUTPUT
        ###----------------------------------

        # Create bins for graphical output
        if t in alist[:]:
            breakdown.append(Counter(cluster)) # No. of cracks based on Phase Interface
            fail_count.append(Counter(failure_modes)) # Failure Modes
            crack_count.append(Counter(crack_types)) # No. of Cracks

            # cluster, failure_modes, crack_types = [], [], []
            if crack_types:
                rose_illustration(fail_angles, fail_modes, t)
                fail_angles, fail_modes = [], []

        # Update progress bar
        # print("Hist\tBroken\tSiesmic\n")
        # print((his_time))
        # print((broken_time))
        # print(((seis_time)))
        # print(time.time()-start_initial_process)
        print_progress(int(t) + 1, len(list_of_files), prefix='Progress:', suffix='Complete')

    datafile.close()
    broken_datafile.close()
    if list_of_files_seismic:
        source_datafile.close()

    print("\nFinished writing processed Stress/Strain data to: %s" % formatting_codes.bold_text(os.path.join(post_processing, 'history.csv')))

    # /// RESULTS VERIFICATION /// #
    # /// CLUSTER THE RESULTS /// #
    # First checks to see if there are any broken joints
    # Moves on to clustering
    # return results to counter (creates bins for histogram plotting)
    # Inserted here to avoid any early errors.

    if not crack_types:
        print("\n\033[1;31;0mNO BROKEN JOINTS TO PROCESS\nTERMINATING\033[0m ")
    else:
        if list_of_files_seismic:
            print("Processing Seismic Clustering")

            ''' Qi's SEISMIC OUTPUT FORMAT!!! '''
            seismicclustering(post_processing)


        # import matplotlib.pyplot as plt
        # print Counter(cluster) # Uncomment if you wish to see the histogram results
        if node_skip == 4:
            # graphical_illustrations.plot_bar_from_counter(post_processing, str(begin), str(end), str(bin_freq),
            #                                               "main_graph")

            seismic_illustrations.load_process_file(post_processing, 10, bin_freq)
            process_b_value.load_process_file(post_processing)
            # import seismic_illustrations

            # graphical_illustrations.plot_bar_from_counter(Counter(cluster), breakdown, fail_count, crack_count, "main_graph")

    # Draw stress/strain graphs
    graphical_illustrations.single_graph(post_processing, str(sim_type))
    # single_graph(post_processing)

    ###----------------------------------
    ### COMPLETED PROCESSING FOLDER
    ###----------------------------------

    global_removal(output, cells, points, data_object)
    print ("\nProcessing Complete: %s" % formatting_codes.bold_text(formatting_codes.calc_timer_values(time.time() - start_initial_process)))
    print (formatting_codes.red_text("-------------------------------"))


'''
PLOT CURVES.
    Inputs:
        - File containing stress/strain information
    Output:
        - Stress/Strain & Stress/Volumetric Stain Curve in various formats.
        - Stress/Strain Curve in various formats.
'''


# def single_graph(sub_dirs):
#     # import graphical_illustrations
#     graphical_illustrations.single_graph(post_processing, str(sim_type))
    # command = "python /hdd/home/aly/Desktop/Dropbox/Python_Codes/irazu_post_processing/graphical_illustrations.py -graphs %s \"%s\"" % (post_processing, str(sim_type))
    # # print(command)
    # os.system(command)

'''
GRAPHICAL REPRESENTATION.
    Inputs:
        - Different broken Phase Interfaces
        - Crack Type (IntER/IntRA granular cracks)
        - Failure Modes (Tensile Dominant / Shear Dominant / Mixed Mode)
    Output:
        - ax1 = Histogram chart of broken Phase Interfaces.
        - ax2 = Stacked Bar Chart different broken Phase Interfaces based on Time Step.
        - ax3 = Stacked Bar Chart of Crack Type (IntER/IntRA granular cracks)
        - ax4 = Stacked Bar Chart for Fracture Mode Type (Tensile Dominant / Shear Dominant / Mixed Mode)
'''


def plot_bar_from_counter(counter, breaks, failures, types, graph_name, ax1=None):
    import json

    # break_series = {}
    # fail_type_series = {}
    # fail_series = {}

    # Temporary files
    # temp_counter = os.path.join(post_processing, 'temp_counter.csv')  # AX1
    # temp_failures = os.path.join(post_processing, 'temp_failures.csv')  # AX2
    # temp_types = os.path.join(post_processing, 'temp_types.csv')  # AX3
    # temp_break = os.path.join(post_processing, 'temp_break.csv')  # AX4

    #Write data from files
    # AX1 Data
    # with open(temp_counter, 'w') as writeFile:
    #     writer = csv.writer(writeFile)
    #     writer.writerows(counter.items())
    # writeFile.close()

    # AX2 Data
    # for key in {key for keys in failures for key in keys}:
    #     fail_series[key] = [(0 if key not in item else item[key]) for item in failures]
    #
    # # with open(temp_failures, 'w') as f:
    # #     json.dump(fail_series, f)
    #
    # # AX3 Data
    # for key in {key for keys in types for key in keys}:
    #     fail_type_series[key] = [(0 if key not in item else item[key]) for item in types]

    # with open(temp_types, 'w') as f:
    #     json.dump(fail_type_series, f)

    # # AX4 Data
    # for i in breaks:
    #     for key in {key for keys in breaks for key in keys}:
    #         break_series[key] = [(0 if key not in item else item[key]) for item in breaks]
    #     # Convert the Keys explicitly to Text to enable json handling
    #     break_series = {str(old_key): val for old_key, val in break_series.items()}

    # with open(temp_break, 'w') as f:
    #     json.dump(break_series, f)

    # import graphical_illustrations
    graphical_illustrations.plot_bar_from_counter(post_processing, str(begin), str(end), str(bin_freq), str(graph_name))

    # command = "python /hdd/home/aly/Desktop/Dropbox/Python_Codes/irazu_post_processing/graphical_illustrations.py -histogram %s %s %s %s %s %s" % (post_processing, str(begin), str(end), str(bin_freq), str(graph_name), str("ax1=None"))
    # os.system(command)


'''
SEISMIC CLUSTERING.
    Inputs:
        - Identifies and sorts similar events by the values of the event energy generated.
        - Output directory.
    Output:
        - The time step start and end of all the events identified during that time window.
        - CSV compliant to the MATLAB script for b/D Value processing.
'''


def seismicclustering(post_processing):
    global rose_angle, rose_mode, initial_azimuth, initial_dip_angle

    ###-------------------------------
    ### SEISMIC CLUSTERING OUTPUT ###
    ###-------------------------------

    # Read the broken joints output and cluster event based on #11 (event energy)
    reader = csv.reader(open(os.path.join(post_processing, "brokenjoints.csv")), delimiter=",")
    for TimeStep, BrokenJointID, Cell1, Cell2, Validation_cell_seismic, MaterialID_1, MaterialID_2, GroupType, CrackType, FailureMode, BrokenJoint, eventenergy, kineticenergyatfailure, kineticenergyatyielding, Angle in reader:
        seismicsortedlist = sorted(reader, key=operator.itemgetter(11), reverse=True)

    s_clust, break_angle = [], []
    rose_mode, rose_angle = [], []
    azim, initial_dip_angle, initial_azimuth = [], [], []
    mx, mn, b_angle, azim_azim = 0 , 0, 0, 0
    if not seismicsortedlist:
        print(formatting_codes.red_text("\nNO BROKEN JOINTS TO PROCESS\n"))
    else:
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

    ###-------------------------------
    ### source_mat (MATLAB) OUTPUT ###
    ###-------------------------------


    source_reader = csv.reader(open(os.path.join(post_processing, "sourcejoints.csv")), delimiter=",")
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

'''
CLEAR MEMORY
    - Clears pvpython Pipeline prior to loading next VTU outputs
    - Acts as garbage collector for python prior to loading next folder
'''


def global_removal(output, cells, points, data_object):
    del output, cells, points, data_object
    gc.collect()
    Delete(pv_grid_reader)
    Delete(pv_grid_reader_broken)
    if list_of_files_seismic:
        Delete(pv_grid_reader_seismic)
    if list_of_files_dfn:
        Delete(pv_grid_reader_dfn)
    print(formatting_codes.bold_text("Memory Cleared"))


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


'''
Calculate Area of Polygon
    Inputs:
        - Point ID of the INITIAL Model extents
        - Specimen Center
    Returns:
        - Updated co-ordinates of the model extents
        - Create a polygon using the co-ordinates
        - Use the Polygon to calculate the Polygon Area (Poly Area)
'''


def area_calculation(output, modelextend, specimen_center):
    # Get the X/Y/Z of the edge co-ordinates
    coordinates = []

    # global poly_area
    for j in modelextend:
        xs, ys, zs = output.GetPoint(j)
        coordinates.append([xs, ys, zs])

    # Sort points in clockwise direction using Origin as sample center
    # Sort on Angle tan(dy/dx); then on radial distance sqrt(dx^2+dy^2)
    coordinates.sort(key=lambda c: (math.atan2(c[1] - specimen_center[1], c[0] - specimen_center[0]),
                                   math.sqrt(
                                       (c[1] - specimen_center[1]) ** 2 + (c[0] - specimen_center[0]) ** 2)) )
    # Create Polygon THEN Calculate Area
    # temp_coor = os.path.join(post_processing, 'temp_coor.csv')
    # with open(temp_coor, 'w') as f:
    #     json.dump(coordinates, f)

    poly_area = geometry_calculator.calculate_polyarea(post_processing, coordinates)

    return poly_area



'''
Function to enable data collection for Rosette (INITIAL)
    Inputs:
        - Depending on condition 
            - CONDITION A) 2D with only one type of DFN property (DATA FROM DFN FILE)
            - CONDITION B) 2D with multi type of DFN property (DATA FROM BROKEN JOINT FILE AND CHECKED FROM BASIC FILE) 
    Returns:
       
'''


def illustration(output, cells, output_broken, broken_joint_info, output_dfn, dfn_joint_info, platen_cells_prop_id, dfn_output_types, out_name):
    rose_angle, damage = [], []
    dx, dy = 0, 0
    # CONDITION A
    # if node_skip == 4 and len(dfn_output_types) == 1:
    #     for i in range(0, output.GetNumberOfCells(), 2):
    #         # Co-relating Broken Joints TO Basic Cell #
    #         # // THIS IS 2D CALCULATIONS! // #
    #         ax, ay, az = point_a = output.GetPoint(i * node_skip)
    #         ax1, ay1, az1 = point_a1 = output.GetPoint((i * node_skip) + 1)
    #         bx, by, bz = point_b = output.GetPoint((i * node_skip) + node_skip)
    #         bx1, by1, bz1 = point_b1 = output.GetPoint((i * node_skip) + node_skip + 1)
    #         cen_ax, cen_ay, cen_az = (ax + bx1) / 2, (ay + by1) / 2, (az + bz1) / 2
    #         cen_bx, cen_by, cen_bz = (bx + ax1) / 2, (by + ay1) / 2, (bz + az1) / 2
    #
    #         # Check to see if BOTH points lie on the boundary
    #         if list(point_a) not in circum_coord and list(point_b) not in circum_coord:
    #             # Lookup the cell of the point
    #             cell_1 = output.FindCell(point_a, None, 0, 1e-4, subId, pcoords, w)
    #             cell_2 = output.FindCell(point_b, None, 0, 1e-4, subId, pcoords, w)
    #             # Lookup the material ID of the cell
    #             material_id_cell_1 = cells.GetArray('material property ID').GetTuple1(cell_1)
    #             material_id_cell_2 = cells.GetArray('material property ID').GetTuple1(cell_2)
    #             # Confirm that the cells do not lie in the platen
    #             if material_id_cell_1 != platen_cells_prop_id or material_id_cell_2 != platen_cells_prop_id:
    #
    #                 ## ROSE DIAGRAM CALCULATIONS
    #                 #  Calculate the slope of the line in degrees
    #                 angle_deg = (math.degrees(math.atan2(cen_by - cen_ay, cen_bx - cen_ax)))
    #                 # print (angle_deg)
    #                 if angle_deg < 0:
    #                     angle_deg = (angle_deg + 180)  # All to Positive
    #
    #                 #  Calculate the length of the line
    #                 cal_len = math.sqrt((cen_by - cen_ay) ** 2 + (cen_bx - cen_ax) ** 2)
    #
    #                 # if cal_len > 1.0:
    #                 # print(i, i * node_skip, (i * node_skip) + node_skip, cal_len, point_a, point_b, cell_1, cell_2)
    #
    #                 dx += math.fabs(cen_bx - cen_ax)
    #                 dy += math.fabs(cen_by - cen_ay)
    #                 rose_angle.append(angle_deg)
    #                 damage.append(cal_len)
    #
    #     # Print Mesh information
    #     print(formatting_codes.green_text("\nDFN Found - Processing Initial DFN"))
    #     print("\tCondition A")
    #     print("\t2D Single DFN - %s" % output.GetNumberOfCells())
    #     print("\tLength of DFN: Total %.2f\tMax %.2f \tMin %.2f\tAverage %.2f" % (
    #     sum(damage), max(damage), min(damage), sum(damage) / len(damage)))

    # CONDITION B
    if node_skip == 4 and len(dfn_output_types) > 0:
        for i in range(0, output_broken.GetNumberOfCells() - 1, 2):
            # Co-relating Broken Joints TO Basic Cell #
            # // THIS IS 2D CALCULATIONS! // #
            ax, ay, az = point_a = output_broken.GetPoint(i * node_skip)
            ax1, ay1, az1 = point_a1 = output_broken.GetPoint((i * node_skip) + 1)
            bx, by, bz = point_b = output_broken.GetPoint((i * node_skip) + node_skip)
            bx1, by1, bz1 = point_b1 = output_broken.GetPoint((i * node_skip) + node_skip + 1)
            cen_ax, cen_ay, cen_az = (ax + bx1) / 2, (ay + by1) / 2, (az + bz1) / 2
            cen_bx, cen_by, cen_bz = (bx + ax1) / 2, (by + ay1) / 2, (bz + az1) / 2
            # Check to see if BOTH points lie on the boundary
            if list(point_a) not in circum_coord and list(point_b) not in circum_coord:
                # Lookup the cell of the point
                cell_1 = output.FindCell(point_a, None, 0, 1e-4, subId, pcoords, w)
                cell_2 = output.FindCell(point_b, None, 0, 1e-4, subId, pcoords, w)
                # Lookup the material ID of the cell
                material_id_cell_1 = cells.GetArray('material property ID').GetTuple1(cell_1)
                material_id_cell_2 = cells.GetArray('material property ID').GetTuple1(cell_2)
                # Confirm that the cells do not lie in the platen
                if material_id_cell_1 != platen_cells_prop_id or material_id_cell_2 != platen_cells_prop_id:
                    ## ROSE DIAGRAM CALCULATIONS
                    #  Calculate the slope of the line in degrees
                    angle_deg = (math.degrees(math.atan2(cen_by - cen_ay, cen_bx - cen_ax)))
                    # print (angle_deg)
                    if angle_deg < 0:
                        angle_deg = (angle_deg + 180)  # All to Positive

                    #  Calculate the length of the line
                    cal_len = math.sqrt((cen_by - cen_ay) ** 2 + (cen_bx - cen_ax) ** 2)

                    # if cal_len > 1.0:
                    # print(i, i * node_skip, (i * node_skip) + node_skip, cal_len, point_a, point_b, cell_1, cell_2)

                    dx += math.fabs(cen_bx - cen_ax)
                    dy += math.fabs(cen_by - cen_ay)
                    rose_angle.append(angle_deg)
                    damage.append(cal_len)
            # else:
            #     print(i, "Point ON boundary")
            #     print(point_a)

        # Print DFN information
        if len(rose_angle) != 0:
            print(formatting_codes.green_text("\nDFN Found - Processing Initial DFN"))
            print("\tCondition B")
            print("\t2D Single Multi-DFN - Total No.: %s" % len(rose_angle))
            print("\tLength of DFN: Total %.2f\tMax %.2f \tMin %.2f\tAverage %.2f" % (sum(damage), max(damage), min(damage), sum(damage) / len(damage)))


    if len(rose_angle) != 0:
        count_initial_dfn = len(rose_angle)
        print(formatting_codes.red_text("\tCreating Rosette (Excluding Platens & Edges)"))
        rose_illustration(rose_angle, damage, "2D_initial_damage_intensity_rosette", dx, dy)



'''
Calculate normal unit vectors
    Inputs:
        - Vector
    Returns:
        - Unit vector of the vector
'''


def unit_vector(vector):
    return vector / numpy.linalg.norm(vector)


'''
Function to return angles between vectors
    Inputs:
        - Two Vectors
    Returns:
        - cos(Angle) between the vectors in degrees 
'''


def angle_between(v1, v2):
    # Calculate the unit vector of each vector
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    # Returns the cos angle in degrees
    return math.degrees(numpy.arccos(numpy.clip(numpy.dot(v1_u, v2_u), -1.0, 1.0)))


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


def rose_illustration(rose_angle, damage, fname, dx=0, dy=0, ax=None):

    # print(rose_angle, damage, fname,)
    if fname == alist[-1]:
        fname1 = "FinalFrame"
    else:
        fname1 = "NIL"

    # illus_var(damage)
    sub_post_processing = os.path.join(post_processing, "Simulation_Rose_Plot")
    Stitched = os.path.join(post_processing, "Stitched")

    if fname != "2D_initial_damage_intensity_rosette" and fname != "2D_mesh_rosette":
        screenshot(post_processing, fname)

    # print(rose_angle)
    # print(damage)
    # exit()
    # temp_rose_angle = os.path.join(post_processing, 'temp_rose_angle.csv')
    # with open(temp_rose_angle, 'w') as writeFile:
    #     for s in rose_angle:
    #         writeFile.write(str(s) + "\n")
    # writeFile.close()
    #
    # temp_str_damage = os.path.join(post_processing, 'temp_str_damage.csv')
    # with open(temp_str_damage, 'w') as writeFile:
    #     for s in damage:
    #         writeFile.write(str(s) + "\n")
    # writeFile.close()


    rose_Illustrations.rose_illustration(str(post_processing), str(width), str(height), str(len(new_overall_lines)), str(reset_scale), str(fname), str(fname1), rose_angle, damage, str("ax=None"), dx, dy )
    # exit()
    # command = (
    #             "python /hdd/home/aly/Desktop/Dropbox/Python_Codes/irazu_post_processing/rose_Illustrations.py %s %s %s %s %s %s %s %s" % (
    #     str(post_processing), str(width), str(height), str(len(new_overall_lines)), str(reset_scale), str(fname), str(fname1), str("ax=None")))

    # os.system(command)
    # os.remove(temp_rose_angle)
    # os.remove(temp_str_damage)


'''
Capture Paraview Screenshots at defined Frequency
# Parts of code obtained using Python Trace in Paraview v5.0.1 
    Inputs:
        - Location of Output Folder
        - Time Step Number (XXXX)
        - Paraview Parameters (Surface/Wireframe)
        - Using Cells and colored by Material Property ID and Failure Mode  
    Returns:
        - Set position of Camera
        - Save Screen shot
'''


def screenshot(fol, fil):
    global ffname

    # set image size
    view.ViewSize = [1000, 1500]  # [width, height]

    # set background color
    view.Background = [1, 1, 1]  # white

    # get display properties
    dp = GetDisplayProperties(pv_grid_reader)
    dp_broken = GetDisplayProperties(pv_grid_reader_broken)

    # change representation type
    dp.SetRepresentationType('Surface')
    dp_broken.SetRepresentationType('Wireframe')

    # set scalar coloring
    ColorBy(dp, ('CELLS', 'material property ID'))
    ColorBy(dp_broken, ('CELLS', 'failure mode'))

    # get color transfer function/color map for 'materialpropertyID'
    # rescale color and/or opacity maps used to include current data range
    failuremodeLUT = GetColorTransferFunction('failuremode')
    materialpropertyIDLUT = GetColorTransferFunction('material property ID')
    materialpropertyIDLUT.InterpretValuesAsCategories = 1
    # Rescale transfer function
    failuremodeLUT.RescaleTransferFunction(0.0, 3.0)
    dp.RescaleTransferFunctionToDataRange(True)

    # Properties modified
    materialpropertyIDLUT.ColorSpace = 'RGB'
    materialpropertyIDLUT.ApplyPreset('Brewer Qualitative Accent', True)
    materialpropertyIDLUT.Annotations = ['0', '', '1', '', '2', '', '3', '', '4', '', '5', '', '5', '']
    failuremodePWF = GetOpacityTransferFunction('failuremode')
    failuremodePWF.ApplyPreset('Preset 4', True)

    # Properties modified on failuremodeLUT
    failuremodeLUT.ColorSpace = 'Diverging'

    # Properties modified on failuremodeLUT
    failuremodeLUT.RGBPoints = [0, 0, 0, 0, 0.99, 0, 0, 0, 1.0, 0.265, 0.0039, 0.328, 1.49, 0.265, 0.0039, 0.328, 1.5, 0.1289, 0.566, 0.5468, 2.0, 0.1289, 0.566, 0.5468, 2.01, 0.9882, 0.9023, 0.1445, 3.0, 0.9882, 0.9023, 0.1445]

    dp_broken.LineWidth = 3.50
    dp.Opacity = 0.75

    # Render
    Render()

    # Screen shot folder
    Screen = os.path.join(fol, "Screenshots")
    # Create subfolder
    if not os.path.exists(Screen):  # Check to see if the folder exists
        os.makedirs(Screen)  # if not then makes the folder
    # Save Screen shot to folder
    ffname = os.path.join(Screen, (str("%0.4d" % int(fil)) + ".png"))
    WriteImage(ffname, magnification = 1, quality = 100)

    return ffname


'''
Calculate Area of Triangle
# http://en.wikipedia.org/wiki/Heron's_formula/home/aly/Desktop/3D_MESH/mesh_study/mesh1_00
    Inputs:
        - Three points of the mesh element (X, Y, Z)
    Returns:
        - Area of the mesh triangle using Heron's formula
'''


def area(a, b, c):
    def distance(p1, p2):
        return math.hypot(p1[0]-p2[0], p1[1]-p2[1])

    side_a = distance(a, b)
    side_b = distance(b, c)
    side_c = distance(c, a)
    s = 0.5 * ( side_a + side_b + side_c)
    return math.sqrt(s * (s - side_a) * (s - side_b) * (s - side_c))


# /// Main function with command line option handling /// #

def main(argv):
    # Parse arguments from user
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', type=str, required=True,
                        help="Base directory containing ParaView output files")
    parser.add_argument('-b', '--batch', action="store_true",
                        help="Batch process of output files in the subdirectories of the given directory")
    # parser.add_argument('-c', '--center', nargs=3, type=float, dest='specimen_center',
    #                     help="Coordinates of the center of the specimen: x y z")
    parser.add_argument('-e', '--element', type=int, dest='platen_cells_prop_id', default = -1,
                        help="Element material property id of the platen cells")
    parser.add_argument('-p', '--point', nargs=2, type=int, dest='platen_points_prop_id',
                        help="Nodal boundary condition ids of the platen points [top platen, bottom platen]")
    parser.add_argument('--gl', type=float, default=0,
                        help="Strain gauge length in mm")
    parser.add_argument('--gw', type=float, default=0,
                        help="Strain gauge width in mm")
    parser.add_argument('-o', '--output', dest='output_file_name', default="history.csv",
                        help="Name of the output csv file")
    parser.add_argument('-f', '--frequency', dest='bin_frequency',
                        help="Frequency of Outputs")
    args = parser.parse_args()

    # If '.' is specified, the current directory is used as the directory path
    if args.dir == '.':
        args.dir = os.getcwd()

    # Placeholder for the list of directory(ies). Change relative path to absolute path.
    directories = [os.path.abspath(args.dir)]

    # If batch processing, do the analysis on all subdirectories of the given directory
    if args.batch:
        print('\nBatch processing...\n')

        directories = findSubdirectories(os.path.abspath(args.dir))

    # Loop over directory(ies)
    if directories is not None and directories:
        for sub_dir in directories:
            processOutputPropID(sub_dir, "vtu",  args.platen_cells_prop_id, args.platen_points_prop_id, args.output_file_name, args.bin_frequency, args.gl, args.gw)


if __name__ == "__main__":
    try:
        sys.settrace(main(sys.argv[1:]))
    except KeyboardInterrupt:
        # print("\n\033[1;31;0mTERMINATED BY USER\n")
        exit("TERMINATED BY USER")