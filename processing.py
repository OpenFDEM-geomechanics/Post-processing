# /////////////////////////////////////////////////////////////// #
# !python3.6
# -*- coding: utf-8 -*-
# Python Script initially created on 2021
# Compiled by Aly @ Grasselli Geomechanics Group, UofT, 2021
# Created using PyCharm
# /////////////////////////////////////////////////////////////// #

import os, re
import collections
import itertools
import pyvista as pv
import formatting_codes

''' 
Current working directory
'''

cwd = os.getcwd()


def findSubdirectories(dir_path):
    '''
    Find all the subdirectories of a given path

    :param dir_path: Starting directory (full path is required)
    :type dir_path: str
    :return: A list of subdirectories
    :rtype: list[str]
    '''

    sub_dirs, listtoremove  = [], []
    for root, dirs, files in os.walk(dir_path):
        for dir_name in dirs:
            sub_dirs.append(os.path.join(root,dir_name))
    sub_dirs = sorted(list(set(sub_dirs))) # Sort directories alphabetically in ascending order
    for i in sub_dirs:
        if 'post_processing' in i:
            listtoremove.append(i)
            listtoremove.append(i.rsplit('/',1)[0])
    print("Found \033[1m%s\033[0m sub-directories" % len(sub_dirs))
    sub_dirs = [x for x in sub_dirs if x not in listtoremove]
    print("After removing already processed file\nFound \033[1m%s\033[0m sub-directories" % len(sub_dirs))
    return sub_dirs


def numericalSort(value):
    '''
    Strip the numerical portion of the file.
    Sort filenames based on natural order (e.g. 1, 2,..., 10, 11, ..., instead of 1, 10, 11, 2, ...)

    :param value: Name of file
    :type value: str
    :return parts: Return the numerical portion of the file
    :rtype: list[str]
    '''
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)  # Split the numerical part of the file
    parts[1::2] = map(int, parts[1::2])  # Return the numerical portion of the file
    return parts


def findOutputFiles(dir_path, file_extension, output_type):
    '''
    Find ParaView output files in a given directory

    :param dir_path: Directory where output files are located (full path is required)
    :param file_extension: File extension of output files ('.vtu'/'.vtp') based on Irazu/openFDEM
    :param output_type: Type of output files, for instance use '_basic_'/'_field_' for the basic/main output files
    :type dir_path: str
    :type file_extension: str
    :type output_type: str
    :return list_of_files: list of output files that match the criteria
    :rtype: list[str]
    '''

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


def processOutputPropID(dir_path):
    '''
     Looks for the various file extensions to decide on processing openFDEM.
     Currently only load BASIC/FIELD VTU/VTP files for Irazu/openFDEM

    :param dir_path: Directory where output files are located
    :type dir_path: str
    :return: list of files based on the defined condition
    :rtype: list
    '''

    ## Checks for IRAZU/openFDEM output files
    global list_of_files
    file_extension, file_name = '.vtp', '_field_'
    list_of_files = findOutputFiles(dir_path, file_extension, file_name)
    if not list_of_files:
        file_extension, file_name = '.vtu', '_basic_'
        list_of_files = findOutputFiles(dir_path, file_extension, file_name)
    print("Processing %s files in %s" % (file_extension, formatting_codes.red_text(dir_path)))
    #  If no output files were found, warn the user and abort
    if list_of_files is None or not list_of_files:
        print(formatting_codes.bold_text("\nERROR: ") +
            "No output files matching extension '%s' containing '%s' were found in %s\n%s\n" % (file_extension, file_name, dir_path, formatting_codes.red_text("ABORTING")))
        return


def model_initial_info(initial_timestep):
    '''
    Obtain information about the model, specifically
        - 2D/3D
        - Identify Irazu/openFDEM output arrays
        - Mineral phase (Mineral IDs)
        - Model/Sample dimensions
        - Type of simulation UCS/BD

    :param initial_timestep: Reads the first time step of the simulation
    :type initial_timestep: str
    :return: rock_elem_ids_max, rock_elem_ids_min, platen_cells_elem_id, width, height
    :rtype: [Tuple[SupportsLessThan, SupportsLessThan, Union[None, vtkStringArray, ndarray, Type[vtkDataArray]], float, float]
    '''

    # rock_elem_ids_max, rock_elem_ids_min, platen_cells_elem_id, width, height
    openfdem_model = pv.read(initial_timestep)

    # Check 2D (3 Points - Triangle) or 3D (4 Points - Tetrahedral) Simulations from the cell vertex.
    number_of_points_per_cell = float(openfdem_model.cell_n_points(0))
    if number_of_points_per_cell == 3: #2D (Triangles)
        print (formatting_codes.green_text("2D Simulation"))
        node_skip = 4
    else: #3D (Tetrahedral)
        print (formatting_codes.green_text("3D Simulation"))
        node_skip = 6
        exit("3D Simulation not supported")

    ''' 
    MODEL PROPERTIES // BOUNDARIES
    '''

    combos = []
    var_dataset = {"openFDEM" : {"mineral_type" : "Property_id",
                           "boundary" : "boundary condition ID",
                            "platen_force" : "force",
                            "platen_displacement" : "displacement",
                           },
                   "IRAZU" : {"mineral_type" : "material property ID",
                            "boundary" : "boundary condition ID",
                            "platen_force" : "force",
                            "platen_displacement" : "displacement",
                            },
                   }
    print("The file contains the following array names\n%s" % openfdem_model.array_names)

    ''' Combinations of the the various material boundaries in the simulation results '''
    # Also returns % composition within the model
    # Create a list of the material types in the model

    global var_data
    try:
        var_data = var_dataset["IRAZU"]
        mineral_bound = openfdem_model.get_array(var_data["mineral_type"])
    except KeyError:
        var_data = var_dataset["openFDEM"]
        mineral_bound = openfdem_model.get_array(var_data["mineral_type"])

    # Return a list of unique material types
    mineral_bound_2 = list(set(mineral_bound))
    const = collections.Counter(mineral_bound)  # counts the frequency of each mineral type
    list_combo = [x for x in itertools.combinations_with_replacement(mineral_bound_2, 2)]  # Creates a list of all possible combinations
    for key, val in enumerate(list_combo):  # Enumerates the listcombo making it "Group #")
        combos.append("Group %d" % key)
        combos.append(val)
    # Combines the Groups, creating a dictionary
    global material_combinations
    material_combinations = dict(
        itertools.zip_longest(*[iter(combos)] * 2, fillvalue=""))
    # Display the various elements in the models
    print(formatting_codes.green_text("Model Element Properties:"))
    for key, val in const.items():
        per = "{0:.2f}".format(float(val) / float(openfdem_model.GetNumberOfCells()) * 100)
        print ("\tElement ID %s - Percentage Composition %s %%" % (formatting_codes.bold_text(key), formatting_codes.bold_text(per)))
    print ("\tPlease be aware of the presence of the Platens in the %'s")

    '''
    Obtain Specimen Center
        - Using rock_model bounds
    '''

    # If Platen Property ID not defined.
    # Lookup cell element ID on the top center and then trace points
    # Using this information, we obtain the platen prop ID

    top_center_point = [openfdem_model.GetCenter()[0], openfdem_model.bounds[3], openfdem_model.bounds[5]]

    top_center_cell = openfdem_model.extract_cells(openfdem_model.find_closest_cell(top_center_point))

    if top_center_cell == -1:
        print(formatting_codes.red_text("Unable to identify Platen ID Correctly. Proceed with caution or use -e flag."))
    platen_cells_elem_id = pv.cell_array(top_center_cell, var_data['mineral_type'])

    print("\tPlaten Material ID found as %s" % formatting_codes.bold_text(platen_cells_elem_id))

    all_elem_id = list(set(openfdem_model.get_array(var_data["mineral_type"])))
    rock_elem_ids = [x for x in all_elem_id if x != platen_cells_elem_id]
    rock_elem_ids_max, rock_elem_ids_min = max(rock_elem_ids), min(rock_elem_ids)

    rock_model = (openfdem_model.threshold([rock_elem_ids_min, rock_elem_ids_max], var_data["mineral_type"]))
    platen = (openfdem_model.threshold([platen_cells_elem_id, platen_cells_elem_id], var_data["mineral_type"]))

    sample_x_min, sample_x_max  = rock_model.bounds[0], rock_model.bounds[1]
    sample_y_min, sample_y_max  = rock_model.bounds[2], rock_model.bounds[3]
    width = sample_x_max - sample_x_min
    height = sample_y_max - sample_y_min
    thickness = rock_model.bounds[4], rock_model.bounds[5]

    check_edge_point = [rock_model.bounds[1], rock_model.bounds[3], 0]
    check_edge_cell = openfdem_model.extract_cells(openfdem_model.find_closest_cell(check_edge_point))
    check_edge_cell = pv.cell_array(check_edge_cell, var_data['mineral_type'])

    if check_edge_cell == -1:
        sim_type = "BD Simulation"
    else:
        sim_type = "UCS Simulation"

    print(sim_type)

    return rock_elem_ids_max, rock_elem_ids_min, platen_cells_elem_id, width, height


def platen_info(pv_cells, platen_boundary_id, var_property):
    '''
    This function thresholds cells based on boundary condition and sums them based on the defined parameter var_property

    :param pv_cells:
    :type pv_cells: pyvista.core.pointset.UnstructuredGrid
    :param platen_boundary_id: boundary id that the threshold should be based on
    :type platen_boundary_id: float
    :param var_property: name of the property (array to b returned)
    :type var_property: Dict[str, str]
    :return: array of the property based on the threshold
    :rtype: ndarray
    '''

    platen_cell_ids = pv_cells.threshold([platen_boundary_id, platen_boundary_id], var_data["boundary"])
    platen_var_prop_list = sum(platen_cell_ids.get_array(var_property))

    if var_property == 'displacement':
        for k in range(0, 3):
            # divide by the number of points per cell  (3 in 2D and 4 in 3D)
            platen_var_prop_list[k] = platen_var_prop_list[k] / (platen_cell_ids.cell_n_points(0) * platen_cell_ids.n_cells)

    return platen_var_prop_list

def UCS_platen_Stress(dir_path):
    '''
    Calculate the stress based on the force in the platens.

    :param dir_path: path of simulation output directory
    :type dir_path: str
    :return: stress (MPa) = force / width of sample
    :rtype: list
    '''

    history_stress = []
    processOutputPropID(dir_path)
    rock_elem_ids_max, rock_elem_ids_min, platen_cells_elem_id, width, height = model_initial_info(list_of_files[0])

    for i in list_of_files:

        openfdem_model_ts = pv.read(i)

        platen = (openfdem_model_ts.threshold([platen_cells_elem_id, platen_cells_elem_id], var_data["mineral_type"]))
        top, bottom = (platen.get_data_range(var_data["boundary"]))

        top_platen_force_list = platen_info(openfdem_model_ts, top, var_data["platen_force"])
        bot_platen_force_list = platen_info(openfdem_model_ts, bottom, var_data["platen_force"])

        # Convert forces from microN to kN and get the average forces & displacements
        avg_platen_force = [0.0, 0.0, 0.0]
        axis_of_loading = 1

        # print(top_platen_force_list, bot_platen_force_list)
        for i in range(0, 3):
            # microN to kN
            avg_platen_force[i] = 0.5 * (abs(top_platen_force_list[i]) + abs(bot_platen_force_list[i])) / 1.0e9

        # stress in MPa (force in kN & area in mm^2)
        stress = avg_platen_force[axis_of_loading] / width * 1.0e3
        history_stress.append(stress)

    return history_stress


def UCS_platen_Strain(dir_path):
    '''
    Calculate the stress based on the strain in the platens.

    :param dir_path: path of simulation output directory
    :type dir_path: str
    :return: strain (%) = displacement / height of sample
    :rtype: list
    '''
    history_strain = []
    processOutputPropID(dir_path)
    rock_elem_ids_max, rock_elem_ids_min, platen_cells_elem_id, width, height = model_initial_info(list_of_files[0])
    for i in list_of_files:
        openfdem_model_ts = pv.read(i)

        platen = (openfdem_model_ts.threshold([platen_cells_elem_id, platen_cells_elem_id], var_data["mineral_type"]))
        top, bottom = (platen.get_data_range(var_data["boundary"]))

        avg_top_platen_disp = platen_info(openfdem_model_ts, top, var_data["platen_displacement"])
        avg_bottom_platen_disp = platen_info(openfdem_model_ts, bottom, var_data["platen_displacement"])

        # Convert average displacements
        avg_platen_disp = [0.0, 0.0, 0.0]
        axis_of_loading = 1

        for i in range(0, 3):
            avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

        strain_from_platen = avg_platen_disp[axis_of_loading] / height * 100.0

        history_strain.append(strain_from_platen)

    return history_strain



