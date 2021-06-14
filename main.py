# /////////////////////////////////////////////////////////////// #
# !python3.6
# -*- coding: utf-8 -*-
# Python Script initially created on 2021
# Compiled by Aly @ Grasselli Geomechanics Group, UofT, 2021
# Created using PyCharm
# /////////////////////////////////////////////////////////////// #

import platform
import sys
import os

if platform.architecture()[0] != "64bit":
    exit("Compatible only on x64")

# Check paraview Version before proceeding
try:
    assert sys.version_info >= (3, 5) and sys.version_info <= (3, 8)
    print("Python Version: %s" % sys.version.split('\n')[0])
except AssertionError:
    exit("Compatible Python Version 3.5+ upto 3.8.x")

import distro
import collections
import itertools
import numpy as np
import pyvista as pv

cwd = os.getcwd()
print(cwd)
import formatting_codes

print("Computer Running %s, Release %s, Version %s" % (platform.system(), platform.release(), platform.version()))

# Check paraview Version before proceeding
try:
    assert sys.version_info >= (3, 5) and sys.version_info <= (3, 8)
    print("Python Version: %s" % sys.version.split('\n')[0])
except AssertionError:
    exit("Compatible Python Version 3.5+ upto 3.8.x")


# def findCenter(threshold_data):
#     x_center = (threshold_data.bounds[0] + threshold_data.bounds[1]) / 2
#     y_center = (threshold_data.bounds[2] + threshold_data.bounds[3]) / 2
#     z_center = (threshold_data.bounds[4] + threshold_data.bounds[5]) / 2
#     threshold_data_bounds = [x_center, y_center, z_center]
#     return (threshold_data_bounds)

def load_file(filename):
    openfdem_model = pv.read(filename)

    # Check 2D (3 Points - Triangle) or 3D (4 Points - Tetrahedral) Simulations from the cell vertex.

    number_of_points_per_cell = float(openfdem_model.cell_n_points(0))
    if number_of_points_per_cell == 3: #2D (Triangles)
        print (formatting_codes.green_text("2D Simulation"))
        node_skip = 4
    else: #3D (Tetrahedral)
        print (formatting_codes.green_text("3D Simulation"))
        node_skip = 6
        exit("3D Simulation not supported")

    print(formatting_codes.green_text("Starting Initialization"))

    ''' 
    MODEL PROPERTIES // BOUNDARIES
    '''

    combos = []
    var_dataset = {"YGEO":{"mineral_type":"Property_id",
                           },
                   "IRAZU":{"mineral_type":"material property ID"
                            },
                   }
    print("The file contains the following array names\n%s" % openfdem_model.array_names)

    ''' Combinations of the the various material boundaries in the simulation results '''
    # Also returns % composition within the model
    # Create a list of the material types in the model

    try:
        var_data = var_dataset["IRAZU"]
        mineral_bound = openfdem_model.get_array(var_data["mineral_type"])
    except KeyError:
        var_data = var_dataset["YGEO"]
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

    # model_center = findCenter(openfdem_model)
    # specimen_center = findCenter(rock_model)

    # print(specimen_center)
    # print(rock_model)
    # print(rock_model.extract_all_edges())
    # print()

    # If Platen Property ID not defined.
    # Lookup cell element ID on the top center and then trace points
    # Using this information, we obtain the platen prop ID

    top_center_point = [openfdem_model.GetCenter()[0], openfdem_model.bounds[3], openfdem_model.bounds[5]]

    top_center_cell = openfdem_model.extract_cells(openfdem_model.find_closest_cell(top_center_point))

    if top_center_cell == -1:
        print(formatting_codes.red_text("Unable to identify Platen ID Correctly. Proceed with caution or use -e flag."))
    platen_cells_elem_id = pv.cell_array(top_center_cell, var_data['mineral_type'])

    print("\tPlaten Material ID found as %s" % formatting_codes.bold_text(platen_cells_elem_id))

    # platen_cells_elem_id = 3
    all_elem_id = list(set(openfdem_model.get_array(var_data["mineral_type"])))
    rock_elem_ids = [x for x in all_elem_id if x != platen_cells_elem_id]
    rock_elem_ids_max, rock_elem_ids_min = max(rock_elem_ids), min(rock_elem_ids)

    rock_model = (openfdem_model.threshold([rock_elem_ids_min, rock_elem_ids_max], var_data["mineral_type"]))


    sample_x_min, sample_x_max  = rock_model.bounds[0], rock_model.bounds[1]
    sample_y_min, sample_y_max  = rock_model.bounds[2], rock_model.bounds[3]
    width = sample_x_max - sample_x_min
    height = sample_y_max - sample_y_min
    thickness = rock_model.bounds[4], rock_model.bounds[5]

    check_edge_point = [rock_model.bounds[1], rock_model.bounds[3], 0]
    print(openfdem_model.find_closest_cell(check_edge_point))
    check_edge_cell = openfdem_model.extract_cells(openfdem_model.find_closest_cell(check_edge_point))
    check_edge_cell = pv.cell_array(check_edge_cell, var_data['mineral_type'])
    print(check_edge_cell)

    if check_edge_cell == -1:
        sim_type = "BD Simulation"
    else:
        sim_type = "UCS Simulation"

    print(sim_type)
    print(openfdem_model.volume)
    exit()
    #
    # mesh = pv.read(filename)
    #
    # print(mesh.get_array('material property ID'))
    # # a = mesh.get_array('material property ID')
    # # print(set(list(a)))
    # print(mesh.bounds)
    # cops = mesh.plot(1)
    # print(mesh.get_data_range('material property ID'))
    # print((mesh.threshold([0,2], 'material property ID')).bounds)
    # exit()
    # for i in range(0,4):
    #     print(i)
    #     platen = (mesh.threshold([i , i], 'material property ID'))
    #     print(platen)
    #     a = platen.get_array('displacement')
    #     print(len(a))
    #     print(platen.bounds)
    #     print(platen.get_data_range('displacement'))
    #     print(a.min(), a.max())
    # print(mesh)
    # # for i in mesh.array_names:
    # #     cpos = mesh.plot(i)
    # #     print(cpos)
    # # cpos = mesh.plot(1)
    # # print(cpos)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # load_file(os.path.join(cwd, "examples", "ygeo", "BD-GUI_field_0.vtp"))
    load_file(os.path.join(cwd, "examples", "Irazu", "01_25red_20170718_10strength_25modes_multi_phase_0MPa_femdem_multi_phase_femdem.r2m_basic_0.vtu"))
