# /////////////////////////////////////////////////////////////// #
# Python Script initially created on 2022-11-30
# Compiled by Aly @ Grasselli Geomechanics Group, UofT, 2022
# Created using PyCharm
# /////////////////////////////////////////////////////////////// #

import openfdem as fd
import time
import os
import matplotlib
import matplotlib.pyplot as plt
import windrose
import math
import numpy
import random, pandas

# START OF EXECUTION
abs_start = time.time()

'''
Default MATPLOTLIB Fonts
'''

# plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams["figure.figsize"] = [5, 5]
matplotlib.rcParams['font.family'] = ['arial']
matplotlib.rcParams['font.size'] = 8

'''
DEMO!!!
'''
ds = False
rose = False
additional = False
modulus = False
plot_SS = False
Model_attributes = False
extract_along_line = False
dim3_model = False

plt = True

'''
3D Models
'''
if plt:
    # print("LOAD 3D Models")
    # for i in ['Run2', 'Run1']:
    #     path = '/external/OUTPUTS/v6_real_dims_0mm_shift'
    #     model = fd.Model(os.path.join(path, i))
    #     # Extract all cells that meet the criteria and split to nodewise data for each time step.
    #     # In this case "BOUNDARY CONDITION" is set to "1" for the threshold with the "FORCE" being extracted at each node.
    #     df = model.extract_threshold_info(thres_id=1, thres_array='boundary', arrays_needed=['platen_force', 'platen_displacement'])
    #     # Sum the X,Y,Z of all nodes for each time step.
    #     df_sum = model.convert_to_xyz_array(df)
    #
    #     # Handle data like a typical DataFrame
    #     ax = df_sum.plot(title = i)
    #     ax.set(xlabel='Time Step', ylabel='Sigma Force')
    #     plt.show()

    import openfdem as fdem
    import matplotlib.pyplot as plt
    data = fdem.Model("/external/Yusuf_PLT/Axial")

    df = data.complete_PLT_stress_strain(load_config="A", platen_id=1)

'''
3D Models
'''
if dim3_model:
    # print("LOAD 3D Models")
    # for i in ['Run2', 'Run1']:
    #     path = '/external/OUTPUTS/v6_real_dims_0mm_shift'
    #     model = fd.Model(os.path.join(path, i))
    #     # Extract all cells that meet the criteria and split to nodewise data for each time step.
    #     # In this case "BOUNDARY CONDITION" is set to "1" for the threshold with the "FORCE" being extracted at each node.
    #     df = model.extract_threshold_info(thres_id=1, thres_array='boundary', arrays_needed=['platen_force', 'platen_displacement'])
    #     # Sum the X,Y,Z of all nodes for each time step.
    #     df_sum = model.convert_to_xyz_array(df)
    #
    #     # Handle data like a typical DataFrame
    #     ax = df_sum.plot(title = i)
    #     ax.set(xlabel='Time Step', ylabel='Sigma Force')
    #     plt.show()

    import openfdem as fdem
    import matplotlib.pyplot as plt
    data = fdem.Model("../example_outputs/Irazu_3D_UCS")

    df = data.complete_UCS_stress_strain(None, False)

    # # # Extract all cells that meet the criteria and split to nodewise data for each time step.
    # # In this case "BOUNDARY CONDITION" is set to "1" for the threshold with the "FORCE" being extracted at each node.
    # df = data.extract_threshold_info(thres_id=1, thres_array='boundary', arrays_needed=['platen_force'])
    # # Sum the X,Y,Z of all nodes for each time step.
    # df_sum = data.convert_to_xyz_array(df)
    # print(df_sum)
    #
    #
    # df = data.extract_threshold_info(thres_id=1, thres_array='boundary', arrays_needed=['platen_force', 'platen_displacement'])
    # # Sum the X,Y,Z of all nodes for each time step.
    # for k,v in df.items():
    #     print(k)
    #     df_sum = data.convert_to_xyz_array(v)
    #     print(df_sum)



'''
DS EXAMPLE
'''

if ds:
    print("LOAD DS")
    model = fd.Model("/external/2D_shear_4mm_profile_normal_load_test")
    model.direct_shear_calculation(1, 'platen_force', progress_bar=True)
    # exit('DS Example')

'''
Draw Rosette Diagram
'''
if rose:
    model = fd.Model("/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/Irazu_UCS")
    model.draw_rose_diagram(t_step=0)
    model.draw_rose_diagram(t_step=5)


    for i in range(0, model.n_timesteps):
        model.draw_rose_diagram(t_step=i, thres_id=0)

    # Draw a pre-defined rosette
    rand_list_ang=[]
    randomlist = []
    n=200
    for i in range(n):
        rand_list_ang.append(random.randint(0,180))
        randomlist.append(random.randint(1, 5))

    ver_data = pandas.DataFrame(list(zip(randomlist, rand_list_ang)),
                   columns =['Mode I', 'Angle'])

    # figure_name = model.draw_rose_diagram(t_step=0, rose_data=ver_data, range='Mode I')
    figure_name = model.draw_rose_diagram(t_step=0, rose_range='Length', thres_id=0, thres_array='boundary')
    figure_name.tight_layout()
    # figure_name.show()
    figure_name.savefig('/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/example.pdf')
    # exit('Rosette')

"""
Additional Commands
"""
if additional:
    model_ucs = fd.Model("/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/Irazu_UCS")
    print(model_ucs.model_dimensions(0))
    print(model_ucs.model_dimensions(1))
    print(model_ucs.rock_sample_dimensions()[3])
    df_1 = model_ucs.complete_UCS_stress_strain(platen_id=None, st_status=True, c_center=[1, 1, 0])

    df_2 = model_ucs.complete_UCS_stress_strain(platen_id=None, st_status=False, c_center=[1, 1, 0])
    print(len(df_1))
    print(len(df_2))

    model = fd.Model('/external/Speed_Cal_Using_Flowstone/UCS/UCS_c_17_5_ts_2_55_GII_90000_v_0_8')
    print(model.model_dimensions(0))
    print(model.model_dimensions(1))
    print(model_ucs.rock_sample_dimensions()[3])
    df_bd = model.complete_UCS_stress_strain(platen_id=None, st_status=True, c_center=[25, 15, 0])
    # exit('Additional Comments')



'''
Calculate Modulus
'''
if modulus:
    model = fd.Model(
        "/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/Irazu_UCS")
    df_1 = model.complete_UCS_stress_strain(st_status=True)
    print(df_1)
    print(df_1.keys())

    print(model.Etan50_mod(df_1)[0])
    print(model.Etan50_mod(df_1, linear_bestfit=False)[0])
    print(model.Etan50_mod(df_1, loc_strain='Gauge Displacement Y', plusminus_range=1))
    print(model.Esec_mod(df_1, 70))
    print(model.Esec_mod(df_1, 0.5))
    # print(model.Eavg_mod(df_1, 20, 30))
    print(model.Eavg_mod(df_1, 0.5, 0.6))
    print(model.Eavg_mod(df_1, 0.5, 0.6, linear_bestfit=False))
    # exit('Modulus')

'''
Extract along line
'''
if extract_along_line:
    model = fd.Model(
        "/external/HC_160000", 'vtp')
    for i in range (0, 11):
        cellid = (model.find_cell([0, i, 0]))
        print(cellid)
        # cellinfo = model.extract_cell_info(cellid, 'mineral_type')
    #     cellinfo2 = model.extract_cell_info(cellid, 'platen_force')
    #     print(cellid, cellinfo, cellinfo2)

    # stress_cellinfo = model.extract_cell_info(2146, 'platen_force')
    # cellinfo2 = model.extract_cell_info(76, 'platen_force')

    stress_cellinfo4 = model.extract_cell_info(cell_id=1683, arrays_needed=['temperature'], progress_bar=True)
    # print(stress_cellinfo)
    print(stress_cellinfo4)
    # import formatting_codes
    # print(formatting_codes.red_text(val="str"))
    #
    # print(model.rotary_shear_calculation(0))
    # exit()

    # print(stress_cellinfo)
    # print(cellinfo2)
    # print(stress_cellinfo3)
    # import pandas as pd
    # stress_cellinfo3.to_csv('/home/aly/Desktop/trial1.csv')
    # print(stress_cellinfo3.dtypes)
    # stress_cellinfo3['mineral_type'] = stress_cellinfo3['mineral_type'].astype(float)
    # print(stress_cellinfo3['mineral_type'][0][0])

    # stress_cellinfo4 = pd.DataFrame()
    # for idx, i in enumerate(list(stress_cellinfo3.columns)):
    #     name_list = []
    #     print(len(stress_cellinfo3[i][0]))
    #     if len(stress_cellinfo3[i][0]) == 1:
    #         temp_list = stress_cellinfo3.explode(i)
    #         df = temp_list[i]
    #     else:
    #         for j in range(1, len(stress_cellinfo3[i][0]) + 1):
    #             name_list.append(i + "_N%s" % j)
    #         print(name_list)
    #         temp_list = stress_cellinfo3[i].to_list()
    #         print(temp_list)
    #
    #         df = pd.DataFrame([pd.Series(x) for x in stress_cellinfo3[i]])
    #         df.columns = name_list
    #
    #     stress_cellinfo4 = pd.concat([stress_cellinfo4, df], axis=1)

    # print(stress_cellinfo3)
    print(stress_cellinfo4)
    print(stress_cellinfo4.dtypes)
    # plt.plot(stress_cellinfo4['mineral_type'])
    # plt.plot(stress_cellinfo4['boundary_N1'])

    '''
    PLOTTING METHOD ONE
    '''
    x, y = [], []
    for i, row in stress_cellinfo4.iterrows():
        x.append(i)
        y.append(row['temperature_N2'])
    plt.plot(x, y, c='red', label='temperature_N2')
    plt.legend()
    plt.show()

    '''
    PLOTTING METHOD TWO
    '''

    model = fd.Model(
        "/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/Irazu_UCS")

    stress_cellinfo4 = model.extract_cell_info(cell_id=1683, arrays_needed=['platen_force'], progress_bar=True)

    lx = stress_cellinfo4['platen_force_N2'].to_list()
    lx1 = list(zip(*lx))
    plt.plot(lx1[0], label='platen_force_N2_x')
    plt.plot(lx1[1], label='platen_force_N2_y')
    plt.plot(lx1[2], label='platen_force_N2_z')

    plt.legend()
    plt.show()
    stress_cellinfo4.to_csv('/home/aly/Desktop/trial2.csv')


'''
PLOTTING FUNCTIONS
'''
if plot_SS:
    model = fd.Model(
        "/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/Irazu_UCS")
    df_1 = model.complete_UCS_stress_strain(st_status=True)

    # ax = df_1.plot('Platen Strain', 'Platen Stress', label="Platen")
    # df_1.plot('Gauge Displacement Y', 'Platen Stress', ax=ax, label="Strain Gauge Y")
    # # df_1.plot('Gauge Displacement X', 'Platen Stress', ax=ax, label="Strain Gauge X")
    # plt.xlabel("Strain")
    # plt.ylabel("Axial Stress (MPa)")
    # plt.show()

    plt.figure(figsize=(5, 5))
    # model.plot_stress_strain(df_1[''], df_1['Platen Stress'], label='SG', color='red')
    # model.plot_stress_strain(df_1['Platen Strain'], df_1['Platen Stress'], label='IRAZU UCS', color='green')
    # model.plot_stress_strain(df_1['Gauge Displacement X'], df_1['Platen Stress'], label='SG-X', color='blue')

    plots = {'Gauge Displacement Y': 'red', 'Gauge Displacement X':'blue', 'Platen Strain':'green'}
    for k,v  in plots.items():
        model.plot_stress_strain(df_1[k], df_1['Platen Stress'], color=v)

    plt.title('Stress-Strain Curve')
    plt.xlim(-1, 0.25)
    plt.legend()
    plt.show()


    # exit('Stress-Strain Plot')

if Model_attributes:
    print("Load Model")
    model = fd.Model("../example_outputs/Irazu_UCS")
    # model = fd.Model("/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/Irazu_UCS")
    # model = fd.Model('/external/Speed_Cal_Using_Flowstone/UCS/UCS_c_17_5_ts_2_55_GII_90000_v_0_8')
    # model = fd.Model('/external/Size_7')

    # Get Model information
    print("Number of Timesteps:\t", model.n_timesteps)
    print("Number of Points in model:\t", model.n_points)
    print("Number of Elements in model:\t", model.n_elements)
    print("Engine Type:\t", model._fdem_engine)

    print("Model is 2D/3D?:\t", model.model_domain() / 2)

    print("Get model dimensions as a tuple3 >>> model.model_dimensions()", model.model_dimensions())
    print("Get the width from the dimensions after running dimensions by >>> model.model_width", model.model_width)

    print("If you want to get the extents of a material based on the material id (0) >>> model.model_dimensions(mat_id=0)", model.model_dimensions(mat_id=0))
    print("If you want to get the extents of a material based on the material id (1) >>> model.model_dimensions(mat_id=1)", model.model_dimensions(mat_id=1))

    print("To load a certain file you can do it either by entering the output timestep as an integer as >>> model[X]", model[1])
    print("OR")
    # print("By the actual integration time as >>> model['20000']", model['20000'])

    print("Get ROCK model dimensions >>> model.rock_sample_dimensions()", model.rock_sample_dimensions())
    print("Get the ROCK width from the dimensions", model.sample_width)

    print("Get Platen forces >>> force = model.platen_force()")
    # ax_force = model.platen_force()
    # print(ax_force)
    # print(model.my_fun())
    # print("Get Platen displacement >>> model.platen_displacement()")
    # disp = model.platen_displacement()
    # print(model.platen_displacement())
    #
    # print("Process UCS >>> model.process_UCS()")
    # print(model.process_UCS())
    #
    # print("Get Tan E Modulus >>> model.E[0]")
    # print(model.E[0])
    #
    # print("You can also get E over a range using >>> model.E_mod[ax_force, disp, 1, 6]")
    # print(model.E_mod(ax_force, disp, 0, 1))
    print('-----')
    df_1 = model.complete_UCS_stress_strain(False, 12, 7)

    import pandas as pd
    # print(df_1)
    # print(type(df_1))
    # print(df_1['Platen Strain'])
    # ax = plt.figure()
    ax = df_1.plot('Platen Strain', 'Platen Stress', label="Platen")
    # df_1.plot('Gauge Displacement Y', 'Platen Stress', ax=ax, label="Strain Gauge Y")
    # df_1.plot('Gauge Displacement X', 'Platen Stress', ax=ax, label="Strain Gauge X")
    plt.xlabel("Strain")
    plt.ylabel("Axial Stress (MPa)")
    plt.show()