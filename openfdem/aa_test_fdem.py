import sys
sys.path.insert(0,'./')

import openfdem as fd
import time
import os
import matplotlib
import matplotlib.pyplot as plt

# START OF EXECUTION
abs_start = time.time()

my_path = os.path.dirname(
    os.path.abspath(__file__))  # Figures out the absolute path for you in case your working directory moves around.

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

print("Load Model")
# model = fd.Model("../example_outputs/Irazu_UCS")
# model = fd.Model("/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/Irazu_UCS")
# model = fd.Model('/external/Speed_Cal_Using_Flowstone/UCS/UCS_c_17_5_ts_2_55_GII_90000_v_0_8')
model = fd.Model("/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/Irazu_UCS")
# model = fd.Model('/external/Size_7')
df_1 = model.complete_BD_stress_strain(True)
# print(df_1)
exit()
'''
Calculate Modulus
'''
print(model.Etan50_mod(df_1))
print(model.Esec_mod(df_1, 50))
print(model.Esec_mod(df_1, 0.5))
print(model.Eavg_mod(df_1, 50, 60))
print(model.Eavg_mod(df_1, 0.5, 0.6))

print(model.Etan50_mod(df_1, 'Gauge Displacement Y'))
print(model.Esec_mod(df_1, 50, 'Gauge Displacement Y'))
print(model.Esec_mod(df_1, 0.5, 'Gauge Displacement Y'))
print(model.Eavg_mod(df_1, 50, 60, 'Gauge Displacement Y'))
print(model.Eavg_mod(df_1, 0.5, 0.6, 'Gauge Displacement Y'))

'''
PLOTTING FUNCTIONS
'''

ax = df_1.plot('Platen Strain', 'Platen Stress', label="Platen")
df_1.plot('Gauge Displacement Y', 'Platen Stress', ax=ax, label="Strain Gauge Y")
# df_1.plot('Gauge Displacement X', 'Platen Stress', ax=ax, label="Strain Gauge X")
plt.xlabel("Strain")
plt.ylabel("Axial Stress (MPa)")
plt.show()

'''
END DEMO!!!
'''
exit()











## Model attributes
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
df_1 = model.complete_stress_strain(False, 12, 7)

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
