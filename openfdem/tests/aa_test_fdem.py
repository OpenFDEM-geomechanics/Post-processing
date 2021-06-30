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

## Model attributes
model = fd.Model("../../example_outputs/Irazu_UCS")

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
print("By the actual integration time as >>> model['20000']", model['20000'])

print("Get ROCK model dimensions >>> model.rock_sample_dimensions()", model.rock_sample_dimensions())
print("Get the ROCK width from the dimensions", model.sample_width)

print("Get Platen forces >>> force = model.platen_force()")
ax_force = model.platen_force()
print(model.platen_force())

print("Get Platen displacement >>> model.platen_displacement()")
disp = model.platen_displacement()
print(model.platen_displacement())

print("Process UCS >>> model.process_UCS()")
print(model.process_UCS())

print("Get Tan E Modulus >>> model.E[0]")
print(model.E[0])

print("You can also get E over a range using >>> model.E_mod[ax_force, disp, 1, 6]")
print(model.E_mod(ax_force, disp, 0, 1))

