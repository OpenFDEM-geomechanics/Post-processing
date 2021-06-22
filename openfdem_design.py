import openfdem as fd

output_dir = 'example_outputs/Irazu_UCS'
model = fd.Model(output_dir)
# Implementation options:
# 1. Read files immediately and store data as arrays, or
# 2. Store file names and initial geometry data but load files on demand when functions require data
#
# Either way, it determines the following model attributes when initialized

## Model attributes
model.n_timesteps
model.n_points
model.n_elements
model.n_cohesive_elements

model.is_mechanical
model.is_hydraulic
model.is_thermal

model.bounds(mat_id=None)
# bounds of entire model

model.bounds(mat_id=0)
# bounds of material 0

## UCS 
# Assumes UCS model
# Have warnings in place if model does not appear as UCS

# List or array of axial stress at platen for each timestep
axial_stress = model.axial_stress()

# List or array of axial strain of UCS model for each timestep
ucs_strain = model.axial_strain()

# List or array of lateral strain at centre of UCS model for each timestep
lateral_strain = model.lateral_strain()

# Volumetric strain estimated by axial and lateral strain calculation
vol_strain = model.vol_strain()

# Volumetric strain estimated considering all elements
vol_strain_exact = model.vol_strain(exact=True)

# Peak strength value
ucs = model.ucs_strength()

# Peak strength value and associated timestep index
ucs, ucs_index = model.ucs_strength(return_index=True)

# Tangent elastic modulus at 50% strength
e_tan = model.e_tan()

# Tangent elastic modulus at 75% strength
e_tan_75 = model.e_tan(strength_threshold=0.75)

# Secant elastic modulus at 50% strength
e_sec = model.e_sec()

# Average elastic modulus at 50% strength
# How is linear portion estimated?
e_avg = model.e_avg()

## BD
# Functions assume BD
# Produce warnings when model does not appear as BD

# List of peak stress per timestep estimated in model from platen force and disk diameter
bd_stress = model.bd_expected_stress()

# List of peak element stress per timestep in model
bd_actual_stress = model.bd_actual_stress()

# Peak stress calculated from platen force
bd_strength = model.bd_strength()

## Strain gauges
# Can be a general processing tool for any model

# List of strain per timestep
strain_gauge = model.strain_gauge(point1, point2)

