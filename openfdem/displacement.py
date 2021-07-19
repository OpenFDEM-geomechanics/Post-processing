import glob
import os
import os.path as path
import re

import pyvista as pv
import time
import concurrent.futures
from multiprocessing import Process
from threading import Thread
import openfdem as fd
## Model attributes

def calc_timer_values(end_time):
    minutes, sec = divmod(end_time, 60)
    if end_time < 60:
        return ("\033[1m%.2f seconds\033[0m" % end_time)
    else:
        return ("\033[1m%d minutes and %d seconds\033[0m." % (minutes, sec))

start = time.time()
# model = fd.Model("../example_outputs/Irazu_UCS")
model = fd.Model("/external/Size_7")
f_names = (model._basic_files)

history_strain = []  # List to hold the stress_history
avg_platen_disp = [0.0, 0.0, 0.0]  # Dummy cell
axis_of_loading = 1  # Axis of loading in Y direction.

## Get rock dimension.
model.rock_sample_dimensions()
## Check UCS Simulation
if model.simulation_type() != "UCS Simulation":
    print("Simulation appears to be not for compressive strength")

# Load each timestep
# backup copy compressed below

def history_strain_func(f_name):
    openfdem_model_ts = pv.read(f_name)
    platen = (openfdem_model_ts.threshold([model.platen_cells_elem_id, model.platen_cells_elem_id],
                                          model.var_data["mineral_type"]))
    top, bottom = (platen.get_data_range(model.var_data["boundary"]))

    avg_top_platen_disp = model.platen_info(openfdem_model_ts, top, model.var_data["platen_displacement"])
    avg_bottom_platen_disp = model.platen_info(openfdem_model_ts, bottom,
                                              model.var_data["platen_displacement"])

    for i in range(0, model.number_of_points_per_cell):
        avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

    strain_from_platen = avg_platen_disp[axis_of_loading] / model.sample_height * 100.0
    history_strain.append(strain_from_platen)

    return strain_from_platen

def main():
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(history_strain_func, f_names)  # is self the list we are iterating over
    return results

if __name__ == '__main__':
    main()
    print(calc_timer_values(time.time() - start))


# platen_displacement(model)
