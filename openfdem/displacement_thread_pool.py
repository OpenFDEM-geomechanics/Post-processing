import glob
import os
import os.path as path
import re
from itertools import repeat

import pyvista as pv
import time
import concurrent.futures
from multiprocessing import Process
import threading
import openfdem as fd


## Model attributes

def calc_timer_values(end_time):
    minutes, sec = divmod(end_time, 60)
    if end_time < 60:
        return ("\033[1m%.2f seconds\033[0m" % end_time)
    else:
        return ("\033[1m%d minutes and %d seconds\033[0m." % (minutes, sec))


# Load each timestep
# backup copy compressed below

def history_strain_func_hdf(f_name, model):
    openfdem_model_ts = model.storage.read_file(f_name)
    platen = (openfdem_model_ts.threshold([model.platen_cells_elem_id, model.platen_cells_elem_id],
                                          model.var_data["mineral_type"]))
    top, bottom = (platen.get_data_range(model.var_data["boundary"]))

    avg_top_platen_disp = model.platen_info(openfdem_model_ts, top, model.var_data["platen_displacement"])
    avg_bottom_platen_disp = model.platen_info(openfdem_model_ts, bottom,
                                               model.var_data["platen_displacement"])

    avg_platen_disp = [0.0, 0.0, 0.0]  # Dummy cell
    axis_of_loading = 1  # Axis of loading in Y direction.

    for i in range(0, model.number_of_points_per_cell):
        avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

    strain_from_platen = avg_platen_disp[axis_of_loading] / model.sample_height * 100.0

    return strain_from_platen


def history_strain_func(f_name, model):
    openfdem_model_ts = pv.read(f_name)
    platen = (openfdem_model_ts.threshold([model.platen_cells_elem_id, model.platen_cells_elem_id],
                                          model.var_data["mineral_type"]))
    top, bottom = (platen.get_data_range(model.var_data["boundary"]))

    avg_top_platen_disp = model.platen_info(openfdem_model_ts, top, model.var_data["platen_displacement"])
    avg_bottom_platen_disp = model.platen_info(openfdem_model_ts, bottom,
                                               model.var_data["platen_displacement"])

    avg_platen_disp = [0.0, 0.0, 0.0]  # Dummy cell
    axis_of_loading = 1  # Axis of loading in Y direction.

    for i in range(0, model.number_of_points_per_cell):
        avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

    strain_from_platen = avg_platen_disp[axis_of_loading] / model.sample_height * 100.0

    return strain_from_platen


if __name__ == '__main__':

    model = fd.Model("../example_outputs/Irazu_GBM")
    f_names = (model._basic_files)

    ## Get rock dimension.
    model.rock_sample_dimensions()
    ## Check UCS Simulation
    if model.simulation_type() != "UCS Simulation":
        print("Simulation appears to be not for compressive strength")
    print("---- PV tests ----")
    print("Start for loop")
    start = time.time()
    result = []
    for fname in f_names:
        result.append(history_strain_func(fname, model))
    print(len(result))
    print(calc_timer_values(time.time() - start))

    print(f"Start concurrent")
    start = time.time()
    max=100
    with concurrent.futures.ThreadPoolExecutor(max) as executor:
        results = list(executor.map(history_strain_func, f_names, repeat(model),
                                    chunksize=len(f_names) //max))  # is self the list we are iterating over
        print(f'the current thread count is: %s' % threading.activeCount())
    print(len(results))
    print(calc_timer_values(time.time() - start))

    print(f"Start concurrent optimized n")
    for i in range(1, 10):
        n = 2 ** i
        print(f"n = {n}, chunksize = {len(f_names) // n}")
        start = time.time()
        with concurrent.futures.ThreadPoolExecutor(n) as executor:
            results = list(executor.map(history_strain_func, f_names, repeat(model),
                                        chunksize=len(f_names) // n))  # is self the list we are iterating over
            print(f'the current thread count is: %s' % threading.activeCount())
        print(len(results))
        print(calc_timer_values(time.time() - start))

    print("chunksize = 1")
    start = time.time()
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(
            executor.map(history_strain_func, f_names, repeat(model)))  # is self the list we are iterating over
    print(len(results))
    print(calc_timer_values(time.time() - start))

    ##--------------------------------------------------------------
    print("---- HDF tests ----")
    print("Start for loop")
    start = time.time()
    result = []
    for fname in f_names:
        result.append(history_strain_func_hdf(fname, model))
    print(len(result))
    print(calc_timer_values(time.time() - start))

    print(f"Start concurrent")
    start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max) as executor:
        results = list(executor.map(history_strain_func_hdf, f_names, repeat(model),
                                    chunksize=len(f_names) // max))  # is self the list we are iterating over
        print(f'the current thread count is: %s' % threading.activeCount())
    print(len(results))
    print(calc_timer_values(time.time() - start))

    print(f"Start concurrent optimized n")
    for i in range(1, 10):
        n = 2 ** i
        print(f"n = {n}, chunksize = {len(f_names) // n}")
        start = time.time()
        with concurrent.futures.ThreadPoolExecutor(n) as executor:
            results = list(executor.map(history_strain_func_hdf, f_names, repeat(model),
                                        chunksize=len(f_names) // n))  # is self the list we are iterating over
            print(f'the current threads are = %s' % threading.activeCount())
        print(len(results))
        print(len(f_names) // n)
        print(calc_timer_values(time.time() - start))

    print("chunksize = 1")
    start = time.time()
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(
            executor.map(history_strain_func_hdf, f_names, repeat(model)))  # is self the list we are iterating over
    print(len(results))
    print(calc_timer_values(time.time() - start))

# platen_displacement(model)
