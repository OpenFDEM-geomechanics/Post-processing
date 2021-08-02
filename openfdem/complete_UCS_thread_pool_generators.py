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
import pandas as pd
## Model attributes

def calc_timer_values(end_time):
    minutes, sec = divmod(end_time, 60)
    if end_time < 60:
        return ("\033[1m%.2f seconds\033[0m" % end_time)
    else:
        return ("\033[1m%d minutes and %d seconds\033[0m." % (minutes, sec))


history_strain, history_stress = [], []
gauge_disp_x, gauge_disp_y = [], []

def history_strain_func(f_name, model, cv, ch):
    openfdem_model_ts = pv.read(f_name)
    platen = (openfdem_model_ts.threshold([model.platen_cells_elem_id, model.platen_cells_elem_id],
                                          model.var_data["mineral_type"]))
    top, bottom = (platen.get_data_range(model.var_data["boundary"]))

    top_platen_force_list = model.platen_info(openfdem_model_ts, top, model.var_data["platen_force"])
    bot_platen_force_list = model.platen_info(openfdem_model_ts, bottom, model.var_data["platen_force"])

    avg_top_platen_disp = model.platen_info(openfdem_model_ts, top, model.var_data["platen_displacement"])
    avg_bottom_platen_disp = model.platen_info(openfdem_model_ts, bottom,
                                               model.var_data["platen_displacement"])

    avg_platen_disp = [0.0, 0.0, 0.0]  # Dummy cell
    avg_platen_force = [0.0, 0.0, 0.0]  # Dummy cell
    axis_of_loading = 1  # Axis of loading in Y direction.

    for i in range(0, model.number_of_points_per_cell):
        # Convert forces from microN to kN and get the average forces
        avg_platen_force[i] = 0.5 * (abs(top_platen_force_list[i]) + abs(bot_platen_force_list[i])) / 1.0e9
        avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

    # stress in MPa (force in kN & area in mm^2)
    stress = avg_platen_force[axis_of_loading] / model.sample_width * 1.0e3
    history_stress.append(stress)

    strain_from_platen = avg_platen_disp[axis_of_loading] / model.sample_height * 100.0
    history_strain.append(stress)

    displacement_y, displacement_x = 0.0, 0.0

    v_strain_gauge = openfdem_model_ts.extract_cells(cv).get_array(model.var_data['gauge_displacement'])
    h_strain_gauge = openfdem_model_ts.extract_cells(ch).get_array(model.var_data['gauge_displacement'])

    for i in range(0, len(h_strain_gauge)):
        # Vertical contraction is assumed positive
        # Horizontal expansion is assumed positive
        if i < 6:
            # Bottom cells of vertical strain gauge
            # Right cells of horizontal strain gauge
            displacement_y += v_strain_gauge[i][1]
            displacement_x -= h_strain_gauge[i][0]
        else:
            # Top cells of vertical strain gauge
            # Left cells of horizontal strain gauge
            displacement_y -= v_strain_gauge[i][1]
            displacement_x += h_strain_gauge[i][0]

    displacement_y = displacement_y / 6.0
    displacement_x = displacement_x / 6.0

    # Calculate strains in percentage (%)
    strain_x = displacement_x / gauge_length * 100.0
    strain_y = displacement_y / gauge_length * 100.0
    # print('strain_x & _y at time step', strain_x, strain_y)
    gauge_disp_x.append(strain_x)
    gauge_disp_y.append(strain_y)


    yield history_stress, history_strain, gauge_disp_x, gauge_disp_y


def set_strain_gauge(gauge_length=None, gauge_width=None):
    '''
    Calculate local strain based on the dimensions of a virtual strain gauge placed at the center of teh model with x/y dimnesions. By default set to 0.25 of the length/width.

    :param gauge_length: length of the virtual strain gauge
    :type gauge_length: float
    :param gauge_width: width of the virtual strain gauge
    :type gauge_width: float
    :return:
    '''

    pv, ph = [], []

    if gauge_width==None:
        gauge_width = 0.25 * model.sample_width
    if gauge_length == None:
        gauge_length = 0.25 * model.sample_height


    specimen_center = model.rock_model.center

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
    for ps in range(0, len(pv)):
        cv.append(model.first_file.find_closest_cell(pv[ps]))
        ch.append(model.first_file.find_closest_cell(ph[ps]))

    if -1 in pv or -1 in ph:
        print("Check Strain Gauge Dimensions")

    print('\tVertical Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (pv, cv))
    print('\tHorizontal Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (ph, ch))

    return ch, cv, gauge_width, gauge_length


if __name__ == '__main__':
    # model = fd.Model("../example_outputs/Irazu_UCS")
    model = fd.Model("/external/Size_7")
    f_names = (model._basic_files)

    ## Get rock dimension.
    model.rock_sample_dimensions()
    ## Check UCS Simulation
    if model.simulation_type() != "UCS Simulation":
        print("Simulation appears to be not for compressive strength")
    print("---- PV tests ----")
    print("Start for loop")
    start = time.time()
    global gauge_width, gauge_length
    cv, ch, gauge_width, gauge_length = set_strain_gauge()

for fname in f_names:
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(history_strain_func, f_names, repeat(model), cv, ch))  # is self the list we are iterating over


for fname_iter in f_names:
    hist=history_strain_func(fname_iter, model, cv, ch)
    hist.__next__()

print(calc_timer_values(time.time() - start))

UCS_df = pd.DataFrame(list(zip(history_stress, history_strain, gauge_disp_x, gauge_disp_y)),
                      columns=['Platen Stress', 'Platen Strain', 'Gauge Displacement X', 'Gauge Displacement Y'])
print(UCS_df)