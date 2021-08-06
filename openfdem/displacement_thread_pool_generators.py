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
    print(strain_from_platen)
    yield strain_from_platen


    def set_strain_gauge(self, gauge_length=None, gauge_width=None):
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
            gauge_width = 0.25 * self.sample_width
        if gauge_length == None:
            gauge_length = 0.25 * self.sample_height

        specimen_center = self.rock_model.center

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
            cv.append(self.first_file.find_closest_cell(pv[ps]))
            ch.append(self.first_file.find_closest_cell(ph[ps]))

        if -1 in pv or -1 in ph:
            print("Check Strain Gauge Dimensions")

        print('\tVertical Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (pv, cv))
        print('\tHorizontal Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (ph, ch))

        gauge_disp_x, gauge_disp_y = [], []

        # Load each timestep
        for openfdem_model_ts in self:
            displacement_y, displacement_x = 0.0, 0.0

            v_strain_gauge = openfdem_model_ts.extract_cells(cv).get_array(self.var_data['gauge_displacement'])
            h_strain_gauge = openfdem_model_ts.extract_cells(ch).get_array(self.var_data['gauge_displacement'])

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

        return gauge_disp_x, gauge_disp_y

if __name__ == '__main__':
    model = fd.Model("../example_outputs/Irazu_UCS")
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
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(history_strain_func, f_names, repeat(model)))  # is self the list we are iterating over
        print(f'the current thread count is: %s' % threading.activeCount())
    print(len(results))
    print(calc_timer_values(time.time() - start))

print(result)

for fname_iter in f_names:
    hist=history_strain_func(fname_iter, model)
    hist.__next__()

