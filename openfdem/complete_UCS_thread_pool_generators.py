import concurrent.futures
import time
from itertools import repeat

import pandas as pd
import pyvista as pv

from . import formatting_codes

# Initialise Variables
history_strain, history_stress = [], []
gauge_disp_x, gauge_disp_y = [], []


def history_strain_func(f_name, model, cv, ch):
    """
    Calculate the axial stress from platens, axial strain from platens and SG as well as lateral strain from SG

    :param f_name: name of vtu file being processed
    :type f_name: str
    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param cv: list of cells at the corner of the vertical strain gauge
    :type cv: list[int]
    :param ch: list of cells at the corner of the horizontal strain gauge
    :type ch: list[int]

    :return: Stress, Platen Strain, SG Strain, SG Lateral Strain
    :rtype: Generator[Tuple[list, list, list, list], Any, None]
    """

    openfdem_model_ts = pv.read(f_name)

    '''STRESS-STRAIN PLATENS'''

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

    # Calculate the stress in MPa (force in kN & area in mm^2)
    stress_from_platen = avg_platen_force[axis_of_loading] / model.sample_width * 1.0e3
    history_stress.append(stress_from_platen)

    # Calculate strains in percentage (%)
    strain_from_platen = avg_platen_disp[axis_of_loading] / model.sample_height * 100.0
    history_strain.append(strain_from_platen)

    '''STRAIN GAUGE ANALYSIS'''

    displacement_y, displacement_x = 0.0, 0.0

    if cv and ch:
        # Extract the data of the cells that cover the extents of the SG
        v_strain_gauge = openfdem_model_ts.extract_cells(cv).get_array(model.var_data['gauge_displacement'])
        h_strain_gauge = openfdem_model_ts.extract_cells(ch).get_array(model.var_data['gauge_displacement'])

        for i in range(0, len(h_strain_gauge)):
            # Vertical contraction is assumed positive
            # Horizontal expansion is assumed positive
            if i < 6:
                # Bottom cells of vertical strain gauge
                # Right cells of horizontal strain gauge
                displacement_y += h_strain_gauge[i][1]
                displacement_x -= v_strain_gauge[i][0]
            else:
                # Top cells of vertical strain gauge
                # Left cells of horizontal strain gauge
                displacement_y -= h_strain_gauge[i][1]
                displacement_x += v_strain_gauge[i][0]

        # Get average displacement in all points.
        displacement_y = displacement_y / 6.0
        displacement_x = displacement_x / 6.0

        # Calculate strains in percentage (%)
        strain_x = displacement_x / g_length * 100.0
        strain_y = displacement_y / g_length * 100.0

        # Append to list
        gauge_disp_x.append(strain_x)
        gauge_disp_y.append(strain_y)

    yield history_stress, history_strain, gauge_disp_x, gauge_disp_y


def set_strain_gauge(model, gauge_length=None, gauge_width=None):
    """
    Calculate local strain based on the dimensions of a virtual strain gauge placed at the center of teh model with
    x/y dimensions. By default set to 0.25 of the length/width.

    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param gauge_length: length of the virtual strain gauge
    :type gauge_length: float
    :param gauge_width: width of the virtual strain gauge
    :type gauge_width: float

    :return: Cells that cover the horizontal and vertical gauges as well as the gauge width and length
    :rtype: [list, list, float, float]
    """

    pv, ph = [], []

    if not gauge_width or gauge_width == 0:
        gauge_width = 0.25 * model.sample_width
    if not gauge_length or gauge_length == 0:
        gauge_length = 0.25 * model.sample_height

    specimen_center = model.rock_model.center

    # Points the constitute the SG
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

    # Cells at the points of the SG
    # These dont change throughout the post-processing
    cv, ch = [], []
    for ps in range(0, len(pv)):
        cv.append(model.first_file.find_closest_cell(pv[ps]))
        ch.append(model.first_file.find_closest_cell(ph[ps]))

    # Verify that SG points are within the domain and return valid cells
    if -1 in pv or -1 in ph:
        print("Check Strain Gauge Dimensions\nWill not process the strain gauges")
        st_status = False
    else:
        print('\tVertical Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (pv, cv))
        print('\tHorizontal Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (ph, ch))

    # Global SG length
    global g_length
    g_length = gauge_length

    return ch, cv, gauge_width, gauge_length


def main(model, platen_id, st_status, gauge_width, gauge_length):
    """
    Main concurrent Thread Pool to calculate the full stress-strain

    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param platen_id: Manual override of Platen ID
    :type platen_id: None or int
    :param st_status: Enable/Disable SG Calculations
    :type st_status: bool
    :param gauge_width: SG width
    :type gauge_width:  float
    :param gauge_length: SG length
    :type gauge_length: float

    :return: full stress-strain information
    :rtype: pd.DataFrame
    """

    # File names of the basic files
    f_names = model._basic_files

    # Get rock dimension.
    model.rock_sample_dimensions(platen_id)

    # Check UCS Simulation
    if model.simulation_type() != "UCS Simulation":
        print("Simulation appears to be not for compressive strength")
        exit("Simulation appears to be not for compressive strength")

    # Global declarations
    start = time.time()

    # Initialise the Strain Gauges
    if st_status:  # Enabled SG st_status == True
        cv, ch, gauge_width, gauge_length = set_strain_gauge(model, gauge_width, gauge_length)
    else:
        cv, ch, gauge_width, gauge_length = [], [], 0, 0

    # Load basic files in the concurrent Thread Pool
    for fname in f_names:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(history_strain_func, fname, repeat(model), cv, ch))  # is self the list we are iterating over

    # Iterate through the files in the defined function
    for fname_iter in f_names:
        hist = history_strain_func(fname_iter, model, cv, ch)
        hist.__next__()

    print(formatting_codes.calc_timer_values(time.time() - start))

    # Merge all into a pandas DataFrame
    if st_status:  # SG Enabled st_status == True
        ucs_df = pd.DataFrame(list(zip(history_stress, history_strain, gauge_disp_x, gauge_disp_y)),
                              columns=['Platen Stress', 'Platen Strain', 'Gauge Displacement X', 'Gauge Displacement Y'])
    else:  # SG Disabled st_status == False
        ucs_df = pd.DataFrame(list(zip(history_stress, history_strain)),
                              columns=['Platen Stress', 'Platen Strain'])
    return ucs_df
