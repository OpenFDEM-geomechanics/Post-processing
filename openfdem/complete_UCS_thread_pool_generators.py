import concurrent.futures
import math
import time
from itertools import repeat

import pandas as pd
import pyvista as pv

try:
    from . import formatting_codes
except ImportError:
    import formatting_codes

loading_dir_dict = {0: "X", 1: "Y", 2: "Z"}

def history_strain_func(f_name, model, cv, ch, axis):
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

    avg_platen_disp = [0.0, 0.0, 0.0] # Dummy cell
    avg_platen_force = [0.0, 0.0, 0.0] # Dummy cell
    load_axis = axis[0]  # Axis of loading in Y direction.

    for i in range(0, 3):
        # Convert forces from microN to kN and get the average forces
        avg_platen_force[i] = 0.5 * (abs(top_platen_force_list[i]) + abs(bot_platen_force_list[i])) / 1.0e9
        avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

    # Calculate the stress in MPa (force in kN & area in mm^2)
    stress_from_platen = avg_platen_force[load_axis] / samp_A * 1.0e3
    history_stress.append(stress_from_platen)

    # Calculate strains in percentage (%)
    strain_from_platen = avg_platen_disp[load_axis] / samp_L * 100.0
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


def set_strain_gauge(model, gauge_length=None, gauge_width=None, c_center=None):
    """
    Calculate local strain based on the dimensions of a virtual strain gauge placed at the center of teh model with
    x/y dimensions. By default, set to 0.25 of the length/width.

    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param gauge_length: length of the virtual strain gauge
    :type gauge_length: float
    :param gauge_width: width of the virtual strain gauge
    :type gauge_width: float
    :param c_center: User-defined center of the SG
    :type c_center: None or list[float, float, float]

    :return: Cells that cover the horizontal and vertical gauges as well as the gauge width and length
    :rtype: [list, list, float, float]
    """

    pv, ph = [], []

    if not gauge_width or gauge_width == 0:
        gauge_width = 0.25 * model.sample_width
    if not gauge_length or gauge_length == 0:
        gauge_length = 0.25 * model.sample_height

    if c_center:
        x_cor, y_cor, z_cor = c_center[0], c_center[1], c_center[2]
    else:
        x_cor, y_cor, z_cor = model.rock_model.center

    # Points that constitute the SG
    pv.append([x_cor + gauge_width / 2.0,
               y_cor - gauge_length / 2.0,
               0.0])
    pv.append([x_cor - gauge_width / 2.0,
               y_cor - gauge_length / 2.0,
               0.0])
    pv.append([x_cor + gauge_width / 2.0,
               y_cor + gauge_length / 2.0,
               0.0])
    pv.append([x_cor - gauge_width / 2.0,
               y_cor + gauge_length / 2.0,
               0.0])
    ph.append([x_cor + gauge_length / 2.0,
               y_cor + gauge_width / 2.0,
               0.0])
    ph.append([x_cor + gauge_length / 2.0,
               y_cor - gauge_width / 2.0,
               0.0])
    ph.append([x_cor - gauge_length / 2.0,
               y_cor + gauge_width / 2.0,
               0.0])
    ph.append([x_cor - gauge_length / 2.0,
               y_cor - gauge_width / 2.0,
               0.0])
    print('\tDimensions of SG are %s x %s' % (gauge_length, gauge_width))

    # Cells at the points of the SG
    # These dont change throughout the post-processing
    cv, ch = [], []
    for ps in range(0, len(pv)):
        cv.append(model.first_file.find_closest_cell(pv[ps]))
        ch.append(model.first_file.find_closest_cell(ph[ps]))

    # Verify that SG points are within the domain and return valid cells
    if -1 in cv or -1 in ch:
        print(formatting_codes.red_text("Check Strain Gauge Dimensions\nWill not process the strain gauges"))
        st_status = False
    else:
        print('\tVertical Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (pv, cv))
        print('\tHorizontal Gauges\n\t\textends between %s\n\t\tcover cells ID %s' % (ph, ch))

    # Global SG length
    global g_length
    g_length = gauge_length

    return ch, cv, gauge_width, gauge_length

def check_loading_direction(model, f1, f2):

    openfdem_model_ts_init = pv.read(f1)
    openfdem_model_ts_final = pv.read(f2)

    '''STRESS-STRAIN PLATENS'''

    platen = (openfdem_model_ts_init.threshold([model.platen_cells_elem_id, model.platen_cells_elem_id],
                                          model.var_data["mineral_type"]))
    top, bottom = (platen.get_data_range(model.var_data["boundary"]))

    disp_init = model.platen_info(openfdem_model_ts_init, top, model.var_data["platen_displacement"])
    disp_final = model.platen_info(openfdem_model_ts_final, top, model.var_data["platen_displacement"])

    avg_platen_init = [0.0, 0.0, 0.0] # Dummy cell
    avg_platen_final = [0.0, 0.0, 0.0] # Dummy cell

    for i in range(0, 3):
        # Convert forces from microN to kN and get the average forces
        avg_platen_init[i] = abs(disp_init[i])
        avg_platen_final[i] = abs(disp_final[i])

    diff_loading = [abs(avg_platen_final[i] - avg_platen_init[i]) for i in range(0, 3)]

    return diff_loading.index(max(diff_loading))

def main(model, platen_id, st_status, axis_of_loading, gauge_width, gauge_length, c_center, user_samp_A=None, user_samp_L=None, progress_bar=False):
    """
    Main concurrent Thread Pool to calculate the full stress-strain

    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param platen_id: Manual override of Platen ID
    :type platen_id: None or int
    :param st_status: Enable/Disable SG Calculations
    :type st_status: bool
    :param axis_of_loading: Enable/Disable SG
    :type axis_of_loading: None or int
    :param gauge_width: SG width
    :type gauge_width:  float
    :param gauge_length: SG length
    :type gauge_length: float
    :param c_center: User-defined center of the SG
    :type c_center: None or list[float, float, float]
    :param user_samp_A: Sample Area
    :type user_samp_A: None or float
    :param user_samp_L: Sample Length
    :type user_samp_L: None or float
    :param progress_bar: Show/Hide progress bar
    :type progress_bar: bool

    :return: full stress-strain information
    :rtype: pd.DataFrame
    """

    # Global declarations
    start = time.time()

    # Initialise Variables
    global history_strain, history_stress, gauge_disp_x, gauge_disp_y
    history_strain, history_stress = [], []
    gauge_disp_x, gauge_disp_y = [], []

    # File names of the basic files
    f_names = model._basic_files

    # Get rock dimension.
    x_dim, y_dim, z_dim, md_extent = model.rock_sample_dimensions(platen_id)
    dim_list = [x_dim, y_dim, z_dim]

    # Check UCS Simulation
    if model.simulation_type() != "UCS Simulation":
        print("Simulation appears to be not for compressive strength")
        exit("Simulation appears to be not for compressive strength")

    if axis_of_loading:
        print("\tPredefined user-defined loading axis [%s] is %s-direction" % (axis_of_loading, loading_dir_dict[axis_of_loading]))
        axis_of_loading = [axis_of_loading]
    elif model.number_of_points_per_cell != 3 and axis_of_loading is None:
        axis_of_loading = [check_loading_direction(model, f_names[0], f_names[-1])]
        print("\t3D Loading direction detected as [%s] is %s-direction" % (axis_of_loading[0], loading_dir_dict[axis_of_loading[0]]))
    else:
        print("\tPredefined loading Axis [1] is Y-direction")
        axis_of_loading = [1]

    # Determine A and L for Stress and Strain.
    global samp_A, samp_L
    # Sample length is the dimension in the direction of loading
    if user_samp_L is None:
        samp_L = dim_list[axis_of_loading[0]]
    else:
        samp_L = user_samp_L
    # del that dimension and keep the remaining ones.
    del dim_list[axis_of_loading[0]]
    temp_dim = dim_list

    # If 2D - the reduced list will be 0 and the width. Get width.
    if user_samp_A is None:
        if model.number_of_points_per_cell == 3:
            samp_A = max(temp_dim)
        # If 3D - locate [x, y, z] at the corner of the model.
        else:
            x_extent, y_extent, z_extent = md_extent[0:2], md_extent[2:4], md_extent[4:6]
            model_extent = [x_extent, y_extent, z_extent]
            f_cell_coord = []
            for idx in range(0, len(model_extent)):
                if idx == axis_of_loading[0]:
                    f_cell_coord.append(0)
                else:
                    f_cell_coord.append(model_extent[idx][0])
            # Check location. If -1 [IndexError] then circle, else square.
            try:
                model.find_cell(f_cell_coord)
                samp_A = math.prod(temp_dim)
            except IndexError:
                if temp_dim[0] == temp_dim[1]:
                    samp_A = 3.142 * temp_dim[0] * temp_dim[0] / 4
    else:
        samp_A = user_samp_A

    print("Values used in calculations are\n\tArea\t%0.2f\n\tLength\t%0.2f" % (samp_A, samp_L))

    # Initialise the Strain Gauges
    if model.number_of_points_per_cell != 3 and st_status:
        print(formatting_codes.red_text('Strain Gauges not supported in 3D\nWill not process the strain gauges'))
        cv, ch, gauge_width, gauge_length = [], [], 0, 0
        st_status = False
    elif st_status and model.number_of_points_per_cell == 3:  # Enabled SG st_status == True
        cv, ch, gauge_width, gauge_length = set_strain_gauge(model, gauge_width, gauge_length, c_center)
    else:
        cv, ch, gauge_width, gauge_length = [], [], 0, 0


    if st_status:  # Enabled SG st_status == True
        cv, ch, gauge_width, gauge_length = set_strain_gauge(model, gauge_width, gauge_length, c_center)
    else:
        cv, ch, gauge_width, gauge_length = [], [], 0, 0

    # Load basic files in the concurrent Thread Pool
    for fname in f_names:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(history_strain_func, fname, repeat(model), cv, ch, axis_of_loading))  # is self the list we are iterating over

    # Iterate through the files in the defined function
    for idx, fname_iter in enumerate(f_names):
        hist = history_strain_func(fname_iter, model, cv, ch, axis_of_loading)
        if progress_bar:
            formatting_codes.print_progress(idx + 1, len(f_names), prefix='Progress:', suffix='Complete')
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
