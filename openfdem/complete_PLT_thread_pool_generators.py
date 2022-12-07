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

def history_strain_func(f_name, model, axis, load_config):
    """
    Calculate the Uncorrected Point Load Strength Index from platens

    :param f_name: name of vtu file being processed
    :type f_name: str
    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param axis: loading axis
    :type axis: [int]
    :param load_config: type of PLT Test. "A" "D" "B"
    :type load_config: [str]

    :return: Uncorrected Point Load Strength Index
    :rtype: Generator[Tuple[list], Any, None]
    """

    openfdem_model_ts = pv.read(f_name)

    '''STRESS-STRAIN PLATENS'''

    platen = (openfdem_model_ts.threshold([model.platen_cells_elem_id, model.platen_cells_elem_id],
                                          model.var_data['basic']["mineral_type"]))
    top, bottom = (platen.get_data_range(model.var_data['basic']["boundary"]))

    top_platen_force_list = model.platen_info(openfdem_model_ts, top, model.var_data['basic']["platen_force"])
    bot_platen_force_list = model.platen_info(openfdem_model_ts, bottom, model.var_data['basic']["platen_force"])

    # avg_top_platen_disp = model.platen_info(openfdem_model_ts, top, model.var_data["platen_displacement"])
    # avg_bottom_platen_disp = model.platen_info(openfdem_model_ts, bottom,
    #                                            model.var_data["platen_displacement"])

    # avg_platen_disp = [0.0, 0.0, 0.0] # Dummy cell
    avg_platen_force = [0.0, 0.0, 0.0] # Dummy cell
    load_axis = axis[0]  # Axis of loading in Y direction.

    for i in range(0, 3):
        # Convert forces from microN to kN and get the average forces
        avg_platen_force[i] = 0.5 * (abs(top_platen_force_list[i]) + abs(bot_platen_force_list[i])) / 1.0e9
        # avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

    # Calculate the stress in MPa (force in kN & area in mm^2)
    stress_from_platen = avg_platen_force[load_axis] / (samp_De_square) * 1.0e3
    history_stress.append(stress_from_platen)

    # # Calculate strains in percentage (%)
    # strain_from_platen = avg_platen_disp[load_axis] / samp_L * 100.0
    # history_strain.append(strain_from_platen)

    yield history_stress


def check_loading_direction(model, f1, f2):

    openfdem_model_ts_init = pv.read(f1)
    openfdem_model_ts_final = pv.read(f2)

    '''STRESS-STRAIN PLATENS'''

    platen = (openfdem_model_ts_init.threshold([model.platen_cells_elem_id, model.platen_cells_elem_id],
                                          model.var_data['basic']["mineral_type"]))
    top, bottom = (platen.get_data_range(model.var_data['basic']["boundary"]))

    disp_init = model.platen_info(openfdem_model_ts_init, top, model.var_data['basic']["platen_displacement"])
    disp_final = model.platen_info(openfdem_model_ts_final, top, model.var_data['basic']["platen_displacement"])

    avg_platen_init = [0.0, 0.0, 0.0] # Dummy cell
    avg_platen_final = [0.0, 0.0, 0.0] # Dummy cell

    for i in range(0, 3):
        # Convert forces from microN to kN and get the average forces
        avg_platen_init[i] = abs(disp_init[i])
        avg_platen_final[i] = abs(disp_final[i])

    diff_loading = [abs(avg_platen_final[i] - avg_platen_init[i]) for i in range(0, 3)]

    return diff_loading.index(max(diff_loading))

def main(model, load_config, platen_id, axis_of_loading, De, progress_bar=False):
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

    # # Check UCS Simulation
    # if model.simulation_type() != "UCS Simulation":
    #     print("Simulation appears to be not for compressive strength")
    #     exit("Simulation appears to be not for compressive strength")

    if axis_of_loading:
        print("\tPredefined user-defined loading axis [%s] is %s-direction" % (axis_of_loading, loading_dir_dict[axis_of_loading]))
        axis_of_loading = [axis_of_loading]
    elif model.number_of_points_per_cell != 3 and axis_of_loading is None:
        axis_of_loading = [check_loading_direction(model, f_names[0], f_names[-1])]
        print("\t3D Loading direction detected as [%s] is %s-direction" % (axis_of_loading[0], loading_dir_dict[axis_of_loading[0]]))
    else:
        print("\tPredefined loading Axis [1] is Y-direction")
        axis_of_loading = [1]

    # Determine De for Stress
    global samp_De_square
    print("Loading Configuration %s" % model._plt_test_types[load_config])
    if De is None:
        samp_L = dim_list[axis_of_loading[0]]
        if load_config == "D":
            samp_De_square = dim_list[axis_of_loading[0]] * dim_list[axis_of_loading[0]]
        elif load_config == "A":
            temp_dim = [x for x in dim_list if x != samp_L]
            samp_De_square = 4 * (math.prod(temp_dim)) / 3.142
        elif load_config == "D":
            print("WD are taken as the max. dimension! Please be cautious")
            temp_dim = [x for x in dim_list if x != samp_L]
            samp_De_square = 4 * (math.prod(temp_dim)) / 3.142
    else:
        print("User defined De_squared")
        samp_De_square = De

    print("De_squared used in calculation is %0.2f" % samp_De_square)

    # samp_L = dim_list[axis_of_loading[0]]
    # temp_dim = [x for x in dim_list if x != samp_L]
    # if model.number_of_points_per_cell == 3:
    #     samp_W = max(temp_dim)
    # else:
    #     x_extent, y_extent, z_extent = md_extent[0:2], md_extent[2:4], md_extent[4:6]
    #     model_extent = [x_extent, y_extent, z_extent]
    #     f_cell_coord = []
    #     for idx in range(0, len(model_extent)):
    #         if idx == axis_of_loading[0]:
    #             f_cell_coord.append(0)
    #         else:
    #             f_cell_coord.append(model_extent[idx][0])
    #     try:
    #         model.find_cell(f_cell_coord)
    #         samp_W = math.prod(temp_dim)
    #     except IndexError:
    #         if temp_dim[0] == temp_dim[1]:
    #             samp_W = 3.142 * temp_dim[0] * temp_dim[0] / 4
    #
    # print(samp_W, samp_L)


    # Load basic files in the concurrent Thread Pool
    for fname in f_names:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(history_strain_func, fname, repeat(model), axis_of_loading, [load_config]))  # is self the list we are iterating over

    # Iterate through the files in the defined function
    for idx, fname_iter in enumerate(f_names):
        hist = history_strain_func(fname_iter, model, axis_of_loading, [load_config])
        if progress_bar:
            formatting_codes.print_progress(idx + 1, len(f_names), prefix='Progress:', suffix='Complete')
        hist.__next__()

    print(formatting_codes.calc_timer_values(time.time() - start))


    plt_df = pd.DataFrame(list(history_stress),
                          columns=['Platen Stress'])

    return plt_df
