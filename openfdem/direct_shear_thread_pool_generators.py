import concurrent.futures
import time
from itertools import repeat

import pandas as pd
import pyvista as pv

try:
    from . import formatting_codes
except ImportError:
    import formatting_codes


def abs_sum_array(f_name, model, array, edge_list):
    """
    Calculates the absolute sum of the array being interrogated in the sample.

    :param f_name: name of vtu file being processed
    :type f_name: str
    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param array: the name of the array to be extracted
    :type array: str
    :param edge_list: dictionary of location:integers that represent the nodes of the synthetic sample.
    :type edge_list: Dict[str, int]

    :return: list of the absolute sum of the array being interrogated.
    :rtype: list[float]

    :raise IndexError: Unknown Location
    """

    openfdem_model_ts = pv.read(f_name)

    for loc, point_ids_edge in edge_list.items():

        processed_file_extracted = openfdem_model_ts.extract_points(point_ids_edge, include_cells=False)

        if loc.startswith(("Top", "Bottom")):
            lx1 = list(zip(*processed_file_extracted.get_array(array)))
            lx2 = sum(abs(i) for i in lx1[1])
        elif loc.startswith(("Right", "Left")):
            lx1 = list(zip(*processed_file_extracted.get_array(array)))
            lx2 = sum(abs(i) for i in lx1[0])
        else:
            raise IndexError("Unknown Location")

        abs_sum_array_data.append(lx2)

    yield abs_sum_array_data


def sub_filter(vtk_data, y_middle):
    """
    Identify the points on the top/bottom half of the sample. Assumes symmetric sample and the middle Y is the center.

    :param vtk_data: Pointset of the data being filtered
    :type vtk_data: pyvista.core.pointset.UnstructuredGrid
    :param y_middle: mid-point of Y based on sample bounds
    :type y_middle: float

    :return: list of integers that represent the nodes on the top and bottom halves of the DS synthetic sample.
    :rtype: list[int]
    """

    extracted_top_list = []
    extracted_bottom_list = []

    for idx, y_coord in enumerate(vtk_data.points[:,1]):
        if y_coord > y_middle :
            extracted_top_list.append(vtk_data['vtkOriginalPointIds'][idx])
        if y_coord < y_middle :
            extracted_bottom_list.append(vtk_data['vtkOriginalPointIds'][idx])

    return extracted_top_list, extracted_bottom_list


def main(model, platen_id, array, progress_bar=True):
    """
    Main concurrent Thread Pool to calculate the absolute sum of data along edge of sample.

    :param model: FDEM Model Class
    :type model: openfdem.openfdem.Model
    :param platen_id: Manual override of Platen ID
    :type platen_id: int
    :param array: the name of the array to be extracted
    :type array: str
    :param progress_bar: Show/Hide progress bar
    :type progress_bar: bool

    :return: Absolute sum of the extracted array split in Top/Bottom ane Left/Rigth sub-set into Top/Bottom.
    :rtype: pd.DataFrame
    """

    # Initialise Variables
    global abs_sum_array_data, ts_data
    # To reset the value everytime the function is called.
    abs_sum_array_data = []  # Stores "array" in each timestep
    ts_data = {}  # Compiles all data into a dictionary

    # File names of the basic files
    f_names = model._basic_files

    # Threshold data by material ID
    # Get the bounds of the threshold dataset (min,max X - Y - Z)
    _, _, _, rock_model_bounds = model.rock_sample_dimensions(platen_id)
    sample_x_min, sample_x_max, sample_y_min, sample_y_max, sample_z_min, sample_z_max = rock_model_bounds

    # Mid Y to get Top/Bottom
    y_middle = (sample_y_min + sample_y_max) / 2

    extracted_left = model.extract_based_coord(model.rock_model, 0, sample_x_min)
    # Filter top/bottom
    left_top, left_bottom = sub_filter(extracted_left, y_middle)
    # Make vtk subset based on filter
    extracted_left_top = model.rock_model.extract_points(left_top, include_cells=False, adjacent_cells=False)
    extracted_left_bottom = model.rock_model.extract_points(left_bottom, include_cells=False, adjacent_cells=False)

    extracted_right = model.extract_based_coord(model.rock_model, 0, sample_x_max)
    # Filter top/bottom
    right_top, right_bottom = sub_filter(extracted_right, y_middle)

    # Make vtk subset based on filter
    extracted_right_top = model.rock_model.extract_points(right_top, include_cells=False, adjacent_cells=False)
    extracted_right_bottom = model.rock_model.extract_points(right_bottom, include_cells=False, adjacent_cells=False)

    extracted_top = model.extract_based_coord(model.rock_model, 1, sample_y_max)
    extracted_bottom = model.extract_based_coord(model.rock_model, 1, sample_y_min)

    # Compile data into a dictionary
    ids_DS = {"Left": extracted_left,
              "Left_Top": extracted_left_top,
              "Left_Bottom": extracted_left_bottom,
              "Right": extracted_right,
              "Right_Top": extracted_right_top,
              "Right_Bottom": extracted_right_bottom,
              "Top": extracted_top,
              "Bottom": extracted_bottom }

    # No. of extracted points
    print("No. of points")
    for loc, vtk_ids in ids_DS.items():
        print("\t%s\t%s" % (loc, vtk_ids.n_points))  # No. of extracted points
        # print(vtk_ids['vtkOriginalPointIds'])  # Ids of extracted points

    ids_of_edge_points = {}

    # Ids of extracted points in each location
    for loc, vtk_ids in ids_DS.items():
        ids_of_edge_points[loc] = vtk_ids['vtkOriginalPointIds']

    # Global declarations
    start = time.time()

    # Load basic files in the concurrent Thread Pool
    for fname in f_names:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(abs_sum_array, fname, repeat(model), array, ids_of_edge_points))  # is self the list we are iterating over

    # Iterate through the files in the defined function
    for idx, fname_iter in enumerate(f_names):
        hist = abs_sum_array(fname_iter, model, array, ids_of_edge_points)
        if progress_bar:
            formatting_codes.print_progress(idx + 1, len(f_names), prefix='Progress:', suffix='Complete')
        hist.__next__()

    # Convert concurrent Thread Pool from list to Dictionary
    no_of_sides = len(ids_DS.keys())
    # Split list and make into a dictionary
    ts_counter = 0
    for idx in range(0, len(abs_sum_array_data), no_of_sides):
        ts_data[ts_counter] = abs_sum_array_data[idx:idx+no_of_sides]
        ts_counter += 1

    # Convert the dictionary into a DataFrame
    ts_DS_df = pd.DataFrame.from_dict(ts_data, orient='index', columns=ids_DS.keys())
    Summation_Top_Bottom = ts_DS_df["Top"] + ts_DS_df["Bottom"]
    ts_DS_df['Summation_Top_Bottom'] = Summation_Top_Bottom

    print(formatting_codes.calc_timer_values(time.time() - start))

    return ts_DS_df

