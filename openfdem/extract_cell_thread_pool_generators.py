import concurrent.futures
import time
from itertools import repeat

import pyvista as pv

import formatting_codes


def history_cellinfo_func(f_name, model, cell_id, array_needed):
    '''

    :param f_name: name of vtu file being processed
    :type f_name: str
    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param cell_id: ID of the cell from which the data needs to be extracted
    :type cell_id: int
    :param array_needed: Name of the property to extract
    :type array_needed: str
    :return: The value of the property from the cell being extracted
    :rtype: Generator[Tuple()]
    '''

    openfdem_model_ts = pv.read(f_name)

    # Extract Data and convert to list and get value
    cell_value = openfdem_model_ts.extract_cells([cell_id][0]).get_array(model.var_data[[array_needed][0]]).tolist()[0]

    # Append value to list
    cell_data.append(cell_value)
    yield cell_data


def main(model, cellid, arrayname):
    '''
    Main concurrent Thread Pool to get value of the property from the cell being extracted

    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param cellid: ID of the cell from which the data needs to be extracted
    :type cellid: int
    :param arrayname: Name of the property to extract
    :type arrayname: str
    :return: list of the values of the property from the cell being extracted
    :rtype: list
    '''

    global cell_data
    cell_data = []  # to reset the value everytime the function is called.

    # File names of the basic files
    f_names = model._basic_files

    # Global declarations
    start = time.time()

    # Load basic files in the concurrent Thread Pool
    for fname in f_names:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(history_cellinfo_func, f_names, repeat(model), [cellid], [arrayname]))  # is self the list we are iterating over

    # Iterate through the files in the defined function
    for fname_iter in f_names:
        hist = history_cellinfo_func(fname_iter, model, cellid, arrayname)
        hist.__next__()

    print(formatting_codes.calc_timer_values(time.time() - start))

    return cell_data
