import concurrent.futures
import time
from itertools import repeat

import pandas as pd
import pyvista as pv

try:
    from . import formatting_codes
except ImportError:
    import formatting_codes

def history_modelthreshold_func(f_name, model, thres_id, thres_array, array_needed):
    """
    Generate a dictionary of the various array being interrogated for the said cell threshold criteria

    :param f_name: name of vtu file being processed
    :type f_name: str
    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param thres_id: ID of the threshold from which the data needs to be extracted
    :type thres_id: int
    :param thres_array: Array name of item to threshold.
    :type thres_array: str
    :param array_needed: Name of the property to extract
    :type array_needed: list[str]

    :return: The value of the property from the cell being extracted
    :rtype: Generator[Tuple()]
    """

    openfdem_model_ts = pv.read(f_name)
    ts_values = []
    # Extract Data and convert to list and get value

    for i_array_needed in array_needed:
        ts_values = openfdem_model_ts.threshold([thres_id, thres_id], model.var_data[thres_array]).get_array(model.var_data[i_array_needed]).tolist()
        dict_array[i_array_needed].append(ts_values)

    yield dict_array


def main(model, thresid, thresarray, arrayname, progress_bar=False):
    """
    Main concurrent Thread Pool to get value of the property from the cell being extracted

    :param model: FDEM Model Class
    :type model:  openfdem.openfdem.Model
    :param thresid: ID of the cell from which the data needs to be extracted
    :type thresid: int
    :param thresarray: Array name of item to threshold. Default "mineral_type".
    :type thresarray: str
    :param arrayname: Name of the property to extract
    :type arrayname: list[str]
    :param progress_bar: Show/Hide progress bar
    :type progress_bar: bool

    :return: DataFrame of the values of the property from the cell being extracted
    :rtype: pandas.DataFrame
    """

    global cell_data, dict_array
    # To reset the value everytime the function is called.
    cell_data = []
    dict_array = {}

    if type(arrayname) == list:
        for i in arrayname:
            dict_array[i] = []
    else:
        dict_array[arrayname] = []
        arrayname = [arrayname]

    # File names of the basic files
    f_names = model._basic_files

    # Global declarations
    start = time.time()

    # Load basic files in the concurrent Thread Pool
    for fname in f_names:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(history_modelthreshold_func, fname, repeat(model), [thresid], thresarray, arrayname))  # is self the list we are iterating over

    # Iterate through the files in the defined function
    for idx, fname_iter in enumerate(f_names):
        hist = history_modelthreshold_func(fname_iter, model, thresid, thresarray, arrayname)
        if progress_bar:
            formatting_codes.print_progress(idx + 1, len(f_names), prefix='Progress:', suffix='Complete')
        hist.__next__()

    print(formatting_codes.calc_timer_values(time.time() - start))

    # Convert the dictionary into a DataFrame
    cell_df = pd.DataFrame.from_dict(dict_array)

    return cell_df
