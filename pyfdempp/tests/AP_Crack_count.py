# /////////////////////////////////////////////////////////////// #
# Python Script initially created on 2022-12-27
# Compiled by Aly @ Grasselli Geomechanics Group, UofT, 2022
# Created using PyCharm
# /////////////////////////////////////////////////////////////// #

import pyfdempp as fdem
print(fdem.__version__)
import matplotlib.pyplot as plt

# ll = {'model_1': 2.0, 'model_2': 1.2, 'model_7': 4.33, 'model_10': 4.0, 'model_12': 0.25, 'model_14': 0.23,
#       'model_15': 1.0, 'model_16': 0.94, 'model_17': 0.88, 'model_18': 0.83, 'model_19': 0.80}
#
# import operator
# import collections
#
# ll = sorted(ll.items(), key=operator.itemgetter(1))
# ll = collections.OrderedDict(ll)

import pandas as pd
#
# dict_Kdata = pd.DataFrame()
# dict_Kdata_UD = pd.DataFrame()
#
#
# def plot_stacked_bar(dict_to_plot):
#     fig, ax = plt.subplots()
#     dict_Kdata_T = dict_to_plot.T
#     dict_Kdata_T.plot.bar(stacked=True, )
#     plt.xticks(rotation=45)
#     plt.xlabel('Value of Anisotropy')
#     plt.ylabel('No. of Cracks')
#     plt.tight_layout()
#     plt.show()


CC_LUT = [0.0,
          1.0,
          1.0,
          1.5,
          1.5,
          2.0,
          2.0,
          4.0,
          6.0,
          6.0,
          6.5,
          6.5,
          7.0,
          7.0,
          9.0]
CC_LUT_NAME = ['edge',
               'pure tensile',
               'tensile dominant',
               '50-50',
               'shear dominant',
               'pure shear',
               'mixed mode',
               'broken DFN',
               'DFN pure tensile',
               'DFN tensile dominant',
               'DFN 50-50',
               'DFN shear dominant',
               'DFN pure shear',
               'DFN mixed mode']

# for k, v in ll.items():
#     data = fdem.Model("/external/OUTPUTS/AP_HF_Models/DFN_Lower_90/%s" % k)
#     model_cracks_SH_TEN_UD = data.crack_failure_mode_clustering(crack_LUT=CC_LUT, crack_LUT_name=CC_LUT_NAME)
#     model_cracks_SH_TEN = data.crack_failure_mode_clustering()
#     K_data = model_cracks_SH_TEN.iloc[-1]
#     K_data_UD = model_cracks_SH_TEN_UD.iloc[-1]
#
#     series_name = "K = %s" % str(v)
#     dict_Kdata[series_name] = K_data
#     dict_Kdata_UD[series_name] = K_data_UD
#
# plot_stacked_bar(dict_Kdata)
# plot_stacked_bar(dict_Kdata_UD)

# data3D = fdem.Model("/external/OUTPUTS/cube_FDEM_wPlaten_m3")
# df_3D = data3D.complete_UCS_stress_strain(axis_of_loading=2, platen_id=1)


# import pyfdempp as fdem
data = fdem.Model("/external/OUTPUTS/cube_FDEM_wPlaten_m3")
list_of_cellids = []
for i in range(0, 39):
    print(i)
    for j in range(0, 39):
        list_of_cellids.append(data.find_cell([i, j, 39.9]))
print("Cells on plane %s" % len(list_of_cellids))
list_of_cellids = list(set(list_of_cellids))
print("Unique cells on plane %s" % len(list_of_cellids))
sliced_data = data.slice_dataset(list_of_cellids, 'force')
cell_df = pd.DataFrame.from_dict(sliced_data)
xyz = data.convert_to_xyz_array(cell_df['force'])
xyz.to_csv('/home/aly/Desktop/SampleFace.csv')


df_id2 = data.extract_threshold_info(thres_id=2, thres_array='boundary', arrays_needed=['platen_force'], progress_bar=True)
df_sum = data.convert_to_xyz_array(df_id2)
df_id3= data.extract_threshold_info(thres_id=3, thres_array='boundary', arrays_needed=['platen_force'], progress_bar=True)
df_sum3 = data.convert_to_xyz_array(df_id3)
df_id2.to_csv('/home/aly/Desktop/df_id2.csv')
df_id3.to_csv('/home/aly/Desktop/df_id3.csv')

# print(df_sum)
# print(df_sum3)
# list_of_cellids = []
# for i in range (0, 39):
#     print(i)
#     for j in range (0, 39):
#         # print(i, j, 39.9)
#         # data.find_cell([i, j, 39.9])
#         list_of_cellids.append(data.find_cell([i, j, 39.9]))
