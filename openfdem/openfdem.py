# Possible import statements to organize the code
# import _seismic_func
# import _fracture_func
# import _modelling_func
# import _data_func

import glob
import os.path as path
import re
import matplotlib.pyplot as plt
import pandas
import pandas as pd
import numpy as np
import math
import itertools
import pyvista as pv
import windrose
import random

# TODO:
#  3D Model - BD processing
#  Process 3D PLT Tests (Check that the platen ID is loaded correctly in all directions)
#  Process 3D PLT Tests (Identify platen ID automatically)
#  Process 3D UCS Tests (Check that the platen ID is loaded correctly in all directions)

# import complete_UCS_thread_pool_generators
# import complete_BD_thread_pool_generators
# from aggregate_storage import aggregate_storage


class Model:
    """Model class collects datafiles into one interface.
    
    Each data array returns as a list ordered by timestep
    Collection of timesteps?
    handles temporal manipulations
    
    :Example:
        >>> import openfdem as fdem
        >>> model = fdem.Model("../example_outputs/Irazu_UCS")
        
    """

    _file_names = {"OpenFDEM": {"basic": "_field_",
                                "broken_joint": "_brokenjoint_",
                                "soften_joint": "_softenjoint_",
                                "principal_stress_direction": "_principalstress_",
                                "acoustic_emission": "acousticemission"
                                },
                   "Irazu": {"basic": "_basic_",
                             "broken_joint": "_broken_joint_",
                             "soften_joint": "_softened_joint_",
                             "principal_stress_direction": "_principal_stress_strain_",
                             "acoustic_emission": "_seismic_event_"
                             }
                   }

    _var_dataset = {"openFDEM": {"mineral_type": "Property_id",
                                 "boundary": "boundary condition ID",
                                 "platen_force": "Force",
                                 "platen_displacement": "Displacement",
                                 "gauge_displacement": "Displacement",
                                 "temperature": 'temperature',
                                 },
                    "IRAZU": {"mineral_type": "material property ID",
                              "boundary": "boundary condition ID",
                              "platen_force": "force",
                              "platen_displacement": "displacement",
                              "gauge_displacement": "displacement",
                              "temperature": 'temperature',
                              },
                    }

    _plt_test_types = {"A": "Axial",
                      "B": "Blocky",
                      "D": "Diametral"}

    def _numericalSort(value):
        """Strip the numerical portion of the file.
        
        Sort filenames based on natural order (e.g. 1, 2,..., 10, 11, ..., instead of 1, 10, 11, 2, ...)
        
        :param value: Name of file
        :type value: str

        :return parts: Return the numerical portion of the file
        :rtype: list[str]
        
        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            >>> data._numericalSort('..\\example_outputs\\Irazu_UCS\\UCS_tutorial-Run_1_femdem.r2m_broken_joint_540000.vtu')
            ['..\\example_outputs\\Irazu_UCS\\UCS_tutorial-Run_',
             1,
             '_femdem.r',
             2,
             'm_broken_joint_',
             540000,
             '.vtu']
            
        """
        numbers = re.compile(r'(\d+)')
        parts = numbers.split(value)  # Split the numerical part of the file
        parts[1::2] = map(int, parts[1::2])  # Return the numerical portion of the file
        return parts

    def _findOutputFiles(self, dir_path, file_extension, output_type):
        """Find the collection of vtk files that belong to a specific group
        
        :param dir_path: Starting directory (full path is required)
        :type dir_path: str

        :return: A list of subdirectories
        :rtype: list[str]
        
        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            >>> data._findOutputFiles("../example_outputs/Irazu_UCS/",tuple(["vtu","vtp"]),"basic")[0:3]
            ['..\\\\example_outputs\\\\Irazu_UCS\\\\UCS_tutorial-Run_1_femdem.r2m_basic_0.vtu', '..\\\\example_outputs\\\\Irazu_UCS\\\\UCS_tutorial-Run_1_femdem.r2m_basic_20000.vtu', '..\\\\example_outputs\\\\Irazu_UCS\\\\UCS_tutorial-Run_1_femdem.r2m_basic_40000.vtu']
        """

        list_of_files = []
        # Get the list of files/directories
        for ext in file_extension:
            list_of_files.extend(
                glob.glob(dir_path + "/*" + Model._file_names[self._fdem_engine][output_type] + "*." + ext))

        # Sort list of files by their numerical value
        list_of_files = sorted(list_of_files, key=Model._numericalSort)

        # Ensure it is in proper system path format
        list_of_files = [path.relpath(vtkfile) for vtkfile in list_of_files]

        return list_of_files

    def __init__(self, folder=None, runfile=None, fdem_engine=None):
        """Create a Model object that finds all the run files and organizes all the file names on creation.

        :raise LookupError: Folder does not appear to exist or does not have valid output files.
        
        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
        """

        self._folder = folder
        self._runfile = runfile
        self._fdem_engine = fdem_engine
        if runfile is not None:
            if fdem_engine is None:
                if runfile.endswith(".y"):
                    self._fdem_engine = "OpenFDEM"
                elif runfile.endswith(".fdem"):
                    self._fdem_engine = "Irazu"
                else:
                    self._fdem_engine = "OpenFDEM"
        else:
            if fdem_engine is None:
                self._fdem_engine = "Irazu"
        if folder is not None:
            extensions = tuple(["vtu", "vtp"])
            self._basic_files = self._findOutputFiles(folder, extensions, "basic")
            if len(self._basic_files) == 0:
                raise LookupError('Folder does not appear to exist or does not have valid output files.')

            self._broken_files = self._findOutputFiles(folder, extensions, "broken_joint")
            self._soft_files = self._findOutputFiles(folder, extensions, "soften_joint")
            self._pstress_dir_files = self._findOutputFiles(folder, extensions, "principal_stress_direction")
            self._acoustic_files = self._findOutputFiles(folder, extensions, "acoustic_emission")

            self.first_file = pv.read(self._basic_files[0])
            self.n_timesteps = len(self._basic_files)
            self.n_points = self.first_file.number_of_points
            self.n_elements = self.first_file.n_cells
            self.number_of_points_per_cell = self.first_file.cell_n_points(0)
            # self.storage = aggregate_storage(folder,verbose=True)

        if self._fdem_engine == "Irazu":
            self.var_data = Model._var_dataset["IRAZU"]
        elif self._fdem_engine == "OpenFDEM":
            self.var_data = Model._var_dataset["openFDEM"]

    def __getitem__(self, key):
        """ Obtain all arrays by timestep access.

        :param key:Timestep accessible by index or string representing timestep
        :type key: Union[int, float]

        :return: Dataset of the timestep
        :rtype: Union[MultiBlock, UnstructuredGrid]

        :raise IndexError: Index outside model data range.

        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            >>> data['380000']
            UnstructuredGrid (0x7f1c7c897228)
                N Cells:	3968
                N Points:	11904
                X Bounds:	-3.217e+01, 3.149e+01
                Y Bounds:	-5.791e+01, 5.791e+01
                Z Bounds:	0.000e+00, 0.000e+00
                N Arrays:	15
            >>> data[0]
            UnstructuredGrid (0x7f1c7c897768)
                N Cells:	3968
                N Points:	11904
                X Bounds:	-2.800e+01, 2.800e+01
                Y Bounds:	-5.800e+01, 5.800e+01
                Z Bounds:	0.000e+00, 0.000e+00
                N Arrays:	15

        """

        # TODO:
        #  how to identify the data set being interrogated (basic / seismic?)

        try:
            if isinstance(key, int):
                self.timestep_data = pv.read(self._basic_files[key])
            if isinstance(key, str):
                indices = [i for i, elem in enumerate(self._basic_files) if "_" + key in elem]
                self.timestep_data = pv.read(self._basic_files[indices[0]])
            return self.timestep_data
        except IndexError:
            raise IndexError('Index outside model data range')

    # def __contains__(self, item):
    # """ Access to timestep in Model

    # def __len__(self):
    # """ Number of timestep files """

    def __scale_range__(self, scale_data):
        """
        Set the range of the scale min/max for the rosette diagram

        :param scale_data: Data
        :type scale_data: list[float]

        :return: Min/Max range of data and the data interval
        :rtype: list[float]
        """

        if (max(scale_data) - min(scale_data)) > 1:
            max_range = int(max(scale_data) + 1)
            min_range = int(min(scale_data))
            inter = 1
        else:
            max_range = (max(scale_data) + 0.1) // 0.1 * 0.1
            min_range = (min(scale_data)) // 0.1 * 0.1
            inter = 0.1
        return min_range, max_range, inter

    def model_dimensions(self, mat_id=None):
        """
        Function to get the "INITIAL" model bounds and returns the width, height, thickness

        :param mat_id: Optional, if a threshold is specific to a material type
        :type mat_id: int

        :return: model width, model height, model thickness
        :type: tuple[float, float, float]

        :Example:
            >>> import openfdem as fdem
            >>> model = fdem.Model("../example_outputs/Irazu_UCS")
            >>> # Returns the overall model dimensions
            >>> model.model_dimensions()
            (56.0, 116.0, 0.0)
            >>> # Returns the model dimensions based on material id 1
            >>> model.model_dimensions(1)
            (56.0, 116.0, 0.0)
            >>> # Error when material is not found
            >>> model.model_dimensions(3)
            IndexError: Material ID for platen out of range.
            Material Range 0-1
        """

        if mat_id is not None:
            self.threshold_bound_check(mat_id)
            self.thresholds_FDEM_output_files = self.first_file.threshold([mat_id, mat_id],
                                                                          self.var_data["mineral_type"])
        else:
            self.thresholds_FDEM_output_files = self.first_file

        sample_x_min, sample_x_max = self.thresholds_FDEM_output_files.bounds[0], self.thresholds_FDEM_output_files.bounds[1]
        sample_y_min, sample_y_max = self.thresholds_FDEM_output_files.bounds[2], self.thresholds_FDEM_output_files.bounds[3]
        sample_z_min, sample_z_max = self.thresholds_FDEM_output_files.bounds[4], self.thresholds_FDEM_output_files.bounds[5]
        self.model_width = sample_x_max - sample_x_min
        self.model_height = sample_y_max - sample_y_min
        self.model_thickness = sample_z_max - sample_z_min

        return self.model_width, self.model_height, self.model_thickness

    def model_domain(self):
        """
        Identifies the model domain by confirming the simulation cell vertex.
            2D (3 Points - Triangle)
            3D (4 Points - Tetrahedral)

        :return: number of nodes to skip in analysis
        :rtype: int

        :raise Warning: Simulation partially supported.

        :Example:
            >>> import openfdem as fdem
            >>> model = fdem.Model("../example_outputs/Irazu_UCS")
            >>> model.model_domain()
            2D Simulation
            4
        """

        if self.number_of_points_per_cell == 3:
            print("2D Simulation")
            node_skip = 4
        else:  # 3D (Tetrahedral)
            print("3D Simulation")
            node_skip = 6
            raise Warning("3D Simulation partially supported")

        return node_skip

    def threshold_bound_check(self, thres_id, thres_array='boundary'):
        """
        Checks the material ID is a valid choice.

        :param thres_id: ID of the item ot be threshold
        :type thres_id: int
        :param thres_array: Array name of the item ot be threshold
        :type thres_array: str

        :return: ID of the material
        :rtype: int

        :raise IndexError: ID out of range.

        :Example:
            >>> import openfdem as fdem
            >>> model = fdem.Model("../example_outputs/Irazu_UCS")
            >>> model.threshold_bound_check(0)
            0
            >>> model.threshold_bound_check(5)
            IndexError: Material ID for platen out of range.
            Material Range 0-1
        """

        min_thres_id, max_thres_id = self.first_file.get_data_range(self.var_data[thres_array])

        if thres_id not in range(min_thres_id, max_thres_id + 1):
            raise IndexError("Threshold ID out of range.\n%s Range %s-%s" % (thres_array, min_thres_id, max_thres_id))
        else:
            return thres_id

    def openfdem_att_check(self, att):
        """
        Checks that the attribute is a valid choice.

        :param: Attribute
        :type: str

        :return: Attribute
        :rtype: str

        :raise KeyError: Attribute does not exist.

        :Example:
            >>> import openfdem as fdem
            >>> model = fdem.Model("../example_outputs/Irazu_UCS")
            >>> model.openfdem_att_check('mineral_type')
            'mineral_type'
            >>> model.openfdem_att_check('material_property')
            KeyError: Attribute does not exist.
            Available options are mineral_type, boundary, platen_force, platen_displacement, gauge_displacement'
        """

        if att not in self.var_data.keys():
            raise KeyError(
                "Attribute does not exist.\nAvailable options are %s" % ", ".join(list(self.var_data.keys())))
        else:
            return att

    def rock_sample_dimensions(self, platen_id=None):
        """
        Lookup cell element ID on the top center and then trace points Using this information, we obtain the platen prop ID.
        Alternatively the user can define the material ID to exclude

        :param platen_id: Manual override of Platen ID
        :type platen_id: None or int

        :return: sample width, sample height, sample thickness
        :rtype: tuple[float, float, float]

        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            >>> # Let the script try to identify the platen material ID
            >>> data.rock_sample_dimensions()
            Script Identifying Platen
                Platen Material ID found as [1]
            (52.0, 108.0, 0.0)
            >>> # Explicitly defined the platen material ID
            >>> data.rock_sample_dimensions(0)
            User Defined Platen ID
                Platen Material ID found as [0]
            (56.0, 116.0, 0.0)
            >>> # Explicitly defined the platen material ID is out of range
            >>> data.rock_sample_dimensions(3)
            IndexError: Material ID for platen out of range.
            Material Range 0-1
        """

        top_center_point = [self.first_file.GetCenter()[0], self.first_file.bounds[3], self.first_file.bounds[5]]

        top_center_cell = self.first_file.extract_cells(self.first_file.find_closest_cell(top_center_point))

        if platen_id is None:
            print("Script Identifying Platen")
            if top_center_cell == -1:
                print("Unable to identify Platen ID Correctly.")
            self.platen_cells_elem_id = pv.cell_array(top_center_cell, self.var_data['mineral_type'])
        else:
            print("User Defined Platen ID")
            self.threshold_bound_check(platen_id)
            self.platen_cells_elem_id = platen_id

        print("\tPlaten Material ID found as %s" % self.platen_cells_elem_id)

        self.all_elem_id = list(set(self.first_file.get_array(self.var_data["mineral_type"])))
        self.rock_elem_ids = [x for x in self.all_elem_id if x != self.platen_cells_elem_id]
        self.rock_elem_ids_max, self.rock_elem_ids_min = max(self.rock_elem_ids), min(self.rock_elem_ids)

        self.rock_model = (self.first_file.threshold([self.rock_elem_ids_min, self.rock_elem_ids_max], self.var_data["mineral_type"]))
        self.platen = (self.first_file.threshold([self.platen_cells_elem_id, self.platen_cells_elem_id], self.var_data["mineral_type"]))

        sample_x_min, sample_x_max = self.rock_model.bounds[0], self.rock_model.bounds[1]
        sample_y_min, sample_y_max = self.rock_model.bounds[2], self.rock_model.bounds[3]
        sample_z_min, sample_z_max = self.rock_model.bounds[4], self.rock_model.bounds[5]
        self.sample_width = sample_x_max - sample_x_min
        self.sample_height = sample_y_max - sample_y_min
        self.sample_thickness = sample_z_max - sample_z_min
        self.rock_model_extents = [sample_x_min, sample_x_max, sample_y_min, sample_y_max, sample_z_min, sample_z_max]

        return self.sample_width, self.sample_height, self.sample_thickness, self.rock_model_extents

    def simulation_type(self):
        """ Identifies the type of simulation running. BD or UCS.
        Checks the top left corner of the model. If it contains material it is assumed as a rectangle.

        :return: Type of simulation. BD/UCS
        :rtype: str

        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            >>> data.rock_sample_dimensions()
            Script Identifying Platen
                Platen Material ID found as [1]
            (52.0, 108.0, 0.0)
            >>> data.simulation_type()
            'UCS Simulation'
        """
        self.check_edge_point = [self.rock_model.bounds[1], self.rock_model.bounds[3], self.rock_model.bounds[5]]
        self.check_edge_cell = self.first_file.extract_cells(self.first_file.find_closest_cell(self.check_edge_point))
        self.check_edge_cell = pv.cell_array(self.check_edge_cell, self.var_data['mineral_type'])

        if self.check_edge_cell == -1:
            self.sim_type = "BD Simulation"
        else:
            self.sim_type = "UCS Simulation"
        return self.sim_type

    def platen_info(self, pv_cells, platen_boundary_id, var_property):
        """
        This function thresholds cells based on boundary condition and sums them based on the defined parameter var_property

        :param pv_cells:
        :type pv_cells: pyvista.core.pointset.UnstructuredGrid or DataSet
        :param platen_boundary_id: boundary id that the threshold should be based on
        :type platen_boundary_id: float
        :param var_property: name of the property (array to b returned)
        :type var_property: str

        :return: array of the property based on the threshold
        :rtype: ndarray
        """

        platen_cell_ids = pv_cells.threshold([platen_boundary_id, platen_boundary_id], self.var_data["boundary"])
        platen_var_prop_list = sum(platen_cell_ids.get_array(var_property))

        if var_property == 'displacement':
            for k in range(0, 3):
                # divide by the number of points per cell  (3 in 2D and 4 in 3D)
                platen_var_prop_list[k] = platen_var_prop_list[k] / (
                        platen_cell_ids.cell_n_points(0) * platen_cell_ids.n_cells)

        return platen_var_prop_list

    # def get_broken(self, mode_id=None):

    def find_cell(self, model_point):
        """
        Identify the nearest cell in the model to the defined point

        :param model_point: x,y,z of a point in the model which
        :type model_point: list[float, float, float]

        :return: the cell nearest to the point
        :rtype: int

        :raise IndexError: Point outside model domain.

        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            >>> data.find_cell([0, 0, 0])
            2167
            >>> data.find_cell([100, 100, 0])
            IndexError: Point outside model domain.
            X=56.0, Y=116.0, Z=0.0

        """

        if self.first_file.find_closest_cell(model_point) == -1:
            raise IndexError("Point outside model domain.\nX=%s, Y=%s, Z=%s" % (
                self.model_dimensions()[0], self.model_dimensions()[1], self.model_dimensions()[2]))

        return self.first_file.find_closest_cell(model_point)

    def unpack_DataFrame(self, packed_cell_info):
        """
        Unpacking of the original array produced by pyvista
        If the array is a point data, the array is suffixed with _Nx where x is the node on that cell.

        :param packed_cell_info:
        :type packed_cell_info: pandas.DataFrame

        :return: Unpacked DataFrame
        :rtype: pandas.DataFrame
        """

        unpacked_DataFrame = pd.DataFrame()

        for idx, i in enumerate(list(packed_cell_info.columns)):
            name_list = []
            if len(packed_cell_info[i][0]) == 1:
                temp_list = packed_cell_info.explode(i)
                df = temp_list[i]
            else:
                for j in range(1, len(packed_cell_info[i][0]) + 1):
                    name_list.append(i + "_N%s" % j)
                df = pd.DataFrame([pd.Series(x) for x in packed_cell_info[i]])
                df.columns = name_list

            unpacked_DataFrame = pd.concat([unpacked_DataFrame, df], axis=1)

        return unpacked_DataFrame

    def extract_cell_info(self, cell_id, arrays_needed, progress_bar=True):
        """
        Returns the information of the cell based on the array requested.
        If the array is a point data, the array is suffixed with _Nx where x is the node on that cell.
        Also shows a quick example on how to plot the information extracted.

        :param cell_id: Cell ID to extract
        :type cell_id: int
        :param arrays_needed: list of array names to extract
        :type arrays_needed: list[str]
        :param progress_bar: Show/Hide progress bar
        :type progress_bar: bool

        :return: unpacked DataFrame
        :rtype: pandas.DataFrame

        :Example:
            >>> import openfdem as fdem
            >>> import matplotlib.pyplot as plt
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            >>> # Extract data platen_force', 'mineral_type' from Cell ID 1683
            >>> extraction_of_cellinfo = data.extract_cell_info(1683, ['platen_force', 'mineral_type'])
            Columns:
                Name: platen_force_N1, dtype=object, nullable: False
                Name: platen_force_N2, dtype=object, nullable: False
                Name: platen_force_N3, dtype=object, nullable: False
                Name: mineral_type, dtype=object, nullable: False
            >>> # For noded information => PLOTTING METHOD ONE
            >>> x, y = [], []
            >>> for i, row in extraction_of_cellinfo.iterrows():
            >>>     x.append(i)
            >>>     y.append(row['platen_force_N2'][0])
            >>> plt.plot(x, y, c='red', label='platen_force_N2_x')
            [<matplotlib.lines.Line2D object at 0x7f08fe98a310>]
            >>> plt.legend()
             <matplotlib.legend.Legend object at 0x7f08fe9854c0>
            >>> plt.show()
            # For noded information => PLOTTING METHOD TWO
            >>> lx = extraction_of_cellinfo['platen_force_N2'].to_list()
            >>> lx1 = list(zip(*lx))
            >>> plt.plot(lx1[0], label='platen_force_N2_x')
            [<matplotlib.lines.Line2D object at 0x7f08fe859b20>]
            >>> plt.plot(lx1[1], label='platen_force_N2_y')
            [<matplotlib.lines.Line2D object at 0x7f08fe859e50>]
            >>> plt.plot(lx1[2], label='platen_force_N2_z')
            [<matplotlib.lines.Line2D object at 0x7f08fe86a160>]
            >>> plt.legend()
            <matplotlib.legend.Legend object at 0x7f08fe86a340>
            >>> plt.show()
            # For non-nonded information
            >>> plt.plot(lx1[0], label='mineral_type')
            [<matplotlib.lines.Line2D object at 0x7f08fe7e39a0>]
            >>> plt.legend()
            <matplotlib.legend.Legend object at 0x7f08fe7e39d0>
            >>> plt.show()

        """

        if not type(arrays_needed) == list:
            self.openfdem_att_check(arrays_needed)
        else:
            for array_needed in arrays_needed:
                self.openfdem_att_check(array_needed)

        try:
            from . import extract_cell_thread_pool_generators
        except ImportError:
            import extract_cell_thread_pool_generators

        packed_df = extract_cell_thread_pool_generators.main(self, cell_id, arrays_needed, progress_bar)

        unpacked_df = self.unpack_DataFrame(packed_df)

        return unpacked_df

    def convert_to_xyz_array(self, node_df):
        """
        Convert extracted node information into summation based on X, Y and Z

        :param node_df: Extracted node information
        :type node_df: pandas.DataFrame

        :return: A DataFrame with summations along X, Y, Z axis. Column names are ["sum_X", "sum_Y", "sum_Z"]
        :rtype: pandas.DataFrame

        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_3D_UCS")
            >>> # # Extract all cells that meet the criteria and split to nodewise data for each time step.
            >>> # In this case "BOUNDARY CONDITION" is set to "1" for the threshold with the "FORCE" being extracted at each node.
            >>> df = data.extract_threshold_info(thres_id=1, thres_array='boundary', arrays_needed=['platen_force'])
            >>> # Sum the X,Y,Z of all nodes for each time step.
            >>> df_sum = data.convert_to_xyz_array(df)
            >>> print(df_sum)
                           sum_X         sum_Y         sum_Z
                0   0.000000e+00  0.000000e+00  0.000000e+00
                1  -1.224291e+05 -2.348118e+09  4.645789e+04
                2  -8.768720e+04 -4.663436e+09  7.953211e+03
                3  -5.580583e+04 -6.948494e+09 -1.039933e+04
                4  -1.602602e+05 -9.240063e+09  1.065935e+04
                5  -1.588623e+05 -1.152608e+10  4.616695e+04
                ...

        """

        if isinstance(node_df, pd.DataFrame):
            list_sum = np.array(node_df.to_numpy().tolist()).sum(axis=1).tolist()
            df_xyz = pandas.DataFrame(list_sum, columns=["sum_X", "sum_Y", "sum_Z"])
        elif type(node_df) is dict:
            raise AttributeError("Please pass a name of the key. Available keys are %s" % ", ".join(list(node_df.keys())))
        else:
            raise AttributeError("Unknown data type.")

        return df_xyz

    def extract_threshold_info(self, thres_id, thres_array, arrays_needed, progress_bar=True):
        """
        Returns the information of the cell based on the array requested.
        If the array is a point data, the array is suffixed with _Nx where x is cell ID.
        Also shows a quick example on how to plot the information extracted.

        :param thres_id: Threshold ID to extract
        :type thres_id: int
        :param thres_array: Array name of item to threshold.
        :type thres_array: str
        :param arrays_needed: list of array names to extract
        :type arrays_needed: list[str]
        :param progress_bar: Show/Hide progress bar
        :type progress_bar: bool

        :return: A DataFrame or a series of DataFrames nested in a dictionary with the key being the name of the array needed
        :rtype: pandas.DataFrame or dict[pandas.DataFrame]

        :Example:
            >>> import openfdem as fdem
            >>> import matplotlib.pyplot as plt
            >>> data = fdem.Model("../example_outputs/Irazu_3D_UCS")
            >>> # Extract all cells that meet the criteria and split to nodewise data for each time step.
            >>> # In this case "BOUNDARY CONDITION" is set to "1" for the threshold with the "FORCE" being extracted at each node.
            >>> df = data.extract_threshold_info(thres_id=1, thres_array='boundary', arrays_needed=['platen_force'])
            Progress: |//////////////////////////////////////////////////| 100.0% Complete
            54.74 seconds.
            >>> # Sum the X,Y,Z of all nodes for each time step.
            >>> df_sum = data.convert_to_xyz_array(df)

        """

        if not type(arrays_needed) == list:
            self.openfdem_att_check(arrays_needed)
        else:
            for array_needed in arrays_needed:
                self.openfdem_att_check(array_needed)

        try:
            from . import extract_threshold_thread_pool_generators
        except ImportError:
            import extract_threshold_thread_pool_generators

        packed_df = extract_threshold_thread_pool_generators.main(self, thres_id, thres_array, arrays_needed, progress_bar)

        unpacked_df = self.unpack_DataFrame(packed_df)

        if type(arrays_needed) == list and len(arrays_needed) > 1:
            df_dict = {}
            for array_needed in arrays_needed:
                df_dict[array_needed] = unpacked_df.filter(regex=array_needed)
            return df_dict
        else:
            return unpacked_df

    # def model_composition(self):
    # if not self._basic_0_loaded:
    # load_basic(0)

    # def broken_joints(self, mode=None,):

    # def seismic_events(self,time_filter = []):

    # def seismic_clustering(self):

    # def set_strain_gauge(self,point,axis):

    def complete_UCS_stress_strain(self, platen_id=None, st_status=False, axis_of_loading=None, gauge_width=0, gauge_length=0, c_center=None, samp_A=None, samp_L=None, progress_bar=True):
        """
        Calculate the full stress-strain curve

        :param platen_id: Manual override of Platen ID
        :type platen_id: None or int
        :param st_status: Enable/Disable SG
        :type st_status: bool
        :param axis_of_loading: Loading Direction
        :type axis_of_loading: None or int
        :param gauge_width: width of the virtual strain gauge
        :type gauge_width: float
        :param gauge_length: length of the virtual strain gauge
        :type gauge_length: float
        :param c_center: User-defined center of the SG
        :type c_center: None or list[float, float, float]
        :param samp_A: Sample Area
        :type samp_A: None or float
        :param samp_L: Sample Length
        :type samp_L: None or float
        :param progress_bar: Show/Hide progress bar
        :type progress_bar: bool

        :return: full stress-strain information
        :rtype: pandas.DataFrame

        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            # Minimal Arguments
            >>> df_wo_SG = data.complete_UCS_stress_strain()
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
            # full stress-strain without SG
            >>> df_wo_SG = data.complete_UCS_stress_strain(None, False)
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
            # full stress-strain with SG and default dimensions
            >>> df_Def_SG = data.complete_UCS_stress_strain(None, True)
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
                Name: Gauge Displacement X, dtype=float64, nullable: False
                Name: Gauge Displacement Y, dtype=float64, nullable: False
            # full stress-strain with SG and user-defined dimensions
            >>> df_userdf_SG = data.complete_UCS_stress_strain(None, True, 10, 10)
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
                Name: Gauge Displacement X, dtype=float64, nullable: False
                Name: Gauge Displacement Y, dtype=float64, nullable: False
        """

        try:
            from . import complete_UCS_thread_pool_generators
        except ImportError:
            import complete_UCS_thread_pool_generators

        return complete_UCS_thread_pool_generators.main(self, platen_id, st_status, axis_of_loading, gauge_width, gauge_length, c_center, samp_A, samp_L, progress_bar)

    def complete_BD_stress_strain(self, st_status=False, gauge_width=0, gauge_length=0, c_center=None, progress_bar=True):
        """
        Calculate the full stress-strain curve

        :param st_status: Enable/Disable SG
        :type st_status: bool
        :param gauge_width: width of the virtual strain gauge
        :type gauge_width: float
        :param gauge_length: length of the virtual strain gauge
        :type gauge_length: float
        :param c_center: User-defined center of the SG
        :type c_center: None or list[float, float, float]
        :param progress_bar: Show/Hide progress bar
        :type progress_bar: bool

        :return: full stress-strain information
        :rtype: pandas.DataFrame

        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/OpenFDEM_BD")
            # full stress-strain without SG
            >>> df_wo_SG = data.complete_BD_stress_strain(False)
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
            # full stress-strain with SG and default dimensions
            >>> df_Def_SG = data.complete_BD_stress_strain(True)
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
                Name: Gauge Displacement X, dtype=float64, nullable: False
                Name: Gauge Displacement Y, dtype=float64, nullable: False
            # full stress-strain with SG and user-defined dimensions
            >>> df_userdf_SG = data.complete_BD_stress_strain(True, 10, 10)
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
                Name: Gauge Displacement X, dtype=float64, nullable: False
                Name: Gauge Displacement Y, dtype=float64, nullable: False
        """

        try:
            from . import complete_BD_thread_pool_generators
        except ImportError:
            import complete_BD_thread_pool_generators

        return complete_BD_thread_pool_generators.main(self, st_status, gauge_width, gauge_length, c_center, progress_bar)

    def complete_PLT_stress_strain(self, load_config, platen_id=None, axis_of_loading=None, De_squared=None, progress_bar=True):
        """
        Calculate the full stress-strain curve
        :param load_config: type of PLT Test. "A" "D" "B"
        :type load_config: str
        :param platen_id: Manual override of Platen ID
        :type platen_id: None or int
        :param axis_of_loading: Loading Direction
        :type axis_of_loading: None or int
        :param De_squared: equivalent core diameter (i.e., the value of De_squared)
        :type De_squared: None or float
        :param progress_bar: Show/Hide progress bar
        :type progress_bar: bool

        :return: full stress-strain information
        :rtype: pandas.DataFrame

        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            # Minimal Arguments
            >>> df_wo_SG = data.complete_UCS_stress_strain()
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
            # full stress-strain without SG
            >>> df_wo_SG = data.complete_UCS_stress_strain(None, False)
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
            # full stress-strain with SG and default dimensions
            >>> df_Def_SG = data.complete_UCS_stress_strain(None, True)
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
                Name: Gauge Displacement X, dtype=float64, nullable: False
                Name: Gauge Displacement Y, dtype=float64, nullable: False
            # full stress-strain with SG and user-defined dimensions
            >>> df_userdf_SG = data.complete_UCS_stress_strain(None, True, 10, 10)
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
                Name: Gauge Displacement X, dtype=float64, nullable: False
                Name: Gauge Displacement Y, dtype=float64, nullable: False
        """

        try:
            from . import complete_PLT_thread_pool_generators
        except ImportError:
            import complete_PLT_thread_pool_generators

        if load_config not in self._plt_test_types.keys():
            raise IndexError("Unknown PLT Simulation. Supported simulations are %s which correspond to %s" % (", ".join(list(self._plt_test_types.keys())),
                                                                                                              ", ".join(list(self._plt_test_types.values())))
                             )

        return complete_PLT_thread_pool_generators.main(self, load_config, platen_id, axis_of_loading, De_squared, progress_bar)


    def plot_stress_strain(self, strain, stress, ax=None, **plt_kwargs):
        """
        Simple plot of the stress-strain curve

        :param strain: X-axis data [Strain]
        :type strain: pandas.DataFrame
        :param stress: Y-axis data [Stress]
        :type stress: pandas.DataFrame
        :param ax: Matplotlib Axis
        :type ax: matplotlib
        :param plt_kwargs: `~matplotlib.Modules` submodules
        :type plt_kwargs:

        :return: Matplotlib AxesSubplots
        :rtype: Matplotlib Axis

        :Example:
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/OpenFDEM_BD")
            # Minimal Arguments
            >>> df_wo_SG = data.complete_UCS_stress_strain()
            Columns:
                Name: Platen Stress, dtype=float64, nullable: False
                Name: Platen Strain, dtype=float64, nullable: False
            >>> data.plot_stress_strain(df_wo_SG['Platen Strain'], df_wo_SG['Platen Stress'], label='stress-strain', color='green')
            <AxesSubplot:xlabel='Strain (-)', ylabel='Axial Stress (MPa)'>
        """

        if ax is None:
            ax = plt.gca()
        ax.plot(strain, stress, **plt_kwargs)  # example plot here
        plt.xlabel('Strain (-)')
        plt.ylabel('Axial Stress (MPa)')

        return ax

    def Etan50_mod(self, ucs_data, linear_bestfit=True, loc_stress='Platen Stress', loc_strain='Platen Strain', plusminus_range=1):
        """
        Tangent Elastic modulus at 50%. Calculates +/- number of datapoint from the 50% Stress. Defaults to +/- 1 datapoint.

        :param ucs_data: DataFrame containing the stress-strain data
        :type ucs_data: pandas.DataFrame
        :param linear_bestfit: Calculate data based on range extents or linear best fit line.
        :type linear_bestfit: bool
        :param loc_stress: Column to obtain stress from. Defaults to Platen Stress
        :type loc_stress: str
        :param loc_strain: Column to obtain strain from. Defaults to Platen Strain
        :type loc_strain: str
        :param plusminus_range: Range over which to calculate the Elastic modulus
        :type plusminus_range: int

        :return: Tangent Elastic modulus at 50% as a slope and Y-Intercept. Y-Intercept = 0 if linear_bestfit is False
        :rtype: list[float]

        :Example:
            >>> data = pv.read("../example_outputs/Irazu_UCS")
            >>> df_1 = data.complete_UCS_stress_strain()
            >>> data.Etan50_mod(df_1)[0]
            51683.94337878284
            >>> data.Etan50_mod(df_1, linear_bestfit=False)[0]
            51639.21679789497
            >>> df_1 = data.complete_UCS_stress_strain(st_status=True)
            >>> data.Etan50_mod(df_1, loc_strain='Gauge Displacement Y', plusminus_range=1)[0]
            51216.33411269702
        """

        # Find the nearest match to the 50% max stress value.
        ucs_data = ucs_data.reset_index(drop=True)
        df_sorttomax = ucs_data.iloc[0: ucs_data[loc_stress].idxmax()]

        df_sort = df_sorttomax.loc[(df_sorttomax[loc_stress] - (max(df_sorttomax[loc_stress]) / 2)).abs().argsort()[:1]]

        # Return the index of the nearest match
        mid_stress_id = df_sort.index[0]

        if linear_bestfit:
            ucs_data_trim = ucs_data.iloc[mid_stress_id - plusminus_range: mid_stress_id + plusminus_range]
            Etan50, b = np.polyfit(ucs_data_trim[loc_strain], ucs_data_trim[loc_stress], deg=1)
        else:
            # Calculate delta stress/strain between at +/- 1 from the index of the 50% max stress value
            delta_stress = (ucs_data[loc_stress][mid_stress_id + plusminus_range] - ucs_data[loc_stress][
                mid_stress_id - plusminus_range])
            delta_strain = (ucs_data[loc_strain][mid_stress_id + plusminus_range] - ucs_data[loc_strain][
                mid_stress_id - plusminus_range])

            Etan50, b = delta_stress / delta_strain, 0

        Etan50, b = Etan50 * 100, b

        return Etan50, b

    def Esec_mod(self, ucs_data, upperrange, loc_stress='Platen Stress', loc_strain='Platen Strain'):
        """
        Secant Modulus between 0 and upperrange. The upperrange can be a % or a fraction.

        :param ucs_data: DataFrame containing the stress-strain data
        :type ucs_data: pandas.DataFrame
        :param upperrange: Range over which to calculate the Secant Modulus
        :type upperrange: float
        :param loc_stress: Column to obtain stress from. Defaults to Platen Stress
        :type loc_stress: str
        :param loc_strain: Column to obtain strain from. Defaults to Platen Strain
        :type loc_strain: str

        :return: Secant Elastic modulus between 0 and upperrange
        :rtype: float

        :Example:
            >>> data = pv.read("../example_outputs/Irazu_UCS")
            >>> df_1 = data.complete_UCS_stress_strain(st_status=True)
            >>> data.Esec_mod(df_1, 0.5)
            51751.010161057035
            >>> data.Esec_mod(df_1, 0.5, loc_strain='Gauge Displacement Y')
            51279.95421163901
        """

        # Convert percentage to decimal
        if upperrange > 1:
            upperrange = upperrange / 100

        # Find the nearest match to the (upperrange * max stress) value.
        ucs_data = ucs_data.reset_index(drop=True)
        df_sorttomax = ucs_data.iloc[0: ucs_data[loc_stress].idxmax()]
        df_sort = df_sorttomax.iloc[
            (df_sorttomax[loc_stress] - (max(df_sorttomax[loc_stress]) * upperrange)).abs().argsort()[:1]]
        # Return the index of the nearest match
        sec_stress_id = df_sort.index[0]

        # Calculate delta stress/strain between 0 and the defined value
        delta_stress = (df_sorttomax[loc_stress][sec_stress_id + 1])
        delta_strain = (df_sorttomax[loc_strain][sec_stress_id + 1])

        Esec = (delta_stress / delta_strain) * 100

        return Esec

    def Eavg_mod(self, ucs_data, upperrange, lowerrange, linear_bestfit=True, loc_stress='Platen Stress',
                 loc_strain='Platen Strain'):
        """
        Average Elastic modulus between two ranges

        :param ucs_data: DataFrame containing the stress-strain data
        :type ucs_data: pandas.DataFrame
        :param upperrange: Upper range to calculate the average
        :type upperrange: float
        :param lowerrange: Lower range to calculate the average
        :type lowerrange: float
        :param linear_bestfit: Calculate data based on range extents or linear best fit line.
        :type linear_bestfit: bool
        :param loc_stress: Column to obtain stress from. Defaults to Platen Stress
        :type loc_stress: str
        :param loc_strain: Column to obtain strain from. Defaults to Platen Strain
        :type loc_strain: str

        :return: Average Elastic modulus
        :rtype: list[float]

        :raise ZeroDivisionError: The range over which to calculate the Eavg is too small. Consider a larger range.

        :Example:
            >>> data = pv.read("../example_outputs/Irazu_UCS")
            >>> df_1 = data.complete_UCS_stress_strain(st_status=True)
            >>> data.Eavg_mod(df_1, 0.5, 0.6)[0]
            51485.33001517835
            >>> data.Eavg_mod(df_1, 0.5, 0.6, 'Gauge Displacement Y')[0]
            50976.62587224803
        """

        # Convert percentage to decimal
        if upperrange > 1:
            upperrange = upperrange / 100
        if lowerrange > 1:
            lowerrange = lowerrange / 100

        # Find the nearest match to the (upperrange * max stress) value.
        ucs_data = ucs_data.reset_index(drop=True)
        df_max = ucs_data[loc_stress].idxmax()
        df_sort = ucs_data.iloc[0:df_max]
        df_sort_upper = df_sort.iloc[
            (df_sort[loc_stress] - (max(df_sort[loc_stress]) * upperrange)).abs().argsort()[:1]]
        # Find the nearest match to the (lowerrange * max stress) value.
        df_sort_lower = df_sort.iloc[
            (df_sort[loc_stress] - (max(df_sort[loc_stress]) * lowerrange)).abs().argsort()[:1]]

        # Return the index of the nearest match
        upper_stress_id = df_sort_upper.index[0]
        lower_stress_id = df_sort_lower.index[0]

        # Error as the range is too small, which wields the same output time step in the script.
        if upper_stress_id == lower_stress_id:
            raise ZeroDivisionError("The range over which to calculate the Eavg is too small. Consider a larger range.")

        if linear_bestfit:
            ucs_data_trim = ucs_data.loc[upper_stress_id: lower_stress_id]
            Eavg, b = np.polyfit(ucs_data_trim[loc_strain], ucs_data_trim[loc_stress], deg=1)
        else:
            # Calculate delta stress/strain between the defined value
            delta_stress = (ucs_data[loc_stress][upper_stress_id] - ucs_data[loc_stress][lower_stress_id])
            delta_strain = (ucs_data[loc_strain][upper_stress_id] - ucs_data[loc_strain][lower_stress_id])

            Eavg, b = delta_stress / delta_strain, 0

        Eavg, b = Eavg * 100, b

        return Eavg, b

    def extract_based_coord(self, thres_model, coord_xyz, location, include_cells=False, adjacent_cells=False):
        """
        Extract the vtkdata set based on the defined coord location in the x=0 y=1 z=2 location.

        :param thres_model: threshold dataset of the material id of the rock
        :type thres_model: pyvista.core.pointset.UnstructuredGrid
        :param coord_xyz: x=0 y=1 z=2
        :type coord_xyz: int
        :param location: Xmin/Xmax/Ymin/Ymax/Zmin/Zmax
        :type location: float
        :param include_cells: If True, extract the cells that contain at least one of the extracted points. If False, extract the cells that contain exclusively points from the extracted points list.
        :type include_cells: bool
        :param adjacent_cells: Specifies if the cells shall be returned or not
        :type adjacent_cells: bool

        :return: Pointset of the data being filtered
        :rtype: pyvista.core.pointset.UnstructuredGrid
        """

        extracted_cells = thres_model.extract_points(thres_model.points[:, coord_xyz] == location, include_cells=False, adjacent_cells=False)

        return extracted_cells

    def direct_shear_calculation(self, platen_id, array, progress_bar=True):
        """

        :param platen_id: Material id of the platen
        :type platen_id: int
        :param array: the name of the array to be extracted
        :type array: str
        :param progress_bar: Show/Hide progress bar
        :type progress_bar: bool

        :return: DataFrame containing the absolute value of the array for each identified corner. Absolute sum of the extracted array split in Top/Bottom ane Left/Rigth sub-set into Top/Bottom.
        :rtype: pandas.DataFrame

        :Example:
        >>> import openfdem as fdem
        >>> data = fdem.Model("/external/2D_shear_4mm_profile_normal_load_test")
        >>> df = data.direct_shear_calculation(platen_id=1, array='platen_force', progress_bar=True)
        User Defined Platen ID
            Platen Material ID found as 1
        No. of points
            Left	158
            Left_Top	78
            Left_Bottom	80
            Right	158
            Right_Top	76
            Right_Bottom	82
            Top	35
            Bottom	38
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(df['Left_Top'], label='Left Top')
        [<matplotlib.lines.Line2D object at 0x7fe71f187320>]
        >>> plt.plot(df['Left_Bottom'], label='Left Bottom')
        [<matplotlib.lines.Line2D object at 0x7f65cf8975f8>]
        >>> plt.plot(df['Left'], label='Left')
        [<matplotlib.lines.Line2D object at 0x7fe71f187390>]
        >>> plt.legend()
        <matplotlib.legend.Legend object at 0x7fe71f187668>
        >>> plt.show()
        """

        try:
            from . import direct_shear_thread_pool_generators
        except ImportError:
            import direct_shear_thread_pool_generators

        return direct_shear_thread_pool_generators.main(self, platen_id, self.var_data[array], progress_bar)

    def model_vertices(self, t_step=0, thres_id=None, thres_array="mineral_type"):
        """
        Returns a list of the vertices in the form of Point1, Point 2

        :param t_step: Time step in model. Default 0
        :type t_step: int
        :param thres_id: ID of item to threshold. Default None.
        :type thres_id: int
        :param thres_array: Array name of item to threshold. Default "mineral_type".
        :type thres_array: str

        :return: list of the verticies in the model and/or the threshold of it.
        :rtype: list[tuples]

        :Example:
        >>> import openfdem as fdem
        >>> data = fdem.Model("../example_outputs/Irazu_UCS")
        >>> len(data.model_vertices(t_step=0, thres_id=0, thres_array='mineral_type'))
        11196
        >>> len(data.model_vertices(t_step=0, thres_id=1, thres_array='boundary'))
        354
        """

        # Load time step data
        vertices = []
        openfdem_model_ts = pv.read(self._basic_files[t_step])

        # Check is user defined threshold based on array and ID
        if thres_id is not None:
            self.threshold_bound_check(thres_id, thres_array)
            openfdem_model_ts = openfdem_model_ts.threshold([thres_id, thres_id], self.var_data[thres_array])

        # Lookup each cell and get the co-ordinates and vertices to obtain element connectivity.
        for elem_id in range(0, openfdem_model_ts.n_cells):
            vertex = openfdem_model_ts.cell_points(elem_id)
            tuples = [tuple(x) for x in vertex]
            elem_connectivity = itertools.combinations(tuples, 2)
            # Sort the connectivity of points small to big
            for jj in elem_connectivity:
                if jj[0] > jj[1]:
                    jj = (jj[1], jj[0])
                vertices.append(tuple(jj))

        return vertices

    def mesh_geometry(self, vertices):
        """
        Returns a unique set of vertices and calculates their length and orientation.

        :param vertices: list of vertices in the model at a given time step
        :type vertices: list[tuples]

        :return: DataFrame of the vertices length and orientation
        :rtype: pandas.DataFrame

        :Example:
        >>> import openfdem as fdem
        >>> data = fdem.Model("../example_outputs/Irazu_UCS")
        >>> vert = data.model_vertices(t_step=0, thres_id=1, thres_array='mineral_type')
        >>> data.mesh_geometry(vert)
               Length       Angle
        0    2.236068   63.434949
        1    2.000000    0.000000
        2    2.363608   59.436301
        3    2.000000    0.000000
        4    2.244731  117.123188
        ..        ...         ...
        409  2.116948    0.287685
        410  2.000000    0.000000
        411  2.000000    0.000000
        412  1.802781   45.829911
        413  2.227619  116.002627
        <BLANKLINE>
        [414 rows x 2 columns]
        """

        ver_angles, ver_linelen = [], []
        # Get unique vector set
        unique_vertices = set(vertices)

        # Unzip into two points (p1 and p2)
        for un_vertex in unique_vertices:
            p1, p2 = un_vertex[0], un_vertex[1]
            # calculate Orientation and Length
            cal_angle = self.__GetAngleOfLineBetweenTwoPoints__(p1, p2)
            cal_len = self.__GetDistanceBetweenTwoPoints__(p1, p2)
            # append to list
            ver_angles.append(cal_angle)
            ver_linelen.append(cal_len)

        # convert to DataFrame
        ver_data = pd.DataFrame(list(zip(ver_linelen, ver_angles)), columns=['Length', 'Angle'])

        return ver_data

    def draw_rose_diagram(self, t_step, rose_data=None, thres_id=None, thres_array="mineral_type", rose_range="Length"):
        """
        Draw a wind rose diagram based on the information passed.

        :param t_step: Time step in model. Default 0
        :type t_step: int
        :param rose_data: User can bypass requirement and pass a DataFrame with the data. Should be 2 columns with the Angle being the 2nd. Default None
        :type rose_data: DataFrame
        :param thres_id: ID of item to threshold. Default None.
        :type thres_id: int
        :param thres_array: Array name of item to threshold. Default "mineral_type".
        :type thres_array: str
        :param rose_range: Range to calculate the windrose bins. Default "length"
        :type rose_range: str

        :return: windrose figure
        :rtype: matplotlib.pyplot

        :Example:
        >>> import openfdem as fdem

        >>> data = fdem.Model("../example_outputs/Irazu_UCS")
        >>> data.draw_rose_diagram(t_step=0)
        <module 'matplotlib.pyplot' from '/usr/local/lib/python3.8/dist-packages/matplotlib/pyplot.py'>
        >>> data.draw_rose_diagram(t_step=0, rose_range='Length', thres_id=0, thres_array='boundary')
        <module 'matplotlib.pyplot' from '/usr/local/lib/python3.8/dist-packages/matplotlib/pyplot.py'>
        >>> # If you want to save the figure to a pyplot format.
        >>> figure_name = data.draw_rose_diagram(t_step=0, rose_range='Length', thres_id=0, thres_array='boundary')
        >>> figure_name.savefig('/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/example.pdf')
        """

        if rose_data is not None:
            rose_data = rose_data
        else:
            if thres_id is not None:
                self.threshold_bound_check(thres_id, thres_array)
                vert = self.model_vertices(t_step, thres_id, thres_array)
            else:
                vert = self.model_vertices(t_step)

            rose_data = self.mesh_geometry(vert)

        min_range, max_range, inter = self.__scale_range__(rose_data[rose_range])

        ''' Initialise Figure '''
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='windrose')
        viridis = plt.get_cmap('viridis')

        # Plot Sector -90 to 90
        # Set locations of X-Axis
        ax.set_thetagrids([90, 45, 0, -45, -90], labels=['0$^o$', '', '90$^o$', '', '0$^o$'])
        ax.set_thetamin(-90)
        ax.set_thetamax(90)
        ax.set_theta_zero_location("N")
        # Location of Radial Labels
        ax.set_rlabel_position(270)

        # Format the grid/plot
        plt.grid(alpha=0.5)
        plt.tight_layout()

        ''' Generate appropriate y-axis division'''
        # Create bins and counters to find mode value. Bins of 10 degrees.
        counts, bins = np.histogram(rose_data["Angle"], bins=19, range=(-5, 185))
        divisions = (int(str(1).ljust(len(str(max(counts))) - 1, '0')))
        count = (int(math.ceil((max(counts) + divisions) / divisions)) * divisions)
        grid_labels = np.linspace(0, count, 11, endpoint=True, dtype=int)

        ax.bar(rose_data["Angle"], rose_data[rose_range], opening=0.8, edgecolor='white', nsector=36, cmap=viridis,
               bins=np.arange(min_range, max_range, inter))

        # Set Legend And Legend Properties
        legend = ax.set_legend(title=rose_range, bbox_to_anchor=(0.5, 0.05), loc='center', ncol=4)
        legend.get_title().set_fontsize('10')
        legend.get_title().set_weight('bold')
        plt.setp(plt.gca().get_legend().get_texts(), fontsize='10')  # legend 'list' fontsize

        # Set locations of Y-Axis & Label
        ax.set_yticks(grid_labels)
        ax.tick_params(axis='y', pad=7.5, rotation=45)
        ax.set_yticklabels(grid_labels, verticalalignment="top", horizontalalignment='center',
                           size=8)  # works ONLY for on matplotlib 2.2.3+
        ax.text(0.75, 0.15, 'Frequency', horizontalalignment='center', verticalalignment='center', size=10,
                weight="bold", transform=ax.transAxes)

        # Figure Information
        plt.suptitle("Mesh Plot at Time Step %s\n%s lines" % (t_step, len(rose_data["Angle"])), fontsize=16, y=0.85, weight="bold")

        return plt

    def __GetDistanceBetweenTwoPoints__(self, p1, p2):
        """  Calculate distance between two points

        :param p1: First co-ordinate
        :type p1: list[float, float, float]
        :param p2: Second co-ordinate
        :type p2: list[float, float, float]
        :return: Edge length
        :rtype: float
        """

        x2, y2, z2 = p1
        x1, y1, z1 = p2
        edg_len = ((float(x2) - float(x1)) ** 2 +
                   (float(y2) - float(y1)) ** 2 +
                   (float(z2) - float(z1)) ** 2) ** 0.5
        return edg_len

    def __GetAngleOfLineBetweenTwoPoints__(self, p1, p2):
        """  Calculate angle between two points

        :param p1: First co-ordinate
        :type p1: list[float, float, float]
        :param p2: Second co-ordinate
        :type p2: list[float, float, float]
        :return: Angle 0-180
        :rtype: float
        """

        xDiff = p2[0] - p1[0]
        yDiff = p2[1] - p1[1]
        ang = math.degrees(math.atan2(yDiff, xDiff))
        if ang < 0:
            ang += 180
        return ang

    # def process_BD(self):

    # def magnitude_histogram(self):

    # def process_BD(self):

    # def magnitude_histogram(self):


class Timestep:
    """A class handling the data of each timestep.
    
    Each data array returns for only the timestep
    handles spatial manipulations
    """

    def __init__(self, file, runfile=None):
        _file = file
        _runfile = runfile

        __basic_file_loaded = False
        __softened_file_loaded = False
        __broken_file_loaded = False
        __principal_file_loaded = False
        __flow_channel_file_loaded = False
        __flow_direction_file_loaded = False

    # def _load_file(self):

    # def generate_rose_diagram(self):

    # def strain(self):

    # def stress(self,material_id):

    # def mat_id(self):

    # def prin_stress(self):

    # def displacement(self):

    # def fluid_press(self):

    # def force(self):

    # def velocity(self):

    # def active_cells(self):

    # def mass(self):

    # def mean_stress(self):

    # def prin_dev_stress(self):

    # def prin_strain(self):

    # def vol_strain(self):

    # def bound_id(self):
