# Possible import statements to organize the code
# import _seismic_func
# import _fracture_func
# import _modelling_func
# import _data_func

import glob
import os
import os.path as path
import re

import pyvista as pv
import time
import concurrent.futures
from multiprocessing import Process
from threading import Thread


class Model:
    """Model class collects datafiles into one interface.
    
    Each data array returns as a list ordered by timestep
    Collection of timesteps?
    handles temporal manipulations
    
    Example:
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
                                "platen_force": "force",
                                "platen_displacement": "displacement",
                                },
                   "IRAZU": {"mineral_type": "material property ID",
                             "boundary": "boundary condition ID",
                             "platen_force": "force",
                             "platen_displacement": "displacement",
                             },
                   }
    def _numericalSort(value):
        """
        Strip the numerical portion of the file.
        
        Sort filenames based on natural order (e.g. 1, 2,..., 10, 11, ..., instead of 1, 10, 11, 2, ...)
        
        :param value: Name of file
        :type value: str
        :return parts: Return the numerical portion of the file
        :rtype: list[str]
        
        Example:
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
        
        Example:
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
        
        Example:: 
            >>> import openfdem as fdem
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            <BLANKLINE>
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

        Example:
            >>> import openfdem as fdem
            >>> timestep = fdem.model['380000']
            <BLANKLINE>
            >>> timestep = fdem.model[0]
            <BLANKLINE>
        """

        # TODO: how to identify the data set being interrogated (basic / seismic?)
        try :
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
    # """
    # def __len__(self):
    # """ Number of timestep files """

    def model_dimensions(self, mat_id=None):
        '''
        Function to get the "INITIAL" model bounds and returns the width, height, thickness

        :param mat_id: Optional, if a threshold is specific to a material type
        :type mat_id: int
        :return: model width, model height, model thickness
        :type: tuple[float, float, float]

        Example:
            >>> data = pv.read("../example_outputs/Irazu_UCS")
            # Returns the overall model dimensions
            >>> data.model_dimensions
            (56.0, 116.0, 0.0)
            # Returns the model dimensions based on material id 1
            >>> data.model_dimensions(1)
            (52.0, 108.0, 0.0)
            # Error when material is not found
            >>> data.model_dimensions(3)
            Material ID specified out of range
            Process finished with exit code 1
        '''

        if mat_id != None:
            min, max = self.first_file.get_data_range(self.var_data["mineral_type"])
            if mat_id in range(min, max + 1):
                self.thresholds_FDEM_output_files = self.first_file.threshold([mat_id, mat_id], self.var_data["mineral_type"])
            else:
                raise IndexError("Material ID specified out of range")
        else:
            self.thresholds_FDEM_output_files = self.first_file

        sample_x_min, sample_x_max  = self.thresholds_FDEM_output_files.bounds[0], self.thresholds_FDEM_output_files.bounds[1]
        sample_y_min, sample_y_max  = self.thresholds_FDEM_output_files.bounds[2], self.thresholds_FDEM_output_files.bounds[3]
        sample_z_min, sample_z_max = self.thresholds_FDEM_output_files.bounds[4], self.thresholds_FDEM_output_files.bounds[5]
        self.model_width = sample_x_max - sample_x_min
        self.model_height = sample_y_max - sample_y_min
        self.model_thickness = sample_z_max - sample_z_min

        return self.model_width, self.model_height, self.model_thickness


    def model_domain(self):
        '''
        Identifies the model domain if 2D (3 Points - Triangle) or 3D (4 Points - Tetrahedral) by confirming the simulation cell vertex.

        :return: number of nodes to skip in analysis
        :rtype: int
        '''

        if self.number_of_points_per_cell == 3:
            print(("2D Simulation"))
            node_skip = 4
        else:  # 3D (Tetrahedral)
            print(("3D Simulation"))
            node_skip = 6
            exit("3D Simulation not supported")

        return node_skip

    # def load_basic(self,timestep=None):

    def filter_material(self, material_id):
        ''' Filter based on the material id

        :param material_id: material id to be threshold
        :type material_id: int
        :return: Dataset of the threshold
        :rtype: Union[MultiBlock, UnstructuredGrid]
        '''
        #TODO: look into making this with several input MIN MAX range of which ARRAY?
        # Seems like it IS NOT for a MultiBlock condition?!
        ax = self._data.threshold([material_id, material_id], self.var_data["mineral_type"])
        return ax


    def rock_sample_dimensions(self):
        '''
        Lookup cell element ID on the top center and then trace points Using this information, we obtain the platen prop ID

        :return: sample width, sample height, sample thickness
        :type: tuple[float, float, float]
        '''

        top_center_point = [self.first_file.GetCenter()[0], self.first_file.bounds[3], self.first_file.bounds[5]]

        top_center_cell = self.first_file.extract_cells(self.first_file.find_closest_cell(top_center_point))

        if top_center_cell == -1:
            print("Unable to identify Platen ID Correctly.")
        self.platen_cells_elem_id = pv.cell_array(top_center_cell, self.var_data['mineral_type'])

        print("\tPlaten Material ID found as %s" % self.platen_cells_elem_id)

        all_elem_id = list(set(self.first_file.get_array(self.var_data["mineral_type"])))
        self.rock_elem_ids = [x for x in all_elem_id if x != self.platen_cells_elem_id]
        self.rock_elem_ids_max, self.rock_elem_ids_min = max(self.rock_elem_ids), min(self.rock_elem_ids)

        self.rock_model = (self.first_file.threshold([self.rock_elem_ids_min, self.rock_elem_ids_max], self.var_data["mineral_type"]))
        self.platen = (self.first_file.threshold([self.platen_cells_elem_id, self.platen_cells_elem_id], self.var_data["mineral_type"]))

        sample_x_min, sample_x_max = self.rock_model.bounds[0], self.rock_model.bounds[1]
        sample_y_min, sample_y_max = self.rock_model.bounds[2], self.rock_model.bounds[3]
        sample_z_min, sample_z_max = self.rock_model.bounds[4], self.rock_model.bounds[5]
        self.sample_width = sample_x_max - sample_x_min
        self.sample_height = sample_y_max - sample_y_min
        self.sample_thickness = sample_z_max - sample_z_min

        return self.sample_width, self.sample_height, self.sample_thickness


    def simulation_type(self):
        self.check_edge_point = [self.rock_model.bounds[1], self.rock_model.bounds[3], 0]
        self.check_edge_cell = self.first_file.extract_cells(self.first_file.find_closest_cell(self.check_edge_point))
        self.check_edge_cell = pv.cell_array(self.check_edge_cell, self.var_data['mineral_type'])

        if self.check_edge_cell == -1:
            self.sim_type = "BD Simulation"
        else:
            self.sim_type = "UCS Simulation"
        return self.sim_type

    def platen_info(self, pv_cells, platen_boundary_id, var_property):
        '''
        This function thresholds cells based on boundary condition and sums them based on the defined parameter var_property

        :param pv_cells:
        :type pv_cells: pyvista.core.pointset.UnstructuredGrid
        :param platen_boundary_id: boundary id that the threshold should be based on
        :type platen_boundary_id: float
        :param var_property: name of the property (array to b returned)
        :type var_property: Dict[str, str]
        :return: array of the property based on the threshold
        :rtype: ndarray
        '''

        platen_cell_ids = pv_cells.threshold([platen_boundary_id, platen_boundary_id], self.var_data["boundary"])
        platen_var_prop_list = sum(platen_cell_ids.get_array(var_property))

        if var_property == 'displacement':
            for k in range(0, 3):
                # divide by the number of points per cell  (3 in 2D and 4 in 3D)
                platen_var_prop_list[k] = platen_var_prop_list[k] / (platen_cell_ids.cell_n_points(0) * platen_cell_ids.n_cells)

        return platen_var_prop_list


    # def get_broken(self, mode_id=None):

    # def find_cell(self, points):

    # def model_composition(self):
    # if not self._basic_0_loaded:
    # load_basic(0)
    # def n_cells(self): #TODO: Is this function even useful?

    # def n_points(self): #TODO: Is this function even useful?

    # def broken_joints(self, mode=None,):

    # def seismic_events(self,time_filter = []):

    # def seismic_clustering(self):

    # def set_strain_gauge(self,point,axis):

    def platen_force(self, material_id=None, boundary_condition_id=None, location=None):
        # TODO: load based on threshold points (boundary condition)
        #   can not do as there is a -1 boundary condition! And
        #   this should be limited only to the platen boundary condition (e.g. confinement)
        # TODO: Why do we need so many input parameters ?
        '''
        Checks if Compression Simulations.
        Calculate the stress based on the force in the platens.

        :return: stress (MPa) = force / width of sample
        :rtype: list

        Example:
            >>> data = pv.read("../example_outputs/Irazu_UCS")
            >>> axial_stress = data.platen_force()
            [0.0, 4.825237206318255, 9.628822864052236, 14.414373164820148, 19.191640765444546, 23.95880093999921, 28.711674236346393, 33.44027814633586, 38.03245444023052, 13.266469715550484, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30]
        '''

        history_stress = []  # List to hold the stress_history
        avg_platen_force = [0.0, 0.0, 0.0]  # Dummy cell
        axis_of_loading = 1  # Axis of loading in Y direction.

        ## Get rock dimension.
        self.rock_sample_dimensions()
        ## Check UCS Simulation
        if self.simulation_type() != "UCS Simulation":
            print("Simulation appears to be not for compressive strength")

        # Load each timestep
        def stress_thresholding(openfdem_model_ts):
            platen = (openfdem_model_ts.threshold([self.platen_cells_elem_id, self.platen_cells_elem_id],
                                                      self.var_data["mineral_type"]))
            top, bottom = (platen.get_data_range(self.var_data["boundary"]))

            top_platen_force_list = self.platen_info(openfdem_model_ts, top, self.var_data["platen_force"])
            bot_platen_force_list = self.platen_info(openfdem_model_ts, bottom, self.var_data["platen_force"])

            for i in range(0, self.number_of_points_per_cell):
                # Convert forces from microN to kN and get the average forces
                avg_platen_force[i] = 0.5 * (abs(top_platen_force_list[i]) + abs(bot_platen_force_list[i])) / 1.0e9

            # stress in MPa (force in kN & area in mm^2)
            return avg_platen_force[axis_of_loading] / self.sample_width * 1.0e3
            #history_stress.append(stress)               

        processes=[]
        num_processes = os.cpu_count()

        if name == "__main__":
            for timestep in self:
                process = Process(target=stress_thresholding,args=(timestep))
                processes.append(process)

            for process in processes:
                process.start()

            for process in processes:
                process.join()

        # for openfdem_model_ts in self:
        #
        #     platen = (openfdem_model_ts.threshold([self.platen_cells_elem_id, self.platen_cells_elem_id],
        #                                           self.var_data["mineral_type"]))
        #     top, bottom = (platen.get_data_range(self.var_data["boundary"]))
        #
        #     top_platen_force_list = self.platen_info(openfdem_model_ts, top, self.var_data["platen_force"])
        #     bot_platen_force_list = self.platen_info(openfdem_model_ts, bottom, self.var_data["platen_force"])
        #
        #     for i in range(0, self.number_of_points_per_cell):
        #         # Convert forces from microN to kN and get the average forces
        #         avg_platen_force[i] = 0.5 * (abs(top_platen_force_list[i]) + abs(bot_platen_force_list[i])) / 1.0e9
        #
        #     # stress in MPa (force in kN & area in mm^2)
        #     stress = avg_platen_force[axis_of_loading] / self.sample_width * 1.0e3
        #     history_stress.append(stress)

        return history_stress

    def platen_displacement(self, material_id=None, boundary_condition_id=None, location=None):
        # TODO: load based on threshold points (boundary condition)
        #   can not do as there is a -1 boundary condition! And
        #   this should be limited only to the platen boundary condition (e.g. confinement)
        '''
        Checks if Compression Simulations.
        Calculate the strain based on the displacement in the platens.

        :return: strain (%) = displacement / height of sample
        :rtype: list

        Example:
            >>> data = pv.read("../example_outputs/Irazu_UCS")
            >>> axial_disp = data.platen_displacement()
            [0.0, 0.009259259235881877, 0.018518518471763754, 0.02777777770764563, 0.03703703694352751, 0.046296296179409384, 0.05555555541529126, 0.06481481465117314, 0.07407407388705502, 0.08333333312293689, 0.09259259235881877, 0.10185185159470064, 0.11111111083058252, 0.1203703700664644, 0.12962962930234628, 0.13888888853822817, 0.14814814777411003, 0.15740740700999192, 0.16666666624587378, 0.17592592548175565, 0.18518518471763754, 0.1944444439535194, 0.2037037031894013, 0.21296296242528318, 0.22222222166116504, 0.23148148089704693, 0.2407407401329288, 0.24999999936881068, 0.25925925860469257, 0.2685185178405744, 0.27777777707645634]

        '''

        history_strain = []  # List to hold the stress_history
        avg_platen_disp = [0.0, 0.0, 0.0]  # Dummy cell
        axis_of_loading = 1  # Axis of loading in Y direction.

        ## Get rock dimension.
        self.rock_sample_dimensions()
        ## Check UCS Simulation
        if self.simulation_type() != "UCS Simulation":
            print("Simulation appears to be not for compressive strength")

        # Load each timestep
        #backup copy compressed below

        def strain_thresholding():
            for openfdem_model_ts in self:

                platen = (openfdem_model_ts.threshold([self.platen_cells_elem_id, self.platen_cells_elem_id],
                                                      self.var_data["mineral_type"]))
                top, bottom = (platen.get_data_range(self.var_data["boundary"]))

                avg_top_platen_disp = self.platen_info(openfdem_model_ts, top, self.var_data["platen_displacement"])
                avg_bottom_platen_disp = self.platen_info(openfdem_model_ts, bottom, self.var_data["platen_displacement"])


                for i in range(0, self.number_of_points_per_cell):
                    avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

                strain_from_platen = avg_platen_disp[axis_of_loading] / self.sample_height * 100.0

                history_strain.append(strain_from_platen)


        threads = []
        num_threads = 100 #we can make num_threads increase up to 4 or 5 digits of threads

        for i in range(num_threads):
            thread = Thread(target=strain_thresholding)
            threads.append(thread)

        for thread in threads:
            thread.start()

        for thread in threads:
            thread.join()

        return history_strain

    # def generate_rose_diagrams(self,output_folder):
    # """ Call generate_rose_diagram for each timestep"""

    # def cluster_cracks(self):

    # def save_csv_outputs(self,output_folder):

    # def generate_basic_csv(self,output_file):

    # def generate_brokenjoints_csv(self,output_file):

    # def generate_seismic_csv(self,output_file):

    def generate_history_csv(self, output_file):
        '''

        :param output_file:
        :return:
        '''

        try:
            self.history_stress
        except AttributeError:
            print("Processing UCS")
            Model.process_UCS(self)


    # def generate_seismic_clustering_csv(self,output_file):

    def process_UCS(self):
        '''
        Process model as UCS model

        :return: UCS, stress_hist, strain_hist, Elastic_moduli
        :rtype: tuple[float, list[float, list[float, list[float]]

        Example:
            >>> data = pv.read("../example_outputs/Irazu_UCS")
            >>> UCS, stress_hist, strain_hist, Elastic_moduli = data.process_UCS()
            (38.03245444023052, [0.0, 4.825237206318255, 9.628822864052236, 14.414373164820148, 19.191640765444546, 23.95880093999921, 28.711674236346393, 33.44027814633586, 38.03245444023052, 13.266469715550484, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30, 2.036136892438222e-30], [0.0, 0.009259259235881877, 0.018518518471763754, 0.02777777770764563, 0.03703703694352751, 0.046296296179409384, 0.05555555541529126, 0.06481481465117314, 0.07407407388705502, 0.08333333312293689, 0.09259259235881877, 0.10185185159470064, 0.11111111083058252, 0.1203703700664644, 0.12962962930234628, 0.13888888853822817, 0.14814814777411003, 0.15740740700999192, 0.16666666624587378, 0.17592592548175565, 0.18518518471763754, 0.1944444439535194, 0.2037037031894013, 0.21296296242528318, 0.22222222166116504, 0.23148148089704693, 0.2407407401329288, 0.24999999936881068, 0.25925925860469257, 0.2685185178405744, 0.27777777707645634], [5.15399101160927, 5.181743019752671])
        '''

        self.history_stress, self.history_strain = [], []  # List to hold the history
        avg_platen_force = [0.0, 0.0, 0.0]  # Dummy cell
        avg_platen_disp = [0.0, 0.0, 0.0]  # Dummy cell
        axis_of_loading = 1  # Axis of loading in Y direction.

        # Get rock dimension.
        self.rock_sample_dimensions()
        # Check UCS Simulation
        if self.simulation_type() != "UCS Simulation":
            print("Simulation appears to be not for compressive strength")

        # Load each timestep
        for openfdem_model_ts in self:

            platen = (openfdem_model_ts.threshold([self.platen_cells_elem_id, self.platen_cells_elem_id],
                                                  self.var_data["mineral_type"]))
            top, bottom = (platen.get_data_range(self.var_data["boundary"]))

            top_platen_force_list = self.platen_info(openfdem_model_ts, top, self.var_data["platen_force"])
            bot_platen_force_list = self.platen_info(openfdem_model_ts, bottom, self.var_data["platen_force"])

            avg_top_platen_disp = self.platen_info(openfdem_model_ts, top, self.var_data["platen_displacement"])
            avg_bottom_platen_disp = self.platen_info(openfdem_model_ts, bottom, self.var_data["platen_displacement"])

            for i in range(0, self.number_of_points_per_cell):
                # Convert forces from microN to kN and get the average forces
                avg_platen_force[i] = 0.5 * (abs(top_platen_force_list[i]) + abs(bot_platen_force_list[i])) / 1.0e9
                avg_platen_disp[i] = abs(avg_top_platen_disp[i]) + abs(avg_bottom_platen_disp[i])

            # stress in MPa (force in kN & area in mm^2)
            stress = avg_platen_force[axis_of_loading] / self.sample_width * 1.0e3
            self.history_stress.append(stress)

            strain_from_platen = avg_platen_disp[axis_of_loading] / self.sample_height * 100.0
            self.history_strain.append(strain_from_platen)

        self.UCS = max(self.history_stress)

        index = self.history_stress.index(max(self.history_stress))
        index = int(index / 2)

        # E Tangent
        E_tan = self.E_mod(self.history_stress, self.history_strain, index-1, index+1)
        # E secant
        E_sec = self.E_mod(self.history_stress, self.history_strain, 0, index)

        self.E = [E_tan, E_sec]

        return self.UCS, self.history_stress, self.history_strain, self.E


    def E_mod(self, stress, strain, range_min, range_max):
        '''
        Calculate Elastic modulus between two extents (data points)

        :param stress: Stress in MPa
        :type stress: list[float]
        :param strain: Strain in %
        :type strain: list[float]
        :param range_min: Index of the lower end of range
        :type range_min: int
        :param range_max: Index of the upper end of range
        :type range_max: int
        :return: Elastic Modulus in GPa = delta stress / delta strain
        :rtype: float
        '''

        # Stress in MPa
        d_stress = ( stress[range_max] - stress[range_min] )
        # Convert from strain % to strain
        d_strain = ( strain[range_max] - strain[range_min] ) * 100

        return d_stress / d_strain

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
