# Possible import statements to organize the code
#import _seismic_func
#import _fracture_func
#import _modelling_func
#import _data_func

import os, re, glob
import os.path as path

import pyvista as pv

class Model:
    """Model class collects datafiles into one interface.
    
    Each data array returns as a list ordered by timestep
    Collection of timesteps?
    handles temporal manipulations
    
    Example:
        >>> import openfdem as fdem
        >>> model = fdem.Model("../example_outputs/Irazu_UCS")
        
    """
    
    _file_names = {    "OpenFDEM":  {"basic":"_field_",
                                     "broken_joint":"_brokenjoint_",
                                     "soften_joint":"_softenjoint_",
                                     "principal_stress_direction":"_principalstress_",
                                     "acoustic_emission":"acousticemission"
                                     },
                        "Irazu":    {"basic":"_basic_",
                                     "broken_joint":"_broken_joint_",
                                     "soften_joint":"_softened_joint_",
                                     "principal_stress_direction":"_principal_stress_strain_",
                                     "acoustic_emission":"_seismic_event_"
                                     }
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

    def _findOutputFiles(self,dir_path,file_extension,output_type):
        """Find the collection of vtk files that belong to a specific group
        
        :param dir_path: Starting directory (full path is required)
        :type dir_path: str
        :return: A list of subdirectories
        :rtype: list[str]
        
        Example:
            >>> data = fdem.Model("../example_outputs/Irazu_UCS")
            >>> data._findOutputFiles("../example_outputs/Irazu_UCS/",tuple(["vtu","vtp"]),"basic")[0:3]
            ['..\\\\example_outputs\\\\Irazu_UCS\\\\UCS_tutorial-Run_1_femdem.r2m_basic_0.vtu', '..\\\\example_outputs\\\\Irazu_UCS\\\\UCS_tutorial-Run_1_femdem.r2m_basic_20000.vtu', '..\\\\example_outputs\\\\Irazu_UCS\\\\UCS_tutorial-Run_1_femdem.r2m_basic_40000.vtu']
        """
        
        list_of_files = []
        # Get the list of files/directories
        for ext in file_extension:
            list_of_files.extend(glob.glob(dir_path+"/*"+Model._file_names[self._fdem_engine][output_type]+"*."+ext))
        
        # Sort list of files by their numerical value
        list_of_files = sorted(list_of_files, key=Model._numericalSort)
        
        # Ensure it is in proper system path format
        list_of_files = [path.relpath(vtkfile) for vtkfile in list_of_files]
        
        return list_of_files    
    
    def __init__(self,folder=None,runfile=None,fdem_engine=None):
        """Create a Model object that finds all the run files and organizes all the file names on creation.
        
        Example:: 
            >>> data = openfdem.Model("Model folder")
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
            extensions = tuple(["vtu","vtp"])
            self._basic_files = self._findOutputFiles(folder,extensions,"basic")
            if len(self._basic_files) == 0:
                raise LookupError('Folder does not appear to exist or does not have valid output files.')

            self._broken_files = self._findOutputFiles(folder,extensions,"broken_joint")
            self._soft_files = self._findOutputFiles(folder,extensions,"soften_joint")
            self._pstress_dir_files = self._findOutputFiles(folder,extensions,"principal_stress_direction")
            self._acoustic_files = self._findOutputFiles(folder,extensions,"acoustic_emission")

            first_file = pv.read(self._basic_files[0])
            self.n_timesteps = len(self._basic_files)
            self.n_points = first_file.points.shape[0]
            self.n_elements = first_file.n_cells

    # def __getitem__(self,key):
        # """ For timestep access
        # Timestep accessible by index or string representing timestep
        # timestep = data[0]
        # timestep = data['180000']
        # """
        
    # def __contains__(self,item):
        # """ Access to timestep in Model
        # """
    # def __len__(self):
        # """ Number of timestep files """
        
    # def load_basic(self,timestep=None):
        
    # def filter_material(self, material_id):
        
    # def get_broken(self, mode_id=None):
    
    # def find_cell(self, points):
    
    # def model_composition(self):
        # if not self._basic_0_loaded:
            # load_basic(0)
    # def n_cells(self):
    
    # def n_points(self):
        
    # def broken_joints(self, mode=None,):
    
    # def seismic_events(self,time_filter = []):
    
    # def seismic_clustering(self):
    
    # def set_strain_gauge(self,point,axis):
    
    # def platen_force(self,material_id=None,boundary_condition_id=None,location=None):
    
    # def platen_displacement(self,material_id=None, boundary_condition_id=None,location=None):
    
    # def generate_rose_diagrams(self,output_folder):
        # """ Call generate_rose_diagram for each timestep"""
        
    # def cluster_cracks(self):
    
    # def save_csv_outputs(self,output_folder):
    
    # def generage_basic_csv(self,output_file):
    
    # def generate_brokenjoints_csv(self,output_file):
    
    # def generate_seismic_csv(self,output_file):
    
    # def generate_history_csv(self,output_file):
    
    # def generate_seismic_clustering_csv(self,output_file):
    
    # def process_UCS(self):
        # """ Process model as UCS model
        
        # Example: 
        # >>> UCS,stress_hist,strain_hist,Elastic_moduli = model.process_UCS()
        # """
    
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
    
    