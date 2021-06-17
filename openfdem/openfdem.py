# Possible import statements to organize the code
#import _seismic_func
#import _fracture_func
#import _modelling_func
#import _data_func

class Model:
	""" Model class collects datafiles into one interface
	Each data array returns as a list ordered by timestep
	Collection of timesteps?
	handles temporal manipulations
	"""
	_file_names = {	"OpenFDEM":{"basic":"_field_",
							},
					"Irazu":{"basic":"_basic_"}
					}
	def _find_output_files(dir_path,file_extension,output_type):
		os.chdir(dir_path) # cd to the base directory path
		list_of_files = []
		# Get the list of files/directories
		files = os.listdir('.')
		files.sort()
		for f in files:
		# Filter only the files
			if os.path.isfile(f):
			# Clean the list to only include files with the given extension that contain output_type
				if f.endswith(file_extension):
					if output_type in f:
						list_of_files.append(os.path.join(dir_path, f))
		# Sort list of files by their numerical value
		list_of_files = sorted(list_of_files, key=numericalSort)
		return list_of_files	
	
	def __init__(self,folder=None,runfile=None,fdem_):
		""" Create a Model object that finds all the run files and organizes all the file names on creation
		
		Example: 
		>>> data = openfdem.Model("Model folder")
		"""
		self._folder = folder
		self._runfile = runfile
		if runfile.endswith(".y"):
			self._fdem_engine = "OpenFDEM"
		else if runfile.endswith(".fdem"):
			self._fdem_engine = ""
		if folder is not None:
			extensions = tuple("vtu","vtp")
			self._basic_files = _find_output_files(folder,extensions,"_basic_")
			self._broken_files = _find_output_files(folder,extensions,"_basic_")
		
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
	""" A class handling the data of each timestep
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
	
	