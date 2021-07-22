import os
import glob
import copy

from vtkmodules.util import numpy_support
from vtkmodules.numpy_interface import dataset_adapter as dsa
import pyvista as pv
import h5py
from tqdm import tqdm



class aggregate_storage:
    def file_group_key(self,vtkfilename):
        basename = os.path.basename(vtkfilename)
        timestep = os.path.splitext(basename)[0].split("_")[-1]
        return timestep+"/"+basename

    def _store_file_init(self,f,vtkfilename):
        """Protected function to store vtk file data into h5 file

        :param f: H5 file handle with write permissions enabled
        :type f: [type]
        :param vtkfilename: Path to vtk file
        :type vtkfilename: str
        """
        groupkey = self.file_group_key(vtkfilename)
        data = pv.read(vtkfilename)
        array_names = data.array_names
        points = data.points
        cells = data.cells
        offsets = data.offset
        celltypes = data.celltypes
        ds = f.create_dataset(groupkey+"/points",points.shape,points.dtype)
        ds[()] = points
        ds = f.create_dataset(groupkey+"/cells",cells.shape,cells.dtype)
        ds[()] = cells
        ds = f.create_dataset(groupkey+"/offsets",offsets.shape,offsets.dtype)
        ds[()] = offsets
        ds = f.create_dataset(groupkey+"/celltypes",celltypes.shape,celltypes.dtype)
        ds[()] = celltypes
        nondataarrays = ["cells","offsets","celltypes",'points']

        for array_name in array_names:
            if not any(x in array_name for x in nondataarrays):
                array = data.get_array(array_name)
                ds = f.create_dataset(groupkey+"/"+array_name,array.shape,array.dtype)
                ds[()] = array

    def store_file(self,vtkfilename):
        f_w = h5py.File(self.h5filename,"w")
        self._store_file_init(f_w,vtkfilename)
    
    def _check_files_stored(self):
        try:
            f = h5py.File(self.h5filename,'r')
        except:
            return False
        timesteps = f.keys()
        files = [os.path.basename(x) for x in self.files]
        
        for ts in timesteps:
            ts_files = f[ts].keys()
            for ts_file in ts_files:
                if ts_file in files:
                    files.remove(ts_file)
                else:
                    return False

        if not files:
            return True
    def __init__(self,file_directory,h5filename=None,overwrite=False,verbose=False):
        if h5filename is None:
            h5filename = os.path.join(file_directory,"aggregated.h5")
        self.h5filename = h5filename
        
        self.file_directory = file_directory
        self.files = glob.glob(file_directory+"/*.vtu") + glob.glob(file_directory+"/*.vtp")
        storage_complete = False
        if os.path.exists(h5filename):
            if overwrite:
                os.remove(h5filename)
            else:
                storage_complete = self._check_files_stored()
        
        if storage_complete:
            if verbose:
                print("Files have already been stored.. Skipping storage process")
        else:
            if verbose:
                print(f"Selected directory {file_directory}... {len(self.files)} files found")
                print(f"Storing data in {h5filename}")
            f_w = h5py.File(h5filename,"w")
            
            if verbose:
                iterable = tqdm(self.files)
            else:
                iterable = self.files

            for file_name in iterable:
                if verbose:
                    iterable.set_description(f"Storing {os.path.basename(file_name)}")
                self._store_file_init(f_w,file_name)

    def read_file(self,filename,verbose=False):
        if verbose:
            print(f"Reading {filename}\n")

        groupkey = self.file_group_key(filename)
        try:
            f = h5py.File(self.h5filename,'r',swmr=True)
        except:
            raise
        group = f[groupkey]

        data = pv.UnstructuredGrid(group['cells'][()],group['celltypes'][()],group['points'][()])
        nondataarrays = ["cells","offsets","celltypes",'points']
        npoints = group['points'].shape[0]
        for name in group.keys():
            if not any(x in name for x in nondataarrays):
                array = group[name]
                if array.shape[0] == npoints:
                    data.point_arrays.append(group[name],name)
                else:
                    data.cell_arrays.append(group[name],name)
        return data 
if __name__ == '__main__':
    # # filename = r"S:\20210623_BD_DIC\20210622_BD_test_DIC_clean_working_4\results\DICe_solution.e"
    source_path = r'X:\20210622_irazu_example\UCS_example'
    filename = r'X:\20210622_irazu_example\UCS_example\UCS_femdem.r2m_basic_2000.vtu'
    # test_model = pv.read(filename)
    # test_model.cell_arrays
    # print(test_model.array_names)
    # print(test_model)
    # # pv.plot(test_model)

    storage = aggregate_storage(source_path,verbose=True)

    reformed_model = storage.read_file(filename)
    reformed_model.save(os.path.join(source_path,"test.vtu"))
    # print(reformed_model.array_names)
    # reformed_model.plot(scalars='displacement')

