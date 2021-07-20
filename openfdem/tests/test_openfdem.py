import sys
sys.path.insert(0,'./')

import numpy as np
import pytest
import openfdem as fd

data1 = fd.Model("example_outputs/Irazu_UCS",fdem_engine="Irazu")
data2 = fd.Model("example_outputs/OpenFDEM_BD/result",fdem_engine="OpenFDEM")
data3 = fd.Model("example_outputs/Y-GEO_fluid_example",fdem_engine="OpenFDEM")

# Function to check if files are loaded
def check_files(data, n_files):
    assert len(data._basic_files) == n_files 
    
    # Since other files currently don't exist for openfdem outputs, allow other filetypes to have none
    assert len(data._broken_files) == n_files or len(data._broken_files) == 0
    assert len(data._soft_files) == n_files or len(data._soft_files) == 0
    assert len(data._pstress_dir_files) == n_files or len(data._pstress_dir_files) == 0
    assert len(data._acoustic_files) == n_files or len(data._acoustic_files) == 0

def test_file_loading_irazu():
    print(sys.path)
    check_files(data1,31)

def test_file_loading_openfdem():
    check_files(data2,20)
    check_files(data3,3)

# Tests if error occurs with path containing now input files
def test_file_not_loaded_openfdem():
    with pytest.raises(LookupError):
        # Invalid path to model
        assert fd.Model("Irazu_UCS",fdem_engine="Irazu")

# Metadata tests
def test_model_timesteps():
    assert data1.n_timesteps == 31
    assert data2.n_timesteps == 20
    assert data3.n_timesteps == 3

def test_model_counts():
    assert data1.n_points == 11904
    assert data2.n_points == 9000
    assert data3.n_points == 726

    assert data1.n_elements == 3968
    assert data2.n_elements == 3000
    assert data3.n_elements == 242

def test_model_simtype():
    assert data1.model_domain()//2 == 2
    assert data2.model_domain()//2 == 2
    assert data3.model_domain()//2 == 2

def test_model_engine():
    assert data1._fdem_engine == "Irazu"
    assert data2._fdem_engine == "OpenFDEM"
    assert data3._fdem_engine == "OpenFDEM"

def test_model_dims_irazu():
    assert data1.model_dimensions() == (56.0,116.0,0.0)
    assert data1.model_dimensions(mat_id=0) == (52.0,108.0,0.0)
    assert data1.model_dimensions(mat_id=1) == (56.0,116.0,0.0)
    with pytest.raises(Exception):
        assert data1.model_dimensions(mat_id=2)

def test_model_dims_openfdem():
    assert data2.model_dimensions() == (100.0,150.0,0.0)
    assert data2.model_dimensions(mat_id=0) == (100.0,100.0,0.0)
    assert data2.model_dimensions(mat_id=1) == (100.0,150.0,0.0)

    with pytest.raises(Exception):
        assert data2.model_dimensions(mat_id=2)

    assert data3.model_dimensions() == (1.0,1.0,0.0)
    assert data3.model_dimensions(mat_id=0) == (1.0,1.0,0.0)

    with pytest.raises(Exception):
        assert data3.model_dimensions(mat_id=1)

def test_model_rock_dim():
    assert data1.rock_sample_dimensions() == (52.0,108.0,0.0)

def test_model_timestep_access():
    assert np.array_equal(data1[1].get_array(data1.var_data["platen_force"]), data1['20000'].get_array(data1.var_data["platen_force"]))

    with pytest.raises(Exception):
        assert np.array_equal(data1[2].get_array(data1.var_data["platen_force"]), data1['20000'].get_array(data1.var_data["platen_force"]))

def test_model_width():
    assert data1.model_dimensions() == (56.0,116.0,0.0)
    assert data2.model_dimensions() == (100.0,100.0,0.0)
    assert data3.model_dimensions() == (1.0,1.0.0,0)

def test_model_platen_force():
    assert data1.platen_force()[0] == 0.0
    assert data1.platen_force()[-1] == 2.036136892438222e-30
    assert data1.platen_force()[9] == 38.03245444023052
# def test_platen_info():
#     data1.platen_info()