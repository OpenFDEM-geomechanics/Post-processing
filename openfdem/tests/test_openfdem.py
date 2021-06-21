import sys
sys.path.insert(0,'./')

import pytest
import openfdem as fd

data1 = fd.Model("example_outputs/Irazu_UCS",fdem_engine="Irazu")
data2 = fd.Model("example_outputs/OpenFDEM_BD/result",fdem_engine="OpenFDEM")
data3 = fd.Model("example_outputs/Y-GEO_fluid_example",fdem_engine="OpenFDEM")

# Function to check if files are loaded
def check_files(data, n_files):
    assert len(data._basic_files) == n_files 
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

def test_file_not_loaded_openfdem():
    with pytest.raises(LookupError) as e_info:
        # Invalid path to model
        assert fd.Model("Irazu_UCS",fdem_engine="Irazu")