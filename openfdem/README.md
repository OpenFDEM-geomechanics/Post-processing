# OpenFDEM Post-Processing Package
### University of Toronto, 2021

### Support
This Python package performs transformations on hybrid finite-discrete element method (FDEM) models with an unstructured grid in vtk/vtu/vtp format. It currently supports arrays of simulation files from Irazu, Y-Geo (and its common derivatives), as well as OpenFDEM.

###Functionality
Functions include:
- Extract information within the FDEM Model based on the name of the array (e.g., Stress, Strain, Temperature, etc...)
- Extract stress-strain information for UCS and BD Simulations. Optional addition of virtual strain gauges.
- Plotting stress vs strain curves
- Calculate the Elastic Modulus of the dataset. Eavg, Esec and Etan can be evaluated. 
- Extract information of a particular cell based on a sequence of array names. This can be extended to extracting information along a line.
- Extract information of a threshold dataset criteria based on a sequence of array names.
- Extraction of information for cells/thresholds extended to 3D Models as well. 

Please refer to the user manual for more in-depth on the various functions and their usage. 

# Installation
1. Ensure that you have a Python version 3.5-3.9 installed on your machine.
2. Execute the following commands (*if running from test pypi only, run one at the time*)
   ```cmd
   pip install pyvista
   pip install pandas
   pip install tqdm
   pip install h5py
   pip install windrose
   ```
3. Navigate to terminal/cmd and execute the following command:
   ```cmd
   pip install -i https://test.pypi.org/simple/ openfdem
   ```
4. Check that installation was successful by running:
   ```cmd
   pip freeze
   ```

## How to use it - Terminal
1. After installation is complete, in a terminal/cmd session where you keep your Python projects, type the following command:
    ```cmd
    python
    ```
2. In terminal, import openfdem module like so: 
    ```python
    import openfdem as fd
    ```
3. Specify the path to folder where your models (.vtk, .vtu files) are located on your machine.
    ```python
    model = fd.Model("abs_model_path_on_machine")
    ```
4. You are ready to start using openfdem! See the operating module for full description of functions available.

## Example OpenFDEM Functions
1. Getting number of points in your model:
   ```py
   >>> import openfdem as fd
   >>> model = fd.Model("abs_model_path_on_machine")
   >>> model.n_points
   11904
   ```

2. Getting the stress-strain for platen:
   ```python
   >>> import openfdem as fd
   >>> model = fd.Model("abs_model_path_on_machine")
   >>> df = model.complete_UCS_stress_strain(progress_bar=True)
   Script Identifying Platen
       Platen Material ID found as [1]
   Progress: |//////////////////////////////////////////////////| 100.0% Complete
   1.51 seconds
   >>> df
        Platen Stress  Platen Strain
    0    0.000000e+00       0.000000
    1    4.825237e+00       0.009259
    2    9.628823e+00       0.018519
    3    1.441437e+01       0.027778
    4    1.919164e+01       0.037037
    ..            ...            ...
    57   2.036137e-30       0.240741
    58   2.036137e-30       0.250000
    59   2.036137e-30       0.259259
    60   2.036137e-30       0.268519
    61   2.036137e-30       0.277778
    
    [62 rows x 2 columns]
   ```


