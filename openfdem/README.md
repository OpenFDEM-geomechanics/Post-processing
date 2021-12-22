# OpenFDEM Post-Processing Package
### University of Toronto, 2021

This Python package performs transformations on FDEM models in vtk/vtu format. 
Functions calculated include:
- Stress
- Strain
- Plotting stress vs strain curves
- Tan modulus
- Shear

# Installation
1. Ensure that you have a Python version 3.5-3.9 installed on your machine.
2. Execute the following commands (*if running from test pypi only, run one at the time*)
   ```python
   pip install pyvista
   ```
   ```python
   pip install pandas
   ```
   ```python
   pip install tqdm
   ```
   ```python
   pip install h5py
   ```   
3. Navigate to terminal/cmd and execute the following command:
   ```python
   pip install -i https://test.pypi.org/simple/ openfdem
   ```
4. Check that installation was successful by running:
   ```python
   pip freeze
   ```

## How to use it - Termimal
1. After installation is complete, in a terminal/cmd session where you keep your Python projects, type the following command:
    ```python
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
    ```python
    model.n_points
    ```
   Output:
   ```python
   11904
   ```

2. Getting the stress-strain for platen:
   ```python
       model.complete_stress_strain(progress_bar=True)
   ```
   Output:
   ```python
       >>>Script Identifying Platen
       >>>Platen Material ID found as [1]
       >>>Progress: |//////////////////////////////////////////////////| 100.0% Complete
       >>>1.51 seconds
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
       >>> 
   ```

