# pyfdempp - An FDEM Visualiser Post-Processing Package

## About

This Python package performs transformations on hybrid finite-discrete element method (FDEM) models with an unstructured grid in vtk/vtu/vtp format. It currently supports arrays of simulation files from the FDEM solvers:
- [Geomechanica's Irazu software](https://www.geomechanica.com/software/),
- [Y-Geo](https://geogroup.utoronto.ca/software/) (and its common derivatives), as well as 
- [OpenFDEM](https://openfdem.com/html/index.html).

The package is heavily dependent on `pyvista` and is limited to `Python >=3.5`. The package is maintained by the [Grasselli's Geomechanics Group](https://geogroup.utoronto.ca/) at the University of Toronto, Canada, and is part of a collaborative effort by the open-source pacakge [OpenFDEM](https://openfdem.com/html/index.html).

## Functionality

This Python package performs transformations on FDEM models in vtk/vtu format. 
Functions calculated include:
- Extract information within the FDEM Model based on the name of the array (e.g., Stress, Strain, Temperature, etc...). Works in 2D and 3D.
- Extract stress-strain information for UCS and BD Simulations (Works in 2D and 3D). Optional addition of virtual strain gauges (Limited to 2D).
- Plotting stress vs strain curve.
- Calculate the Elastic Modulus of the dataset. Eavg, Esec and Etan can be evaluated. Works in 2D and 3D.
- Extract information of a particular cell based on a sequence of array names. This can be extended to extracting information along a line. Works in 2D and 3D.
- Extract information of a threshold dataset criteria based on a sequence of array names. Works in 2D and 3D.
- Extract mesh information and plot rosette/polar plots. Limited to 2D.
- Automatic detection/ User-defined assigment of loading direction when analysing mechanical simulations, namely UCS, BD, and PLT, in both 2D and 3D simulations.

## How to use it - Terminal

1. After installation is complete, in a terminal/cmd session where you keep your Python projects, type the following command:
    ```python
    python
    ```
2. In terminal, import the module like so:
    ```python 
    import pyfdempp as fd
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
       >>> Script Identifying Platen
       >>> Platen Material ID found as [1]
       >>> Progress: |//////////////////////////////////////////////////| 100.0% Complete
       >>> 1.51 seconds
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



