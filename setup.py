import setuptools
import pathlib
import os

HERE = pathlib.Path("pyfdempp/README.md")
README = (HERE).read_text()

setuptools.setup(
    name="pyfdempp",
    author="Grasselli's Geomechanics Group - University of Toronto",
    author_email="",
    description="Performs transformations on hybrid finite-discrete element method (FDEM) models with an unstructured grid in vtk/vtu/vtp format.",
    long_description=README,
    long_description_content_type="text/markdown",
    keywords="example documentation tutorial",
    url="https://openfdem.com/",
    packages=setuptools.find_packages(),
    python_requires='>=3.5',
    version=0.1,
    install_requires=[
        "pandas>=0.0",
        "numpy>=1.0",
        "pyvista>=0.20",
        "joblib",
        "tqdm"
    ]

)