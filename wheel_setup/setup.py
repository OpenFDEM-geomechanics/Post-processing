import setuptools
import pathlib
import os

HERE = pathlib.Path("About.md")
README = (HERE).read_text()

setuptools.setup(
    name="fdem_visualizer",
    author="University of Toronto, 2022",
    description="Wheel file of fdem visualizer modules",
    long_description=README,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    python_requires='>=3.5,<=3.9',
    version=1.0,
    install_requires=[
        "pandas>=0.0",
        "numpy>=1.0",
        "h5py>=2",
        "pyvista>=0.20",
        "joblib",
        "tqdm"
    ]

)