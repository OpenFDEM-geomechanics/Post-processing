import setuptools
import pathlib
import os

HERE = pathlib.Path("/Users/ekaterinaossetchkina/PycharmProjects/OpenFDEM-Post-Processing/openfdem/README.md")
README = (HERE).read_text()

setuptools.setup(
    name="openfdem",
    author="University of Toronto, 2021",
    description="Wheel file of openfdem modules",
    long_description=README,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    python_requires='>=3.5,<=3.8',
    version=3.5,
    install_requires=[
        "pandas>=0.0",
        "numpy>=1.0",
        "h5py>=2",
        "pyvista>=0.20",
        "joblib",
        "tqdm"
    ]

)