# /////////////////////////////////////////////////////////////// #
# !python3.6
# -*- coding: utf-8 -*-
# Python Script initially created on 6/16/2021
# Compiled by Aly @ Grasselli's Geomechanics Group, UofT, 2021
# Created using PyCharm // Tested on Spyder
# Current Version 06 - Dated August 21, 2018
# /////////////////////////////////////////////////////////////// #

import platform
import sys
import os

'''
VERIFICATION CODE BLOCK
    - Verifies x64 Architecture
    - Check Python compatibility
    - Display system information
'''
if platform.architecture()[0] != "64bit":
    exit("Compatible only on x64")

# Check python compatibility before proceeding
try:
    assert sys.version_info >= (3, 5) and sys.version_info <= (3, 8)
    print("Python Version: %s" % sys.version.split('\n')[0])
except AssertionError:
    exit("Compatible Python Version 3.5+ upto 3.8.x")

cwd = os.getcwd()

print("Computer Running %s, Release %s, Version %s" % (platform.system(), platform.release(), platform.version()))


import processing

fold = os.path.join(cwd, "example_outputs", "Irazu_UCS")

'''
STRESS FROM PLATENS
'''
platen_stress = processing.UCS_platen_Stress(fold)

'''
STRAIN FROM PLATENS
'''
platen_strain = processing.UCS_platen_Strain(fold)

'''
SIMPLE PLOT FUNCTION
'''

import matplotlib.pyplot as plt

plt.plot(platen_strain, platen_stress)
plt.show()