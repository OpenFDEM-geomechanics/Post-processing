import platform
import sys

'''
VERIFICATION CODE BLOCK
    - Verifies x64 Architecture
    - Check Python compatibility
    - Display system information
'''

__version__ = "0.0.1"
if platform.architecture()[0] != "64bit":
    exit("Compatible only on x64")

# Check python compatibility before proceeding
try:
    assert sys.version_info >= (3, 5) and sys.version_info <= (3, 9)
    print("Python Version: %s" % sys.version.split('\n')[0])
except AssertionError:
    print("Python Version: %s" % sys.version.split('\n')[0])
    exit("Compatible Python Version 3.5+ upto 3.8.x")

from .openfdem import (
    Model,
    Timestep
    )