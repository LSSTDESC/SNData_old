from __future__ import absolute_import
import os
from .version import __version__
try:
    from . import filters
except:
    print('This requires a THROUGHPUTS directory to be setup as an env var\n')
    print('This may be achieved, for example, by setting up the LSST stack\n')
    print('or by cloning the throughputs directory from git-lfs, \n')
    print('and setting up env vars. The function is to provide the LSST bandpass\n') 
    print('files\n')
from .aliases import *
from .lightcurve import *
from .photometry import *
from .snanaio import *
from .analyzelcFits import *

here = __file__
basedir = os.path.split(here)[0]
example_data = os.path.join(basedir, 'example_data')
