from __future__ import absolute_import
import os
try:
    from . import filters
except:
    print('This requires a THROUGHPUTS directory to be setup as an env var\n')
    print('This may be achieved, for example, by setting up the LSST stack\n')
    print('or by cloning the throughputs directory from git-lfs, \n')
    print('and setting up env vars. The function is to provide the LSST bandpass\n') 
    print('files\n')
from .version import __version__
from .aliases import *
from .lightcurve import *
from .photometry import *
from .snanaio import *
from . import aliases, lightcurve, snanaio
here = __file__
basedir = os.path.split(here)[0]
example_data = os.path.join(basedir, 'example_data')
__all__ = aliases.__all__ + lightcurve.__all__ + snanaio.__all__ \
        + ['example_data']
