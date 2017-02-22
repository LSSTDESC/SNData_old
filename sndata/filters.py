#!/usr/bin/env python
"""
Module executed an install and setup to register both LSST and MEGACAM filters
as SNCosmo bands.

This allows SNCosmo to automatically use the 
    LSST bandpass objects through the string 'lsst_b' where b is one of 'ugrizy'
    Megacam bandpass objects through the string 'megab' where b is one of 'ugrizy'

Note : The megacam bands added are the 'average' megacam bands using in ugriz.
1. Post 2007 (June), the i band has been changed to i2. This is important for
    SNLS5, but SNLS3 was taken prior to that. 
2. For the precise analysis of MEGACAM SN, one needs to take into account the
    dependence of the radial position of the SN. This is not done here.
"""
from __future__ import absolute_import

import os
import numpy as np
from astropy.units import Unit
import sncosmo

lsstbandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
lsstbanddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'baseline')
megacamPassList = 'ugriz'
megacambanddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'megacam')
# lsstbands = list()
# lsstbp = dict()

for band in lsstbandPassList:

    # setup sncosmo bandpasses
    bandfname = lsstbanddir + "/total_" + band + '.dat'


    # register the LSST bands to the SNCosmo registry
    # Not needed for LSST, but useful to compare independent codes
    # Usually the next two lines can be merged,
    # but there is an astropy bug currently which affects only OSX.
    numpyband = np.loadtxt(bandfname)
    sncosmoband = sncosmo.Bandpass(wave=numpyband[:, 0],
                                   trans=numpyband[:, 1],
                                   wave_unit=Unit('nm'),
                                   name='lsst' + band)

    sncosmo.registry.register(sncosmoband, force=True)

for band in megacamPassList:
    bandfname = os.path.join(megacambanddir, band + 'Mega.fil.txt')
    numpyband = np.loadtxt(bandfname)
    sncosmoband = sncosmo.Bandpass(wave=numpyband[:, 0],
                                   trans=numpyband[:, 1],
                                   wave_unit=Unit('nm'),
                                   name='megacam' + band)
    sncosmo.registry.register(sncosmoband, force=True)
