from __future__ import absolute_import, print_function, division

import sndata as ans
import os
import pandas as pd
from pandas.util.testing import assert_frame_equal
MegacamBandNames = list(x.encode() for x in 'ugriz')
registeredMegaCamBands = tuple('megacam' + band.decode() for band in MegacamBandNames)
def test_load():
    headFile = os.path.join(ans.__path__[0], 'example_data',
                            'snana_fits_HEAD.FITS')
    photFile = os.path.join(ans.__path__[0], 'example_data',
                            'snana_fits_PHOT.FITS')
    sne = ans.SNANASims(headFile=headFile, photFile=photFile,
                        coerce_inds2int=False,
                        SNANABandNames=MegacamBandNames,
                        registeredBandNames=registeredMegaCamBands)
    print(sne.bandNames)
    assert sne.bandNames == MegacamBandNames
    assert sne.newbandNames == ('megacamu', 'megacamg', 'megacamr', 'megacami',
                                'megacamz')
    assert len(sne.bandNameDict.keys()) == 5 
    assert len(sne.headData) == 2
    assert len(sne.get_SNANA_photometry(snid='03d1aw').lightCurve) > 0
    assert len(sne.get_SNANA_photometry(snid='03d1ax').lightCurve) > 0
def test_snanadatafiles():
    headFile = os.path.join(ans.__path__[0], 'example_data',
                            'snana_fits_HEAD.FITS')
    photFile = os.path.join(ans.__path__[0], 'example_data',
                            'snana_fits_PHOT.FITS')
    loc = os.path.join(ans.__path__[0], 'example_data')
    testheadFile = ans.SNANASims.snanadatafile(snanafileroot='snana_fits',
                                               filetype='head',
                                               location=loc)
    testphotFile = ans.SNANASims.snanadatafile(snanafileroot='snana_fits',
                                               filetype='pHot',
                                               location=loc)
    assert testheadFile == headFile
    assert testphotFile == photFile
def test_fromSNANAfileroot():
    loc = os.path.join(ans.__path__[0], 'example_data')
    headFile = os.path.join(ans.__path__[0], 'example_data', 'snana_fits_HEAD.FITS')
    photFile = os.path.join(ans.__path__[0], 'example_data', 'snana_fits_PHOT.FITS')

    sne = ans.SNANASims(headFile=headFile, photFile=photFile, coerce_inds2int=False)
    test_sne = ans.SNANASims.fromSNANAfileroot(snanafileroot='snana_fits',
                                               location=loc)
    assert_frame_equal(test_sne.headData, sne.headData)
