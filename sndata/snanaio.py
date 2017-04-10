from __future__ import absolute_import, print_function, division
import fitsio
import pandas as pd
import numpy as np
import os
import gzip
from .lightcurve import LightCurve
from astropy.table import Table, Column
from collections import defaultdict

__all__ = ['SNANASims','SNChalSims']

lsst_bandNames = 'ugrizY'
lsst_bandpassNames = tuple('lsst' + band
                           for band in lsst_bandNames.lower())
class SNANASims(object):
    """
    Class to represent SNANA simulations of a particular class of objects, ie.
    Ia or Non_Ia
    """
    def __init__(self, headFile, photFile, coerce_inds2int=True,
                 SNANABandNames=lsst_bandNames,
                 registeredBandNames=lsst_bandpassNames):
        """
    Parameters
    ---------
    headFile : string, mandatory
        absolute path to head file of simulation
    photFile : string, mandatory
        absolute path to phot file of simulation
    coerce_inds2int : Bool, optional, defaults to True
        if true, converts SNID from string to int
        SNANABandNames : iterable of strings/characters, optional, defaults to LSST
            characters used to denote the bandpass in the SNANA simulations
        registeredBandNames : iterable of strings, optional, defaults to LSST
            names of the bands registered in SNCosmo
    """
        self.headFile = headFile
        self.photFile = photFile
        self.headData = self.get_headData(self.headFile,
                      coerce_inds2int=coerce_inds2int)
        self.phot = fitsio.FITS(photFile)
        self.bandNames = SNANABandNames
        self.newbandNames = registeredBandNames
        self.bandNameDict = dict(zip(self.bandNames, self.newbandNames)) 

    @classmethod
    def fromSNANAfileroot(cls, snanafileroot, location='./',
                          coerce_inds2int=False, 
                          SNANABandNames=lsst_bandNames,
                          registeredBandNames=lsst_bandpassNames):
        """
        Class constructor from a root file and a location

        Parameters
        ----------
        snanafileroot : string, mandatory
            root file name for the SNANA which is the prefix to
            '_HEAD.FITS', or '_PHOT.FITS'
        location : string, optional defaults to current working directory './' 
            Relative or absolute path to the directory where the head and phot
            files are located
        snids : integer/string, optional defaults to None
            if not None, only SN observations corresponding to SNID snid
            are loaded
        n : Integer, defaults to None
            if not None, only the first n SN light curves are loaded
        """

        headfile = cls.snanadatafile(snanafileroot, filetype='head',
                                     location=location)
        photfile = cls.snanadatafile(snanafileroot, filetype='phot',
                                     location=location)
        return cls(headFile=headfile, photFile=photfile,
                   coerce_inds2int=coerce_inds2int, SNANABandNames=SNANABandNames,
                   registeredBandNames=registeredBandNames)

    @staticmethod
    def snanadatafile(snanafileroot, filetype='head', location='./'):
        '''
        obtain the name of the head or phot file of an SNANA simulation
        and dataset

        Parameters
        ----------
        snanafileroot : string, mandatory
            root file name for the SNANA which is the prefix to
            '_HEAD.FITS', or '_PHOT.FITS'
        filetype : string, optional defaults to 'head'
            'head' or 'phot' depending on whether a summary file or a photometry
            file is being used.
        location : string, optional defaults to current working directory './' 
            relative or absolute path to the directory in which the file is
            located

        Returns
        -------
            string : absolute path to the SNANA file 

        '''

        desiredfiletype = ['head', 'phot']
        filetype = filetype.lower()
        if not filetype in desiredfiletype:
            raise ValueError(
                'filetype should be one of "head" or "phot"', filetype)
        location = os.path.abspath(location)
        suffix = '_HEAD.FITS'
        if filetype.lower() == 'phot':
            suffix = '_PHOT.FITS'
        fname = snanafileroot + suffix
        return os.path.join(location, fname)

    @staticmethod 
    def get_headData(headFile, coerce_inds2int=False):
        """
    read the headData of a SNANA simulation and return a dataframe
    representing the simulation
        
        Parameters
    ----------
    headFile :
    coerce_inds2int :
    """
        _head = Table.read(headFile)
        if _head['SNID'].dtype.type is np.string_:
            data = _head['SNID'].data
            name = _head['SNID'].name
            dtype = _head['SNID'].dtype
            if coerce_inds2int:
                arr = list(np.int(x) for x in data)
                dtype=int
            else:
                arr = list(x.strip().lower() for x in data)
            col = Column(data=arr, name=name, dtype=dtype)
            _head.remove_columns('SNID')
            _head.add_column(col, index=0)
        return  _head.to_pandas().set_index('SNID')
        
    def get_photrows(self, row=None, snid=None):
        """
    return rows of the photometry table corresponding to a SN as listed
    in the head table.

    Parameters
    ----------
    row :
    snid :
    """
        if row is not None:
            ptrs = self.headData.iloc[row][['PTROBS_MIN', 'PTROBS_MAX']]
        elif snid is not None:
            ptrs = self.headData.ix[snid][['PTROBS_MIN', 'PTROBS_MAX']]
        else:
            raise ValueError('Both {0} and {1} cannot be None'
                             'simulataneously'.format('snid', 'row'))
        ptrs = ptrs.astype('int').values
        ptrs[0] -= 1
        return ptrs

    def get_SNANA_photometry(self, snid=None, ptrs=None):
        """
    return the photometry table corresponding to a SN with snid (from the
           head table) or the photometry table within the range of row numbers
    indicated by ptrs

        Parameters
        ----------
        snid : 
        ptrs :
    """
        if ptrs is not None:
            assert np.shape(ptrs) == (2,)
        elif snid is not None:
            ptrs = self.get_photrows(snid=snid.strip().lower())
        else:
            raise ValueError('Both {0} and {1} cannot be None'
                             'simulataneously'.format('snid', 'row'))
        lcData = self.phot[1][ptrs[0]: ptrs[1]].byteswap().newbyteorder()
        lcdf = pd.DataFrame(lcData)
        lcdf['zpsys'] = 'ab'
        lcdf['zp'] = 27.5
 
        return LightCurve(lcdf, bandNameDict=self.bandNameDict, ignore_case=True)

class SNChalSims(SNANASims):
    def __init__(self, headData, phot, coerce_inds2int=True,
                 SNANABandNames=lsst_bandNames,
                 registeredBandNames=lsst_bandpassNames):
        """
    Parameters
    ---------
    headFile : string, mandatory
        absolute path to head file of simulation
    photFile : string, mandatory
        absolute path to phot file of simulation
    coerce_inds2int : Bool, optional, defaults to True
        if true, converts SNID from string to int
        SNANABandNames : iterable of strings/characters, optional, defaults to LSST
            characters used to denote the bandpass in the SNANA simulations
        registeredBandNames : iterable of strings, optional, defaults to LSST
            names of the bands registered in SNCosmo
    """
        self.headFile = None
        self.photFile = None
        self.headData = headData
        self.phot = phot
        self.bandNames = SNANABandNames
        self.newbandNames = registeredBandNames
        self.bandNameDict = dict(zip(self.bandNames, self.newbandNames)) 

    @classmethod
    def fromSNChal(cls,snanafileroot, location='./', zipped=True,
                   coerce_inds2int=False, 
                   SNANABandNames=lsst_bandNames,
                   registeredBandNames=lsst_bandpassNames):

        location = os.path.abspath(location)
        filename = os.path.join(location, snanafileroot)

        dictparams=defaultdict(list)
        all_lightcurves = []
        cls._not_available_names={'pixsize': -9.0,'ptrobs_min': -9,'ptrobs_max':-9,'nxpix': 0,'nypix': -9,'mwebv_err': -9, 'redshift_helio': -9,
                'redshift_helio_err': -9,'hostgal_snsep': -9.0,'hostgal_logmass': -9.0,'hostgal_logmass_err': -9.0,
                'hostgal_mag_g': -9.0,'hostgal_mag_r': -9.0,'hostgal_mag_i': -9.0,'hostgal_mag_z': -9.0,
                'hostgal_sb_fluxcal_g': -9.0,'hostgal_sb_fluxcal_r': -9.0,'hostgal_sb_fluxcal_i': -9.0,'hostgal_sb_fluxcal_z': -9.0,
                'peakmjd': -9.0,'search_type': -9}

        if zipped:
            operator = gzip.open
            lcoperator = gzip.GzipFile
        else:
            operator = open
            lcoperator = str

        with operator(filename, 'rb') as f:
            current_lc = list()
            current_head = list()
            for line in f:

                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                if line.startswith('END:'):
                    snid = cls._process_header(current_head, dictparams)
                    cls._process_lc(snid,current_lc, all_lightcurves)
                    current_lc = list()
                    current_head = list()

                elif line.startswith('OBS:'):
                    current_lc.append(line[4:])
                
                else:
                    current_head.append(line)

        params = pd.DataFrame.from_dict(dictparams).set_index('snid')

        lightcurve = pd.DataFrame.from_records(np.array(all_lightcurves, dtype=np.dtype([('snid', np.int), ('mjd', np.float), ('band', '|S16'), ('flux', np.float), ('fluxerr', np.float), ('zp', np.float), ('zpsys','|S8')])))

        return cls(headData=params, phot=lightcurve,
                   coerce_inds2int=coerce_inds2int, SNANABandNames=SNANABandNames,
                   registeredBandNames=registeredBandNames)

    @classmethod
    def _process_header(cls, header_lines, dictparams):
        for line in header_lines:
            k, _, v = line.partition(':')
            if not _:
                continue
            k = k.strip()
            v = v.strip()

            if k == 'SNID':
                dictparams['snid'].append(int(v))
                snid = int(v)
            elif k == 'IAUC':
                dictparams['iauc'].append(v)
            elif k == 'SNTYPE':
                dictparams['sn_type'].append(int(v))
            elif k == 'FILTERS':
                filters = v
            elif k == 'RA':
                dictparams['ra'].append(float(v.split('deg')[0].strip()))
            elif k == 'DECL':
                dictparams['decl'].append(float(v.split('deg')[0].strip()))
            elif k == 'FAKE':
                dictparams['fake'].append(int(v.split('(')[0].strip()))
            elif k == 'MWEBV':
                dictparams['mwebv'].append(float(v.split('MW')[0].strip()))
            elif k == 'REDSHIFT_SPEC':
                dictparams['hostgal_specz'].append(float(v.split('+-')[0].strip()))   
                dictparams['hostgal_specz_err'].append(float(v.split('+-')[1].strip()))
            elif k == 'HOST_GALAXY_GALID':
                dictparams['hostid'].append(int(v))
            elif k == 'HOST_GALAXY_PHOTO-Z':
                dictparams['redshift_final'].append(float(v.split('+-')[0].strip()))   
                dictparams['hostgal_photoz'].append(float(v.split('+-')[0].strip()))   
                dictparams['redshift_final_err'].append(float(v.split('+-')[1].strip()))   
                dictparams['hostgal_photoz_err'].append(float(v.split('+-')[1].strip()))   
            elif k == 'NOBS':
                dictparams['nobs'].append(int(v))

        for k,v in cls._not_available_names.iteritems():
            dictparams[k].append(v)

        return snid

    @staticmethod
    def _process_lc(snid, lc_lines, all_lightcurves):
        for line in lc_lines:
            items = line.strip().split()
            del items[2]
            items.insert(0, snid)
            for k in (1,3,4): items[k]=float(items[k])
            items.append(27.5)
            items.append('ab')
            all_lightcurves.append(tuple(items))

    def get_SNANA_photometry(self, snid=None, ptrs=None):
        """
    return the photometry table corresponding to a SN with snid (from the
           head table) or the photometry table within the range of row numbers
    indicated by ptrs

        Parameters
        ----------
        snid : 
        ptrs :
    """
        if ptrs is not None:
            return self.phot.loc[ptrs[0]:ptrs[1]]
        elif snid is not None:
            return self.phot.query('snid ==' + str(snid)).reset_index(drop=True)
        