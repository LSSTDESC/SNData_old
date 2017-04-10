from __future__ import absolute_import, print_function, division
import fitsio
import pandas as pd
import numpy as np
import os
import gzip
from .lightcurve import LightCurve
from astropy.table import Table, Column


__all__ = ['SNANASims']

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
    
    @classmethod
    def SNChalParser(cls,snanafileroot, zipped=True, location='./',
                          coerce_inds2int=False, 
                          SNANABandNames=lsst_bandNames,
                          registeredBandNames=lsst_bandpassNames):

        '''
        Reads and returns supernovae data and metadata.
        * filename is a string containing the path to the supernovae light curve data
        * survey is a string containing the survey name
        * snid is an integer containing the supernova ID
        * ra is a float containing the RA of the supernova
        * dec is a float containing the Dec of the supernova
        * mwebv is a float describing the dust extinction
        * hostid is an integer containing the host galaxy ID
        * hostz is an array of floats containing the photometric redshift of the galaxy and the error on the measurement
        * spec is an array of floats containing the redshift
        * sim_type is a string containing the supernova type
        * sim_z is a float containing the redshift of the supernova
        * obs is a pandas dataframe containing [observation time, filter, fluxes, fluxe errors]
        '''

        location = os.path.abspath(location)
        filename = os.path.join(location, snanafileroot)

        survey = snid = iauc = ra = decl = fake = mwebv = hostid = hostz = spec = sim_type = sim_z = filters =    None
        if zipped:
            operator = gzip.open
            lcoperator = gzip.GzipFile
        else:
            operator = open
            lcoperator = str

        with operator(filename, 'rb') as f:
            for lineno,line in enumerate(f):
                s = line.split(':')
                if len(s) <= 0:
                    continue
                if s[0] == 'SURVEY':
                    survey = s[1].strip()
                elif s[0] == 'SNID':
                    snid = int(s[1].strip())
                elif s[0] == 'IAUC':
                    iauc = s[1].strip()
                elif s[0] == 'FILTERS':
                    filters = s[1].strip()
                elif s[0] == 'SNTYPE':
                    sn_type = int(s[1].strip())
                elif s[0] == 'RA':
                    ra = float(s[1].split('deg')[0].strip())
                elif s[0] == 'DECL':
                    decl = float(s[1].split('deg')[0].strip())
                elif s[0] == 'MWEBV':
                    mwebv = float(s[1].split('MW')[0].strip())
                elif s[0] == 'FAKE':
                    fake = int(s[1].split('(')[0].strip())
                elif s[0] == 'REDSHIFT_SPEC':
                    spec = float(s[1].split('+-')[0].strip()), float(s[1].split('+-')[1]    .strip())
                elif s[0] == 'HOST_GALAXY_GALID':
                    hostid = int(s[1].strip())
                elif s[0] == 'HOST_GALAXY_PHOTO-Z':
                    hostz = float(s[1].split('+-')[0].strip()), float(s[1].split('+-')[1    ].strip())
                elif s[0] == 'VARLIST':
                    break

        lightcurve = pd.read_table(lcoperator(filename), sep='\s+', header=0, 
                            skiprows= lineno, skipfooter=1,
                            comment='DETECTION', engine='python',
                            usecols=[1,2,3,4,5])
        lightcurve['SNR'] = lightcurve['FLUXCAL']/lightcurve['FLUXCALERR']
        lightcurve['MAG'] = -2.5*np.log10(lightcurve['FLUXCAL'])+27.5
        lightcurve['MAGERR'] = 2.5*np.log10(lightcurve['FLUXCALERR'])
        lightcurve['ZP'] = [27.5 for i in range(len(lightcurve['MAG']))]
        lightcurve['ZPSYS'] = ['AB' for i in range(len(lightcurve['MAG']))]
         
        params = {'survey': survey, 'snid': snid, 'sn_type': sn_type,
                  'ra': ra, 'decl': decl, 'iauc': iauc,
                  'mwebv': mwebv, 'hostid': hostid, 'hostz': hostz, 'spec': spec, 
                  'filename': filename}

         return lightcurve,params

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
