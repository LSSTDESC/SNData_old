"""
Class for holding Light Curve Data
"""
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
import abc
import numpy as np
import pandas as pd
from astropy.table import Table
import sncosmo
from .aliases import aliasDictionary


__all__ = ['BaseLightCurve', 'LightCurve']

class BaseLightCurve(with_metaclass(abc.ABCMeta, object)):
    """
    Abstract Base Class for Light Curve Data showing methods that need to be
    implemented.
    """
    @abc.abstractproperty
    def props(self):
        pass
    @abc.abstractmethod
    def __init__(self):
        pass

    @abc.abstractproperty
    def lightCurve(self):
        """
        `pd.DataFrame` holding the lightCurve information. There can be more
        columns, but the following columns are mandatory:
        ['mjd', 'band', 'flux', 'fluxerr', 'zp', 'zpsys']
        """
        pass

    @abc.abstractmethod
    def snCosmoLC(self, coaddTimes=None):
        pass

    @abc.abstractmethod
    def coaddedLC(self, coaddTimes=None, timeOffset=0., timeStep=1.0, *args, **kwargs):
        pass

    @abc.abstractmethod
    def remap_filters(names, bandNameDict, ignore_case):
        pass

    @abc.abstractmethod
    def missingColumns(self, lcdf):

        notFound = self.mandatoryColumns - set(lcdf.columns)
        return notFound

    @property
    def mandatoryColumns(self):
        """
        A list of mandatory columns in the light curve dataFrame with
        possible aliases in `self.mandatoryColumnAliases`.

        mjd : time
        band : string
        flux : model flux
        """
        reqd = set(['mjd', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'])
        return reqd

    @property
    def columnAliases(self):
        """
        dictionary that maps standard names as keys to a possible set of
        aliases
        """
        aliases = {}
        aliases['zp'] = ['zp']
        aliases['mjd'] = ['time', 'expmjd', 'date']
        aliases['zpsys'] = ['magsys']
        aliases['band'] = ['filter', 'filtername', 'bandname', 'bands', 'flt']
        aliases['flux'] = ['fluxcal']
        aliases['fluxerr'] = ['flux_err', 'flux_errs', 'fluxerror', 'fluxcalerr']
        return aliases

class LightCurve(BaseLightCurve):
    """
    A Class to represent light curve data.  Light curve data is often available
    with different kinds of column names. This class homogenizes them to a set
    of standard names, and allows simple calculations to be based on the same
    variable names 'mjd', 'band', 'flux', 'fluxerr', 'zp', 'zpsys' which denote
    the time of observation, bandpass of observation, the flux and flux
    uncertainty of the observation.

    zp represents the zero point to convert the flux value to the phsyical flux
    using the zero point system zpsys.
    """

    def __init__(self, lcdf, bandNameDict=None, ignore_case=True, propDict=None):
        """
        Instantiate Light Curve class

        Parameters
        ----------
        lcdf : `pd.DataFrame`, mandatory
            light curve information, must contain columns `mjd`, `band`, `flux`,
            `flux_err`, `zp`, `zpsys`
        bandNameDict : dictionary, optional, default to None
            dictionary of the values in the 'band' column or its alias, and
            values that it should be mapped to.
        ignore_case : bool, optional, defaults to True
            ignore the case of the characters in the strings representing
            bandpasses
        propDict : Dictionary, optional, defaults to None
            a dictionary of properties associated with the light curve
        Example
        -------
        >>> from analyzeSN import LightCurve
        >>> ex_data = sncosmo.load_example_data()
        >>> lc = LightCurve(ex_data.to_pandas()) 
        """

        aliases = self.columnAliases
        standardNamingDict = aliasDictionary(lcdf.columns, aliases)
        if len(standardNamingDict) > 0:
            lcdf.rename(columns=standardNamingDict, inplace=True)

        missingColumns = self.missingColumns(lcdf)
        if len(missingColumns) > 0:
            raise ValueError('light curve data has missing columns',
                             missingColumns)

        self.bandNameDict = bandNameDict
        self._lightCurve  = lcdf
        self.ignore_case = ignore_case
        self._propDict = propDict

    @property
    def props(self):
        return self._propDict

    @classmethod
    def fromSALTFormat(cls, fname):
        _lc = sncosmo.read_lc(fname, format='salt2')
        lc = _lc.to_pandas()
        lc.MagSys = 'ab'
        def filtername(x):
            if 'megacam' in x.lower():
                return 'megacam'
            else:
                return x[:-3].lower()
        banddict = dict((key.lower(), filtername(key) + key[-1])
                        for key in lc.Filter.unique())
        return cls(lc, bandNameDict=banddict, ignore_case=True, propDict=_lc.meta)


    def missingColumns(self, lcdf):
        """
        return a set of columns in the light curve dataframe that are missing
        from the mandatory set of columns

        Parameters
        ----------
        lcdf : `pd.dataFrame`
            a light curve represented as a pandas dataframe
        """

        notFound = self.mandatoryColumns - set(lcdf.columns)
        return notFound

    @staticmethod
    def remap_filters(name, nameDicts, ignore_case=True):
        """
        """
        try:
            if ignore_case:
                _nameDicts = dict((key.lower(), value)
                                  for (key, value) in nameDicts.items()) 
                return _nameDicts[name.lower()]
            else:
                return nameDicts[name]
        except:
            raise NotImplementedError('values for old filter {} not implemented',
                                       name)



    @property
    def lightCurve(self):
        """
        The lightcurve in native format
        """
        # light curve
        _lc = self._lightCurve.copy()

        # return the light curve
        _lc.band = _lc.band.apply(lambda x: x.strip())
        if self.bandNameDict is not None:
            _lc.band = _lc.band.apply(lambda x:
                                      self.remap_filters(x, self.bandNameDict,
                                                    self.ignore_case))
        return _lc

    def snCosmoLC(self, coaddTimes=None, mjdBefore=0., minmjd=None):
        lc = self.coaddedLC(coaddTimes=coaddTimes, mjdBefore=mjdBefore,
                            minmjd=minmjd).rename(columns=dict(mjd='time'))
        return Table.from_pandas(lc)

    @staticmethod
    def sanitize_nan(lcs):
        """
        .. note:: These methods are meant to be applied to photometric tables
        as well
        """
        
        lcs = lcs.copy()
        # Stop gap measure to deal with nans
        avg_error  = lcs.fluxerr.mean(skipna=True)
        lcs.fillna(dict(flux=0., fluxerr=avg_error), inplace=True)
        return lcs
    @staticmethod 
    def discretize_time(lcs, timeOffset=0., timeStep=1.0):
        """
        .. note:: These methods are meant to be applied to photometric tables
        as well
        """
        lcs['night'] = (lcs.mjd - timeOffset) // timeStep 
        lcs.night =  lcs.night.astype(np.int)
        return lcs
    @staticmethod
    def add_weightedColumns(lcs, avg_cols=('mjd', 'flux', 'fluxerr', 'zp'),
                         additional_cols=None,
                         copy=False):
        avg_cols = list(tuple(avg_cols))
        if additional_cols is not None:
            avg_cols += list(additional_cols)
        if copy:
            lcs = lcs.copy()
    
        if 'weights' not in lcs.columns:
            if 'fluxerr' not in lcs.columns:
                raise ValueError("Either fluxerr or weights must be a column in the dataFrame")
            lcs['weights'] = 1.0 / lcs['fluxerr']**2
        for col in avg_cols:
            if col != 'fluxerr':
                #lcs['weighted_' + col] = lcs[col] * lcs[col] * lcs['weights'] *lcs['weights']
                #lcs['weights_squared'] = lcs['weights'] *lcs['weights']
                #else:
                lcs['weighted_' + col] = lcs[col] *lcs['weights']
        return lcs
    @staticmethod
    def coadd(preProcessedlcs, include_snid=True,
              cols=('mjd', 'flux', 'fluxerr', 'zp','zpsys'),
              additionalAvgCols=None,
              additionalColsKept=None,
              keepAll=False,
              keepCounts=True):
        """
        additionalAvgCols : list of strings
        
        .. note:: These methods are meant to be applied to photometric tables
        as well
        """
        grouping=['band', 'night']
        if include_snid:
            grouping = ['snid'] + grouping
        
        default_avg_cols = ['mjd', 'flux', 'zp']
        avg_cols = default_avg_cols
        
        if additionalAvgCols is not None:
            avg_cols += additionalAvgCols
            
        default_add_cols = ['zpsys']
        add_cols = default_add_cols
        if additionalColsKept is not None:
            add_cols += additionalColsKept
        
        lcs = preProcessedlcs
        
        
        #lcs = _preprocess(lcs, cols=cols, timeStep=timeStep, timeOffset=timeOffset)
        grouped = lcs.groupby(grouping)
        aggdict = dict(('weighted_' + col, np.sum) for col in avg_cols)
        aggdict['weights'] = np.sum
        if keepCounts:
            lcs['numExpinCoadd'] = lcs.mjd.copy()
            aggdict['numExpinCoadd'] = 'count'
        for col in add_cols:
            aggdict[col] = 'first'
        
        x = grouped.agg(aggdict)
    
        weighted_cols = list(col for col in x.reset_index().columns
                             if (col.startswith('weighted') and col != 'weighted_fluxerr') )
        yy = x.reset_index()[weighted_cols].apply(lambda y: y/x.weights.values, axis=0)
        yy['weighted_fluxerr_coadded'] = 1.0 / np.sqrt(x.reset_index()['weights'])#/x.reset_index()['weighted_fluxerr']**2)#.apply(lambda y: y/x.weights_squared.values)#/x.weights_squared.values)
        yy.rename(columns=dict((col, col.split('_')[1]) for col in yy.columns), inplace=True)
        keepcols = grouping + add_cols
        if keepCounts:
            keepcols += ['numExpinCoadd']
        return x.reset_index()[keepcols].join(yy)

    def coaddedLC(self,
                  coaddTimes=None,
                  minmjd=None,
                  coaddedValues=['mjd', 'flux', 'fluxerr', 'zp'],
		  additionalValues=['zpsys'],
                  mjdBefore=None,
                  sanitize=True):
        """
        """
        # How should we coadd? group observation in steps of coaddTimes and
        # offsets described by minmjd
        if minmjd is None:
            if mjdBefore is None:
                minmjd = 0.
            else:
                minmjd = self.lightCurve.mjd.min() - mjdBefore

	# Does the light curve have `snid` 
	include_snid = 'snid' in self.lightCurve.columns

        # preprocess the light curve for coaddition
        if not sanitize:
            raise NotImplementedError('nan sanitization must be used for coadds\n')
        lc = self.sanitize_nan(self.lightCurve)
        lc = self.discretize_time(lc, timeOffset=minmjd, timeStep=coaddTimes)
        lc = self.add_weightedColumns(lc,
                                      avg_cols=coaddedValues,
                                      additional_cols=None,
                                      copy=True)
	lc = self.coadd(lc, 
			include_snid=include_snid,
		        cols=coaddedValues,
			additionalAvgCols=None,
			additionalColsKept=None,
			keepAll=False,
			keepCounts=True)
	return lc
    def _coaddedLC(self, coaddTimes=None, mjdBefore=None, minmjd=None):
        """
        return a coadded light curve
        """
        if coaddTimes is None:
            return self.lightCurve

        # otherwise perform coadd
        # minmjd provides an offset for calculating discrete times
        if minmjd is None:
            if mjdBefore is None:
                mjdBefore = 0.
            minmjd = self.lightCurve.mjd.min() - mjdBefore

        lc = self.lightCurve.copy()
        lc['discreteTime'] = (lc['mjd'] - minmjd) // coaddTimes
        lc['discreteTime'] = lc.discreteTime.astype(int)


        aggregations = {'mjd': np.mean,
                        'flux': np.mean,
                        'fluxerr': lambda x: np.sqrt(np.sum(x**2))/len(x), 
                        'discreteTime': 'count',
                        'zp': np.mean,
                        'zpsys': 'first'}
        groupedbynightlyfilters = lc.groupby(['discreteTime','band'])
        glc = groupedbynightlyfilters.agg(aggregations)
        glc.reset_index('band', inplace=True)
        glc.rename(columns=dict(discreteTime='numCoadded'), inplace=True) 
        glc['CoaddedSNR'] = glc['flux'] / glc['fluxerr']
        return glc

