"""
Class for holding Light Curve Data
"""
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
import abc
import numpy as np
import pandas as pd
from astropy.table import Table
from .aliases import aliasDictionary


__all__ = ['BaseLightCurve', 'LightCurve']


class BaseLightCurve(with_metaclass(abc.ABCMeta, object)):
    """
    Abstract Base Class for Light Curve Data showing methods that need to be
    implemented.
    """
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
    def coaddedLC(self, coaddTimes=None, format=None, *args, **kwargs):
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

    def __init__(self, lcdf, bandNameDict=None, ignore_case=True):
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
        Example
        -------
        >>> from analyzeSN import LightCurve
        >>> ex_data = sncosmo.load_example_data()
        >>> lc = LightCurve(ex_data.to_pandas()) 
        """
        self.bandNameDict = bandNameDict
        self._lightCurve  = lcdf
        self.ignore_case = ignore_case
        _ = self.lightCurve


    def missingColumns(self, lcdf):

        notFound = self.mandatoryColumns - set(lcdf.columns)
        return notFound

    @staticmethod
    def remap_filters(name, nameDicts, ignore_case=True):
        """
        """
        try:
            if ignore_case:
                return nameDicts[name.lower()]
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

        # Rename columns to standard names if necessary
        aliases = self.columnAliases
        standardNamingDict = aliasDictionary(_lc.columns, aliases)
        if len(standardNamingDict) > 0:
            _lc.rename(columns=standardNamingDict, inplace=True)

        # If all  mandatory columns exist return the light curve, or
        # raise ValueError citing missing columns
        missingColumns = self.missingColumns(_lc)
        if len(missingColumns) > 0:
            raise ValueError('light curve data has missing columns',
                             missingColumns)
        else:
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

    def coaddedLC(self, coaddTimes=None, mjdBefore=None, minmjd=None):
        """
        return a coadded light curve
        """
        if coaddTimes is None:
            return self.lightCurve

        # otherwise perform coadd
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
