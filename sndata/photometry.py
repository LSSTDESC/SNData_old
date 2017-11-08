"""
Class for holding Photometry Tables
"""
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass

__all__ = ['PhotTables']

import abc
import numpy as np
import pandas as pd
from astropy.table import Table
from .aliases import aliasDictionary
from .lightcurve import LightCurve

class PhotTables(object):
    """
    Data Structures for representing photometry tables or collections of light
    curves. The minimal requirement is that this has all the columns of a
    supernova light curve, but also an index to identify the SN.
    """
    def __init__(self, df, sanitize_nans=True):
        """
        Instantiate the photometry table
        Parameters
        ----------
        df: `pd.DataFrame`
            table of photometry containing at least the following columns
            ('snid', 'mjd', 'flux', 'fluxerr', 'zp', 'zpsys') or aliases
            thereof. Currently, the aliases are standardized using `LightCurve`
            functionality
        sanitize_nans: `Bool`, defaults to True
            if `True`, `nans` in the table are replaced using
            `LightCurve.sanitize_nan`
        """
        # standardized names
        self._lcs = LightCurve(df)
        self._lcs = self._lcs.lightCurve
        self.nan_sanitized = sanitize_nans
        lcs = self._lcs.copy()
        if self.nan_sanitized:
            lcs = LightCurve.sanitize_nan(lcs)
        self.lcs = lcs


    @property
    def mandatoryColumns(self):
        """
        A list of mandatory columns in photTables dataFrame with
        possible aliases in `self.mandatoryColumnAliases`.

        We need to have the mandatory columns for a light curve and `snid`
        """
        new = set(['snid'])

        # columns from a light curve
        reqd = set(['mjd', 'band', 'flux', 'fluxerr', 'zp', 'zpsys']).union(new)
        return reqd

    def coaddedTable(self, timeOffset=0., timeStep=1.0, 
                     avg_cols=('mjd', 'flux', 'fluxerr', 'zp'),
                     additionalAvgCols=None,
                     additionalColsKept=('tileID', 'fieldID', 'zpsys'),
                     additionalAggFuncs=('first', 'first', 'first'),
                     prepend_colNames='coadd_'):
        """
        returns the photometry table coadded over a timeStep of `timeStep` and
        offset of `timeOffset`. 

        Parameters
        ----------
        timeOffset : float, unit of days, defaults to 0.
            offset used in discretization of time for coaddition
        timeStep : float, units of days, defaults to 1.0 
            time period over which ovservations are coadded
        avg_cols : tuple of column names, defaults to light curve defaults
            columns in the photometry table over which inverse variance
            weighted averages are desired. While `fluxerr` does not satisfy
            the above definition it is included too.
        additionalAvgCols : tuple, deprecated, will probably be removed
            meant to be additional column names which will also be averaged with
            the same weights
        additionalColsKept : tuple of strings
            Columns in the photometry table that are desired to be included in
            the coadded result but do not require weighted averages.
        additionalAggFuncs : aggregate function or tuple thereof
            single aggregate function common to all the
        prepend_colNames : string, defaults to 'coadd_'
            string to be prepended to column names aside from `snid`, if left
            as `None`, then no prepending will happen
        """
        include_snid = 'snid' in self.lcs.columns
        if not include_snid:
            raise ValueError('the photTable does not include a column for SNID\n')

        lcs = self.lcs.copy()
        lcs = LightCurve.discretize_time(lcs, timeOffset=timeOffset, timeStep=timeStep)
        lcs = LightCurve.add_weightedColumns(lcs,
                                             avg_cols=avg_cols,
                                             additional_cols=additionalAvgCols,
                                             copy=True)
        weightedcols = list(avg_cols)
        if additionalAvgCols is not None:
            weightedcols += list(additionalAvgCols)
 
        lcs = LightCurve.coaddpreprocessed(lcs, include_snid=include_snid,
                                           cols=weightedcols,
                                           additionalColsKept=additionalColsKept,
                                           additionalAggFuncs='first',
                                           keepAll=False,
                                           keepCounts=True)
        if prepend_colNames is not None:
            coldict = dict((col, prepend_colNames + col) for col in lcs.columns
                           if col not in  ('snid', 'band'))
            lcs.rename(columns=coldict, inplace=True)

        return lcs

    def summary(self,
                coadd=True,
                coaddTimeStep=1.0,
                coaddTimeOffset=0.0,
                paramsdf=None):
        """
        return a summarized multiband light curve with reasonable autogenerated
        summary names

        Parameters
        ----------
        coadd : Bool, defaults to True
            Calculate the same summary parameters if the coadds are included
        coaddTimeStep : float, units of days, default to 1.0
            Time steo for coaddition used if `coadd` parameter is True
        coaddTimeOffset : float, defaults to 0.
            Time offset used to define the grouping of exposures for coaddition
        paramsdf: `pd.DataFrame`, default to None
            contains truth and other metadata about the astrophysical object
            involved. if not `None`, it is joined to the summary  
        """
        summary = LightCurve.summarize(self.lcs, paramsdf=paramsdf)
        if coadd:
            tmp = self.coaddedTable(timeStep=coaddTimeStep,
                                    timeOffset=coaddTimeOffset,
                                    prepend_colNames='')

            nightlySummary = LightCurve.summarize(tmp,
                                                  summary_prefix='coadd_')
            summary = summary.join(nightlySummary)

        # Make sure that some dtypes are converted into ints 
        intcols = list(col for col in summary.columns if 'obs' in col.lower())
        summary[intcols] = summary[intcols].fillna(0)
        if 'tileID' in summary.columns:
            intcols += ['tileID']
        try:
            summary[intcols] = summary[intcols].astype(np.int)
        except:
            pass

        return summary


class BasePhotometry(with_metaclass(abc.ABCMeta, object)):
    def __init__(self, lcs, maxObsHistID, singleLCProps):
        """
        Instantiate class with a collection of or single light curve.

        Parameters
        ----------
        lcs : `pd.DataFrame` having mandatory columns needed to instantiate a
        LightCurve class. In addition, it must have the column `snid` and may
        or may not have the column `PPID` (but might be the index).
        """
        pass

    @abc.abstractmethod
    def append(self, lcs, recalculatePPID=True):
        """
        given a dataframe representing a single or multiple
        `LightCurve.lightCurve` objects, add it as new rows
        to the pandas dataFrame `lcs`.

        Parameters
        ----------
        lcs : dataframe representing single or multiple lightcurves

        recalculatePPID : Bool, defaults to True
            if True, recalculates `PPID` for the dataframe.
        """
        pass

    def singleLCProperties(self):
        """
        immutable sequence of properties that are columns in the light curve
        corresponding to which single light curve properties will be
        constructed
        """
        return self._singleLCProperties

    @property
    def singleBandLC(self):
        """
        return the `pd.dataFrame.groupby` object with each group representing
        a single band light curve of a supernova.
        """
        return self.lightCurve.groupby([['snid', 'band']])

    @staticmethod
    def statTable(PropTuple, callableTuple,
                  grouped=None,
                  dataframe=None,
                  groupTuple=None):
        """
        Parameters
        ----------
        PropTuple :
        callableTuple :
        grouped :
        dataframe :
        groupTuple :
        """
        if grouped is None:
            grouped = dataframe.groupby(list(groupTuple))
        callable_strings = list(s.__name__ for s in callableTuple)
        xx = grouped.agg(dict(zip(PropTuple, callableTuple)))
        xx.columns = ('_'.join(name) for name in zip(callable_strings, xx.columns))
        xx = xx.reset_index().pivot_table(index='snid', columns='band')
        xx.columns = ['_'.join(col).strip() for col in xx.columns.values]
        return xx 

    @abc.abstractmethod
    def pair_method(self, obsHistID, snid, maxObsHistID):
        """
        Combine the obsHistID and snid to form a single index.
        """
        pass

class Photometry(BasePhotometry):
    def __init__(self,
                 lcs,
                 maxObsHistID=10000000,
                 singleLCProps=None):
        """
        """
        self._lcs = lcs.reset_index()
        if singleLCProps is None:
            self._singleLCProperties = ('ModelFlux', 'SNR')
        else:
            self._singleLCProperties = singleLCProps
        self.maxObsHistID = maxObsHistID


    def inverse_pair(self, photID, maxObsHistID=10000000):
        """
        """
        snid = np.floor_divide(photID, maxObsHistID)
        obsHistID = np.mod(photID, maxObsHistID)
        return snid, obsHistID
        
    @staticmethod
    def pair_method(obsHistID, snid, maxObsHistID):
        """
        Combine the obsHistID and snid to form a single index.

        Parameters
        ----------
        snid : int, `np.ndarray`
            unique ID representing each object
        obsHistID : int, `np.ndarray`
            unique ID representing a pointing
        maxObsHistID : int (scalar)
            max value of obsHistID
        """
        return snid * maxObsHistID + obsHistID

    def singleLCProperties(self):
        """
        immutable sequence of properties that are columns in the light curve
        corresponding to which single light curve properties will be
        constructed
        """
        return self._singleLCProperties

    @property
    def singleBandLC(self):
        """
        return the `pd.dataFrame.groupby` object with each group representing
        a single band light curve of a supernova.
        """
        return self.lightCurve.groupby([['snid', 'band']])

    @staticmethod
    def statTable(PropTuple, callableTuple,
                  grouped=None,
                  dataframe=None,
                  groupTuple=None):
        """
        Parameters
        ----------
        PropTuple :
        callableTuple :
        grouped :
        dataframe :
        groupTuple :
        """
        if grouped is None:
            grouped = dataframe.groupby(list(groupTuple))
        callable_strings = list(s.__name__ for s in callableTuple)
        xx = grouped.agg(dict(zip(PropTuple, callableTuple)))
        xx.columns = ('_'.join(name) for name in zip(callable_strings, xx.columns))
        xx = xx.reset_index().pivot_table(index='snid', columns='band')
        xx.columns = ['_'.join(col).strip() for col in xx.columns.values]
        return xx 
