"""
Class for holding Light Curve Data
"""
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
__all__ = ['BasePhotometry', 'Photometry']

import abc
import numpy as np
import pandas as pd
from astropy.table import Table
from .aliases import aliasDictionary



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
