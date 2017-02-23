#!/usr/bin/env python
"""
classes and methods to analyze the outputs of light curve characterization
routines, particularly for SALT2 light curves
"""
from __future__ import absolute_import
import numpy as np
import pandas as pd
from copy import deepcopy
from . import cov_utils as cutils

__all__ = ['ResChar']
class ResChar(object):
    """
    A class to hold results of characterizing a light curve fit. Mostly
    this will be constructed from `sncosmo.utils.res` instances
    """

    def __init__(self,
                 vparam_names,
                 param_names,
                 parameters,
                 covariance,
                 errors,
                 samples=None,
                 weights=None,
                 sncosmoModel=None):
        """
        constructor for class

        Parameters
        ----------
        vparam_names : list of strings, mandatory
            model parameters inferred in the characterization
        param_names : list of strings, mandatory
            model parameters (complete list)
        parameters : list of floats, mandatory
            values of all model parameters in same order as param_names
        covariance : `numpy.ndarray`, mandatory
            Covariance of all the varied parameters. Should be a
            symmetric positive definite array of shape(n, n) where
            n = len(vparam_names)
        samples : `numpy.ndarray`, optional, defaults to None
            samples of shape(num, n) where num = number of samples, and
            n = len(vparam_names). May not be independent
        weights : `np.array` , optional defaults to None
            must have length equal to len(samples) if present. For
            mcmc_lc, this is usually used as `np.ones`. For nest_lc
            the weights are used to calculate quantities
        sncsomoModel : `sncosmo.Model`, optional, defaults to None
            model returned from sncsomo estimate
        """

        self.vparam_names = vparam_names
        self.param_names = param_names
        self._parameters = parameters
        self._covariance = covariance
        self.errors = errors
        self._samples = samples
        self.weights = weights
        self.sncosmoModel = sncosmoModel
        self.sample_names = deepcopy(self.vparam_names)

    @classmethod
    def fromSNCosmoRes(cls, SNCosmoRes):
        """
        Instantiate this object from an instance of `sncosmo.utils.res`

        Parameters
        ----------
        SNCosmoRes : instance of `sncosmo.utils.res
        """

        # samples if the method was mcmc/nest_lc but not if max_lc
        # weights makes sense for mcmc methods

        res, model = SNCosmoRes

        samples = None
        if 'samples' in res.keys():
            samples = res['samples']
            weights = np.ones(len(samples))
        else:
            samples = None
            weights = None

        # if method was nest_lc
        if 'weights' in res.keys():
            weights = res['weights']


        return cls(vparam_names=res.vparam_names,
                   param_names=res.param_names,
                   parameters=res.parameters,
                   covariance=res.covariance,
                   errors=res.errors,
                   samples=samples,
                   weights=weights,
                   sncosmoModel=model)

    @property
    def parameters(self):
        """
        return the model parameters as a `pd.Series`
        """
        return pd.Series(self._parameters, index=self.param_names)

    @property
    def vparams(self):
        """
        return the values of the varied parameters as a `pd.Series`
        """
        vparameters = [self._parameters[self.param_names.index(v)]
                       for v in self.vparam_names]
        return pd.Series(vparameters, index=self.vparam_names)

    @property
    def covariance(self):
        """
        return the covariance as a `pd.DataFrame`
        """
        return cutils.covariance(self._covariance, paramNames=self.vparam_names)

    @property
    def samples(self):
        """
        returns the samples if this is the result of mcmc / nest_lc
        characterization
        TODO: Not really thinking about nestlc yet
        """
        if self._samples is None:
            return None
        return pd.DataFrame(self._samples, columns=self.sample_names)

    def salt_covariance_linear(self, x0Truth=None):
        """
        """
        x0 = self.parameters.ix['x0']
        if x0Truth is not None:
            x0 = x0Truth

        factor = - 2.5 /np.log(10)
        # drop other parameters like t0
        cov = self.covariance.copy()
        cov = cutils.subcovariance(covariance=cov,
                                   paramList=['x0', 'x1', 'c'],
                                   array=False)
        covariance = cutils.log_covariance(cov, paramName='x0',
                                    paramValue=x0, factor=factor)

        covariance.rename(columns={'x0':'mB'}, inplace=True)
        covariance['name'] = covariance.columns
        covariance.set_index('name', inplace=True)
        covariance.index.name = None

        return covariance


    def mu_variance_linear(self, alpha=0.14, beta=3.1):
        """
        """
        A = np.array([1.0, alpha, -beta])
        _cov = self.salt_covariance_linear()
        print(_cov, A)
        sc = cutils.subcovariance(_cov, paramList=['mB', 'x1', 'c'],
                                  array=True)
        return cutils.expAVsquare(sc, A)

        return 0.
    def salt_samples(self, alpha=0.14, beta=-3.1, MDelta=0., trimMDelta=True):
        """
        return the samples for SALT2 variables

        Parameters
        ----------
        alpha : float, optional, defaults to 0.14
        beta : float, optional, defaults to 3.1
        MDelta : float, or array-like of size len(samples), optional, defaults
                to 0.
            Offset or environment dependent additive value to obtain the
            correct absolute luminosity/ distance modulus.
        TrimMDelta : remove MDelta and x0 columns. Used for easier plotting
            of posteriors.
        Returns
        -------
        `pd.dataFrame` with extra column names ['mB', 'MDelta', 'mu']

        .. note :
            1. mu = mB + alpha * x1 - beta * c - MB, Eqn. 4 of JLA paper
            While mB is the peak rest frame BessellB magnitude, this is
            related -2.5 *np.log10(x0) through an additive correction. To a
            very good approximation, it is independent of x1, c values due to
            the defining constraints of the SALT2 model.
            2. MB is often taken to be environment dependent. Here we absorb
            the constant part of MB and the additive correction mentioned
            above as MDelta
            3. The constant part of MDelta does not matter if the cosmology
            likelihood is marginalized analytically over a constant additive
            term to the distance modulus.
            4. The values alpha = 0.14 and beta = 3.1 are taken from
            Table. 10 of the JLA paper. This is merely for convenince. For
            actual calculations, these values must be inferred as well.
            JLA Paper : http://arxiv.org/pdf/1401.4064v2.pdf; DeltaM



        """

        samples = self.samples.copy()
        samples['mB'] = -2.5 * np.log10(samples['x0'])
        samples['MDelta'] = MDelta
        samples['mu'] = samples['mB'] + alpha * samples['x1'] \
            + beta * samples['c'] + samples['MDelta']
        del samples['MDelta']
        del samples['x0']
        return samples
