#!/usr/bin/env python
#
# Estimate z0 and z/d simultaneously from single-level micrometeorological
# measurements using the approaches given in Graf et al., 2011.
# 
# Zhan Li, zhanli@gfz-potsdam.de
# Created: Tue Oct  6 17:31:52 CEST 2020
# 
# Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014. Intercomparison of
# Methods for the Simultaneous Estimation of Zero-Plane Displacement and
# Aerodynamic Roughness Length from Single-Level Eddy-Covariance Data.
# Boundary-Layer Meteorol 151, 373–387.
# https://doi.org/10.1007/s10546-013-9905-z

import warnings

import numpy as np
from sklearn.linear_model import LinearRegression

from micro import VON_KARMAN, PHI_M_CONST_U, PHI_M_CONST_S
from micro import _psi_m, _monin_obukhov_z0

from utils import check_is_fitted

class SurfaceAerodynamicFPRE():
    """Flux-Profile-based Regression estimator of surface aerodynamic
    parameters. 

    This estimator takes single-level micrometeorological measurements to
    simultaneously estimate two surface aerodynamic parameters including
    effective/aerodynamic measurement height z and aerodynamic roughness length
    z0 using the regression method based on flux-profile similarity theory
    [1]_. 

    Parameters
    ----------
    solver : str {'univariate', 'bivariate'}
        If 'univariate':
          Use the Eq. (9) and (10) in [1]_ to solve a univariate linear
          regression.
        If 'bivariate':
          Use the Eq. (9) and (11) in [1]_ to solve a bivariate linear
          regression. 

    min_nobs : float
        Minimum number of observations that meet the applicability criteria
        given by Table 1 in [1]_ to carry out the estimates.

    selector : None or callable
        A function to select input data for the estimation. If None, use
        default function that select data according to the Table 1 in [1]_. If
        a user-defined callable, it should be in the following form,
        ``user_function(data, z0max, zmax)``.

    Attributes
    ----------
    N_ : integer
        Number of valid observations that meet the applicability criteria given
        by Table 1 in [1]_.

    Nin_ : integer
        Number of input observations after excluding NaN and Inf.

    References
    ----------
    .. [1] Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014.
    Intercomparison of Methods for the Simultaneous Estimation of Zero-Plane
    Displacement and Aerodynamic Roughness Length from Single-Level
    Eddy-Covariance Data.  Boundary-Layer Meteorol 151, 373–387.
    https://doi.org/10.1007/s10546-013-9905-z
    """

    def __init__(
            self, 
            solver='univariate', 
            min_nobs=10, 
            selector=None):
        self.solver = solver
        self.min_nobs = min_nobs
        self.selector=selector

    def _cleaner(self, data):
        """Clean out NaN and Inf from data.

        Parameters
        ----------
        data : ndarray of shape (n_obs, 3)
            Input ``n_obs`` observations of 3 variables from micrometeorology
            measurements in 3 columns of ``data``: alongwind speed, friction
            velocity, and Monin-Obukohv length, in that order. 
 
        Returns
        -------
        sdata : ndarray of shape (n_clean, 3)
            Selected ``n_clean`` observations of 3 variables after cleaning. 
        """
        data = data.copy()

        sflag = np.logical_or(np.isnan(data), np.isinf(data))
        sflag = np.logical_not(sflag)
        sflag = np.all(sflag, axis=1)
        data = data[sflag, :]
        return data

    def _default_selector(self, data, z0max, zmax):
        """Filter data according to the applicabilit criteria in Table 1 in
        [1]_.

        Parameters
        ----------
        data : ndarray of shape (n_obs, 3)
            Input ``n_obs`` observations of 3 variables from micrometeorology
            measurements in 3 columns of ``data``: alongwind speed, friction
            velocity, and Monin-Obukohv length, in that order. 
 
        z0max : float
            Upper bound to the expected z0 values (aerodynamic surface
            roughness length, e.g. 10% of canopy height for vegetated surface),
            in meters.

        zmax : float
            Upper bound to the expected z values (effective/aerodynamic
            measurement height, in meters.
       
        Returns
        -------
        sdata : ndarray of shape (N_, 3)
            Selected ``N_`` observations of 3 variables after filtering. 
        """
        sflag_arr = [ 
                1 / data[:, 2] < -0.103 / zmax, 
                1 / data[:, 2] < -0.084 / z0max, 
                1 / data[:, 2] >  0.037 / z0max, 
                1 / data[:, 2] >      1 / zmax, 
                ]
        sflag = np.logical_not(np.any(np.vstack(sflag_arr).T, axis=1))
        data = data[sflag, :]
        return data
 
    def _filter(self, data, z0max, zmax):
        """Filter data according to the applicabilit criteria in Table 1 in
        [1]_.

        Parameters
        ----------
        data : ndarray of shape (n_obs, 3)
            Input ``n_obs`` observations of 3 variables from micrometeorology
            measurements in 3 columns of ``data``: alongwind speed, friction
            velocity, and Monin-Obukohv length, in that order. 
 
        z0max : float
            Upper bound to the expected z0 values (aerodynamic surface
            roughness length, e.g. 10% of canopy height for vegetated surface),
            in meters.

        zmax : float
            Upper bound to the expected z values (effective/aerodynamic
            measurement height, in meters.
       
        Returns
        -------
        sdata : ndarray of shape (N_, 3)
            Selected ``N_`` observations of 3 variables after filtering. 
        """
        data = self._cleaner(data)
        self.Nin_ = data.shape[0]

        selector2use = self._default_selector
        if self.selector is not None and callable(self.selector):
            selector2use = self.selector
        data = selector2use(data, z0max, zmax)
        self.N_ = data.shape[0]
        return data
   
    def _calc_xy(self, data):
        """Calculate response variable y and explanatory variable X for the
        linear regression approach. 

        Returns
        -------
        X : ndarray of shape (N_, 1 or 2)
            If the solver is 'univariate', X has the shape of (N_, 1). Refer to
            Eq.  (10) in [1]_. If the solver is 'bivariate', X has the shape
            (N_, 2). Refer to Eq. (11) in [1]_, with the 1st column of X being
            :math:`u_*/\kappa`, and the 2nd column of X being :math:`(u_*
            \\beta)/(\kappa L)`

        y : ndarray of shape (N_, 1)
            Refer to Eq. (10) or (11) in [1]_.

        References
        ----------
        .. [1] Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014.
        Intercomparison of Methods for the Simultaneous Estimation of
        Zero-Plane Displacement and Aerodynamic Roughness Length from
        Single-Level Eddy-Covariance Data.  Boundary-Layer Meteorol 151,
        373–387.  https://doi.org/10.1007/s10546-013-9905-z
        """
        if self.solver == 'univariate':
            X = PHI_M_CONST_S / data[:, 2] 
            y = data[:, 0] * VON_KARMAN / data[:, 1] 
            X, y = X[:, np.newaxis], y[:, np.newaxis]
        elif self.solver == 'bivariate':
            X = [
                    data[:, 1] / VON_KARMAN, 
                    data[:, 1]*PHI_M_CONST_S / (VON_KARMAN*data[:, 2])]
            X = np.vstack(X).T
            y = data[:, 0][:, np.newaxis]
        else:
            raise ValueError('Unrecognized solver={0:s}'.format(self.solver))
        return X, y

    def fit_transform(self, data, z0max, zmax):
        """Build the estimator from single-level micrometeorolgoical data and
        apply it to estimate surface aerodynamic parameters.

        Parameters
        ----------
        data : ndarray of shape (n_obs, 3)
            Input ``n_obs`` observations of 3 variables from micrometeorology
            measurements in 3 columns of ``data``: alongwind speed, friction
            velocity, and Monin-Obukohv length, in that order. 

        z0max : float
            Upper bound to the expected z0 values (aerodynamic surface
            roughness length, e.g. 10% of canopy height for vegetated surface),
            in meters.

        zmax : float
            Upper bound to the expected z values (effective/aerodynamic
            measurement height, in meters.

        Returns
        -------
        z : float 
            Estimated z (effective/aerodynamic measurement height). 

        z0 : float 
            Estimated z0 (surface aerodynamic roughness length). 
        """
        self.z0max_ = z0max
        self.zmax_ = zmax

        data = self._filter(data, z0max, zmax)

        if self.N_ < self.min_nobs: 
            warnings.warn(
                    'too few data points left after validation.'
                    ' return nan')
            z, z0 = np.nan, np.nan
            self.reg_ = None
        else:
            X, y = self._calc_xy(data)
            if self.solver == 'univariate':
                # 1st regression model (univariate linear regression)
                # Eq. 10 in Graf et al. 2014
                reg = LinearRegression().fit(X, y)
                z = reg.coef_[0, 0]
                z0 = z / np.exp(reg.intercept_[0])
            elif self.solver == 'bivariate':
                # 2nd regression model (bivariate linear regression)
                # Eq. 11 in Graf et al. 2014
                reg = LinearRegression(fit_intercept=False).fit(X, y)
                z = reg.coef_[0, 1]
                z0 = z / np.exp(reg.coef_[0, 0])
            else:
                raise ValueError(
                        'Unrecognized solver={0:s}'.format(self.solver))
            self.reg_ = reg
        return z, z0

class SurfaceAerodynamicFPIT():
    """Flux-Profile-based Iterative estimator of surface aerodynamic
    parameters. 

    This estimator takes single-level micrometeorological measurements to
    simultaneously estimate two surface aerodynamic parameters including
    effective/aerodynamic measurement height z and aerodynamic roughness length
    z0 using the iterative method based on flux-profile similarity theory
    [1]_, [2]_. 

    Parameters
    ----------
    solver : str {'sigma-s', 'sigma-s-approx'}
        If 'sigma-s':
          Use the Eq. (6) and (7) in [1]_ to solve the minimization of the
          exact version of :math:`\Sigma_{S}`.
        If 'sigma-s-approx':
          Use the Eq. (6) and (8) in [1]_ to solve the minimization of the
          approximate version of :math:`\Sigma_{S}`.

    min_nobs : float
        Minimum number of observations that meet the applicability criteria
        given by Table 1 in [1]_ to carry out the estimates.

    selector : None or callable
        A function to select input data for the estimation. If None, use
        default function that select data according to the Table 1 in [1]_. If
        a user-defined callable, it should be in the following form,
        ``user_function(data, z0max, zmax)``.

    Attributes
    ----------
    N_ : integer
        Number of valid observations that meet the applicability criteria given
        by Table 1 in [1]_.

    Nin_ : integer
        Number of input observations after excluding NaN and Inf.

    References
    ----------
    .. [1] Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014.
    Intercomparison of Methods for the Simultaneous Estimation of Zero-Plane
    Displacement and Aerodynamic Roughness Length from Single-Level
    Eddy-Covariance Data.  Boundary-Layer Meteorol 151, 373–387.
    https://doi.org/10.1007/s10546-013-9905-z

    .. [2] Martano, P., 2000. Estimation of Surface Roughness Length and
    Displacement Height from Single-Level Sonic Anemometer Data. J. Appl.
    Meteor. 39, 708–715.
    https://doi.org/10.1175/1520-0450(2000)039<0708:EOSRLA>2.0.CO;2

    """
    def __init__(
            self, 
            solver='sigma-s', 
            min_nobs=10, 
            selector=None):
        self.solver = solver
        self.min_nobs = min_nobs
        self.selector = selector

    def _cleaner(self, data):
        data = data.copy()

        sflag = np.logical_or(np.isnan(data), np.isinf(data))
        sflag = np.logical_not(sflag)
        sflag = np.all(sflag, axis=1)
        data = data[sflag, :]
        return data

    def _default_selector(self, data, z0max, zmax):
        sflag_arr = [ 
                data[:, 0] < 1.5, 
                1 / data[:, 2] < -0.084 / z0max, 
                1 / data[:, 2] >  0.037 / z0max, 
                1 / data[:, 2] >      1 / zmax, 
                ]
        sflag = np.logical_not(np.any(np.vstack(sflag_arr).T, axis=1))
        data = data[sflag, :]
        return data

    def _filter(self, data, z0max, zmax):
        """Filter data according to the applicabilit criteria in Table 1 in
        [1]_.

        Parameters
        ----------
        data : ndarray of shape (n_obs, 3)
            Input ``n_obs`` observations of 3 variables from micrometeorology
            measurements in 3 columns of ``data``: alongwind speed, friction
            velocity, and Monin-Obukohv length, in that order. 

        z0max : float
            Upper bound to the expected z0 values (aerodynamic surface
            roughness length, e.g. 10% of canopy height for vegetated surface),
            in meters.

        zmax : float
            Upper bound to the expected z values (effective/aerodynamic
            measurement height, in meters.

        Returns
        -------
        sdata : ndarray of shape (N_, 3)
            Selected ``N_`` observations of 3 variables after filtering. 
        """
        data = self._cleaner(data)
        self.Nin_ = data.shape[0]

        selector2use = self._default_selector
        if self.selector is not None and callable(self.selector):
            selector2use = self.selector
        data = selector2use(data, z0max, zmax)
        self.N_ = data.shape[0]
        return data       

    def _eval_objective_func(self, data, zv):
        zv_mat, lm_mat = np.meshgrid(zv, data[:, 2])
        _, u_mat = np.meshgrid(zv, data[:, 0])
        _, ustar_mat = np.meshgrid(zv, data[:, 1])
        z0v = _monin_obukhov_z0(zv_mat, u_mat, ustar_mat, lm_mat)

        z0_values = np.nanmean(z0v, axis=0)

        if self.solver == 'sigma-s':
            psim = _psi_m(zv_mat, lm_mat)
            sv = VON_KARMAN * u_mat / ustar_mat - psim
            objective_values = np.nanstd(sv, axis=0)
        elif self.solver == 'sigma-s-approx':
            objective_values = np.nanstd(z0v, axis=0)/ z0_values
        else:
            raise ValueError(
                    'Unrecognized solver={0:s}'.format(self.solver))
        return objective_values, z0_values

    def fit_transform(self, data, z0max, zmax, zv):
        """Build the estimator from single-level micrometeorolgoical data and
        apply it to estimate surface aerodynamic parameters.

        Parameters
        ----------
        data : ndarray of shape (n_obs, 3)
            Input ``n_obs`` observations of 3 variables from micrometeorology
            measurements in 3 columns of ``data``: alongwind speed, friction
            velocity, and Monin-Obukohv length, in that order. 

        z0max : float
            Upper bound to the expected z0 values (aerodynamic surface
            roughness length, e.g. 10% of canopy height for vegetated surface),
            in meters.

        zmax : float
            Upper bound to the expected z values (effective/aerodynamic
            measurement height, in meters.

        zv : ndarray of shape (n_possible, )
            List of ``n_possible`` possible z values in meters for numerical
            search of optimal z and z0.
       
        Returns
        -------
        z : float 
            Estimated z (effective/aerodynamic measurement height). 

        z0 : float 
            Estimated z0 (surface aerodynamic roughness length). 
        """
        self.z0max_ = z0max
        self.zmax_ = zmax

        data = self._filter(data, z0max, zmax)
        if self.N_ < self.min_nobs: 
            warnings.warn(
                    'too few data points left after validation.'
                    ' return nan')
            z, z0 = np.nan, np.nan
        else:
            objective_values, z0_values = self._eval_objective_func(data, zv)
            ix = np.nanargmin(objective_values)
    
            if ix == 0 or ix == len(zv)-1:
                z, z0 = np.nan, np.nan
            else:
                z, z0 = zv[ix], z0_values[ix]

        return z, z0

class SurfaceAerodynamicFVIT():
    """Flux-Variance-based Iterative estimator of surface aerodynamic
    parameters. 

    This estimator takes single-level micrometeorological measurements to
    simultaneously estimate two surface aerodynamic parameters including
    effective/aerodynamic measurement height z and aerodynamic roughness length
    z0 using the iterative method based on flux-variance similarity theory
    [1]_, [2]_. 

    Parameters
    ----------
    solver : str {'sigma-w', 'sigma-t', 'sigma-w-reg'}
        If 'sigma-w':
          Estimate z (effective/aerodynamic height) using the Eq. (12) in [1]_
          to solve the minimization of the objective function concerning
          vertical wind velocity (w). Then use Eq. (6) in [1]_ to estimate z0
          (aerodynamic roughness lenght) based on log wind profile and
          Monin-Obukhov similarity theory.
        If 'sigma-t':
          Estimate z (effective/aerodynamic height) using the Eq. (13) in [1]_
          to solve the minimization of the objective function concerning
          virtual/sonic temperature (w). Then use Eq. (6) in [1]_ to estimate
          z0 (aerodynamic roughness lenght) based on log wind profile and
          Monin-Obukhov similarity theory.
        If 'sigma-w-reg':
          Estimate z (effective/aerodynamic height) using the Eq. (12) in [1]_
          to solve the minimization of the objective function concerning
          vertical wind velocity (w). Then use Eq. (14) in [1]_ in a linear
          regression to estimate z0 (aerodynamic roughness lenght).

    min_nobs : float
        Minimum number of observations that meet the applicability criteria
        given by Table 1 in [1]_ to carry out the estimates.

    selector : None or callable
        A function to select input data for the estimation. If None, use
        default function that select data according to the Table 1 in [1]_. If
        a user-defined callable, it should be in the following form,
        ``user_function(data, z0max, zmax)``.

    Attributes
    ----------
    N_ : integer
        Number of valid observations that meet the applicability criteria given
        by Table 1 in [1]_.

    Nin_ : integer
        Number of input observations after excluding NaN and Inf.

    References
    ----------
    .. [1] Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014.
    Intercomparison of Methods for the Simultaneous Estimation of Zero-Plane
    Displacement and Aerodynamic Roughness Length from Single-Level
    Eddy-Covariance Data.  Boundary-Layer Meteorol 151, 373–387.
    https://doi.org/10.1007/s10546-013-9905-z

    .. [2] Panofsky, H.A., 1984. Vertical variation of roughness length at the
    Boulder Atmospheric Observatory. Boundary-Layer Meteorol 28, 305–308.
    https://doi.org/10.1007/BF00121309
    """
    def __init__(
            self, 
            solver='sigma-w', 
            min_nobs=10, 
            selector=None, 
            selector_reg=None):
        self.solver = solver
        self.min_nobs = min_nobs
        self.selector = selector
        self.selector_reg = selector_reg
        # estimates for the universal constants, as used for this method by
        # Toda & Sugita 2003 and Panofsky & Dutton 1984.
        self._C1 = 1.3
        self._C2 = 2.0
        self._C3 = 0.99
        self._C4 = 0.06;

    def _cleaner(self, data):
        data = data.copy()

        sflag = np.logical_or(np.isnan(data), np.isinf(data))
        sflag = np.logical_not(sflag)
        sflag = np.all(sflag, axis=1)
        data = data[sflag, :]
        return data

    def _default_selector(self, data, z0max, zmax):
        sflag_arr = [
                # remove positive L and negative near-zero L, we need
                # moderately unstable condition.
                zmax / data[:, 2] > 0, 
                zmax / data[:, 2] < -99999, 
                # exclude near-zero denominator in the left-hand side of Eq.
                # (12) and (13) in Graf et al., 2014. 
                data[:, 1] <  0.05, 
                np.abs(data[:, 5]) < 0.3, 
                ]

        sflag = np.logical_not(np.any(np.vstack(sflag_arr).T, axis=1))
        data = data[sflag, :]
        return data

    def _default_selector_reg(self, data, z0max, zmax):
        sflag_arr = [
                1 / data[:, 2] < -0.4, 
                1 / data[:, 2] >  0.4, 
                data[:, 0] < 1.5, 
                ]

        sflag = np.logical_not(np.any(np.vstack(sflag_arr).T, axis=1))
        data = data[sflag, :]
        return data

    def _filter(self, data, z0max, zmax, reg=False):
        """Filter data according to the applicabilit criteria in Table 1 in
        [1]_.

        Parameters
        ----------
        data : ndarray of shape (n_obs, 6)
            Input ``n_obs`` observations of 3 variables from micrometeorology
            measurements in 6 columns of ``data``: alongwind speed, friction
            velocity, Monin-Obukohv length, standard deviation of vertical
            wind, standard deviation of sonic temperature, and friction
            temperature, in that order. 

        z0max : float
            Upper bound to the expected z0 values (aerodynamic surface
            roughness length, e.g. 10% of canopy height for vegetated surface),
            in meters.

        zmax : float
            Upper bound to the expected z values (effective/aerodynamic
            measurement height, in meters.

        reg : boolean
            Whether data filtering is to use the linear regression of Eq. (14)
            in [1]_ to esitmate z0. If so, more criteria are applied to filter
            input data. 

        Returns
        -------
        sdata : ndarray of shape (N_, 3)
            Selected ``N_`` observations of 3 variables after filtering. 

        References
        ----------
        .. [1] Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014.
        Intercomparison of Methods for the Simultaneous Estimation of
        Zero-Plane Displacement and Aerodynamic Roughness Length from
        Single-Level Eddy-Covariance Data.  Boundary-Layer Meteorol 151,
        373–387.  https://doi.org/10.1007/s10546-013-9905-z
        """
        data = self._cleaner(data)
        self.Nin_ = data.shape[0]

        if reg:
            selector2use = self._default_selector_reg
            if self.selector_reg is not None and callable(self.selector_reg):
                selector2use = self.selector_reg
        else:
            selector2use = self._default_selector
            if self.selector is not None and callable(self.selector):
                selector2use = self.selector
        
        data = selector2use(data, z0max, zmax)
        self.N_ = data.shape[0]
        return data       

    def _eval_objective_func(self, data, zv):
        C1, C2, C3, C4 = self._C1, self._C2, self._C3, self._C4

        zv_mat, lm_mat = np.meshgrid(zv, data[:, 2])
        _, stdw_mat = np.meshgrid(zv, data[:, 3])
        _, ustar_mat = np.meshgrid(zv, data[:, 1])
        _, stdt_mat = np.meshgrid(zv, data[:, 4])
        _, tstar_mat = np.meshgrid(zv, data[:, 5])

        _, u_mat = np.meshgrid(zv, data[:, 0])
        z0v = _monin_obukhov_z0(zv_mat, u_mat, ustar_mat, lm_mat)
        z0_values = np.nanmean(z0v, axis=0)

        if self.solver == 'sigma-w' or self.solver == 'sigma-w-reg':
            lhs = stdw_mat / ustar_mat
            rhs = C1 * (1 - C2 * zv_mat / (lm_mat))**(1./3)
            objective_values = np.nanmean((rhs - lhs)**2, axis=0)
        elif self.solver == 'sigma-t':
            lhs = stdt_mat / np.abs(tstar_mat)
            rhs = C3 * (C4 - zv_mat / (lm_mat))**(-1./3)
            objective_values = np.nanmean(((rhs - lhs) / rhs)**2, axis=0)
        else:
            raise ValueError(
                    'Unrecognized solver={0:s}'.format(self.solver))
        
        return objective_values, z0_values

    def _calc_xy(self, data):
        """Calculate response variable y and explanatory variable X for the
        linear regression approach given by Eq. (14) in [1]_ to estimate z0
        (aerodynamic roughness length). Only valid if the ``solver`` ==
        'sigma-w-reg'.

        Returns
        -------
        X : ndarray of shape (N_, 1)
            Refer to Eq. (14) in [1]_.

        y : ndarray of shape (N_, 1)
            Refer to Eq. (14) in [1]_.

        References
        ----------
        .. [1] Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014.
        Intercomparison of Methods for the Simultaneous Estimation of
        Zero-Plane Displacement and Aerodynamic Roughness Length from
        Single-Level Eddy-Covariance Data.  Boundary-Layer Meteorol 151,
        373–387.  https://doi.org/10.1007/s10546-013-9905-z
        """
        if self.solver == 'sigma-w-reg':
            X = data[:, 0][:, np.newaxis]
            y = data[:, 3][:, np.newaxis]
        else:
            raise ValueError('Invalid for solver={0:s}'.format(self.solver))
        return X, y

    def fit_transform(self, data, z0max, zmax, zv):
        """Build the estimator from single-level micrometeorolgoical data and
        apply it to estimate surface aerodynamic parameters.

        Parameters
        ----------
        data : ndarray of shape (n_obs, 6)
            Input ``n_obs`` observations of 3 variables from micrometeorology
            measurements in 6 columns of ``data``: alongwind speed, friction
            velocity, Monin-Obukohv length, standard deviation of vertical
            wind, standard deviation of sonic temperature, and friction
            temperature, in that order. 

        z0max : float
            Upper bound to the expected z0 values (aerodynamic surface
            roughness length, e.g. 10% of canopy height for vegetated surface),
            in meters.

        zmax : float
            Upper bound to the expected z values (effective/aerodynamic
            measurement height, in meters.

        zv : ndarray of shape (n_possible, )
            List of ``n_possible`` possible z values in meters for numerical
            search of optimal z and z0.
       
        Returns
        -------
        z : float 
            Estimated z (effective/aerodynamic measurement height). 

        z0 : float 
            Estimated z0 (surface aerodynamic roughness length). 
        """
        self.z0max_ = z0max
        self.zmax_ = zmax

        sdata = self._filter(data, z0max, zmax)
        if self.N_ < self.min_nobs: 
            warnings.warn(
                    'too few data points left after validation.'
                    ' return nan')
            z, z0 = np.nan, np.nan
        else:
            objective_values, z0_values = self._eval_objective_func(sdata, zv)
            ix = np.nanargmin(objective_values)
    
            if ix == 0 or ix == len(zv)-1:
                z, z0 = np.nan, np.nan
            else:
                z, z0 = zv[ix], z0_values[ix]

        if self.solver == 'sigma-w-reg':
            sdata = self._filter(data, z0max, zmax, reg=True)
            if self.N_ < 1:
                warnings.warn(
                        'No observations meet applicability criteria'
                        ' for {0:s}'.format(self.solver))
                z0 = np.nan
            else:
                C1 = self._C1
                X, y = self._calc_xy(sdata)
                reg = LinearRegression(fit_intercept=False).fit(X, y)
                slope = reg.coef_[0, 0]
                z0 = z / np.exp(VON_KARMAN * C1 / slope)

        return z, z0
