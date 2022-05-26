% Ensemble outputs of z0 and d from the six approaches for d (z_arr, in the
% script) and the seven approaches for z0 (z0_arr in the script) that are
% described in (Graf et al., 2014).
% 
% Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014. Intercomparison of
% Methods for the Simultaneous Estimation of Zero-Plane Displacement and
% Aerodynamic Roughness Length from Single-Level Eddy-Covariance Data.
% Boundary-Layer Meteorol 151, 373–387.
% https://doi.org/10.1007/s10546-013-9905-z
%
% Zhan Li, zhanli@gfz-potsdam.de Created: Wed Sep  9 15:14:47 CEST 2020

function [z_arr, z0_arr, nobs_arr] = ...
    z0d_ensemble(data, ...
                 z0max, ...
                 zmax, ...
                 zv, ...
                 wd_binsize, ...
                 half_wd_win, ...
                 min_nobs)
    %{
    Ensemble outputs of z0 and z from the approaches described in (Graf et
    al., 2014)

    Parameters
    ----------
    data : matrix of shape (nobs, 6)
        Input data matrix, one observation per row and one variable per column.
        The columns provide 6 variables in the following order: wind_speed,
        ustar, monin_obukhov_length, std_w, std_sonic_temperature, Tstar,
        wind_direction

    z0max : float
        Upper bound to the expected z0 values (e.g. 0.1h), in meters.

    zmax : float
        Upper bound to the expected z values, in meters.

    zv : vector of shape (nz,)
        List of possible z values in meters for numerical search of optimal z
        and z0 in the FP-IT approaches and FV-IT approaches.  

    wd_binsize : float
        Bin size of wind direction, in degrees. The z and z0 will be estimated
        per each bin to account for possible spatial variation along wind
        direction. If wd_binsize <=0, no data binning based on wind direction.
        All the input observations will be used in the estimation together. 

    half_wd_win : float
        Half window size of wind directions to increase observation counts per
        each bin of wind direction and also smooth estimates along wind
        direction. Observations in each bin of wind direction plus and mius
        this half window size will be used in the estimation of z and z0 for
        this bin. 

    min_nobs : integer
        Minimum number of valid observations after flitering by the criteria of
        an estimation approach for a reliable estimate.

    Returns
    -------
    z_arr : vector of shape (nobs, 6)
        Output list of z estimates by the six approaches (labeled as in (Graf
        et al., 2014)) in the following order: FP-RE-1, FP-RE-2, FP-IT-1,
        FP-IT-2, FV-IT-1, FV-IT2

    z0_arr : vector of shape (nobs, 7)
        Output list of z0 estimates by the seven approaches (labeld as in (Graf
        et al., 2014)) in the following order: FP-RE-1, FP-RE-2, FP-IT-1,
        FP-IT-2, FV-IT-1-Eq6, FV-IT-2-Eq6, FV-IT-1-Eq14

    References
    ----------
    Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014.
    Intercomparison of Methods for the Simultaneous Estimation of Zero-Plane
    Displacement and Aerodynamic Roughness Length from Single-Level
    Eddy-Covariance Data. Boundary-Layer Meteorol 151, 373–387.
    https://doi.org/10.1007/s10546-013-9905-z
    %}

    % min_nobs = 10;

    u_ = data(:, 1);
    ustar_ = data(:, 2);
    L_ = data(:, 3);
    stdw_ = data(:, 4);
    stdT_ = data(:, 5);
    Tstar_ = data(:, 6);
    wd_ = data(:, 7);

    z_arr = zeros(size(data, 1), 6) + nan;
    z0_arr = zeros(size(data, 1), 7) + nan;
    nobs_arr = zeros(size(data, 1), 7);

    if wd_binsize <= 0
        nbins = 1;
        bin_beg = [0];
        bin_end = [360];
        half_wd_win = 180;
    else
        nbins = ceil(360./wd_binsize);
        bin_beg = wd_binsize*(0:nbins-1);
        bin_end = wd_binsize*(1:nbins); bin_end(end) = 360;
    end
    for i = 1:nbins
        wd_low = bin_beg(i) - half_wd_win;
        wd_wrap = wd_ - wd_low;
        wd_wrap(wd_wrap<0) = wd_wrap(wd_wrap<0) + 360;
        wd_wrap(wd_wrap>360) = wd_wrap(wd_wrap>360) - 360;
        i_sflag = wd_wrap >= 0 & wd_wrap < 2*half_wd_win;
        o_sflag = wd_ >= bin_beg(i) & wd_ < bin_end(i);
        if sum(i_sflag) < min_nobs || sum(o_sflag) < 1
            continue;
        end
        % 1st model set - linear simple and 3d regression 
        [z_arr(o_sflag, 1), z0_arr(o_sflag, 1), z_arr(o_sflag, 2), z0_arr(o_sflag, 2), nobs_arr(o_sflag, 1), _] = ...
            z0d_reg([u_(i_sflag) ustar_(i_sflag) L_(i_sflag)], ...
                z0max, zmax, min_nobs);
        nobs_arr(o_sflag, 2) = nobs_arr(o_sflag, 1);

        % 2nd set of models - iterative Martano   
        [z_arr(o_sflag, 4), z0_arr(o_sflag, 4), z_arr(o_sflag, 3), z0_arr(o_sflag, 3), nobs_arr(o_sflag, 4), _] = ...
            z0d_mart([u_(i_sflag) ustar_(i_sflag) L_(i_sflag)], ...
                z0max, zmax, zv, min_nobs);
        nobs_arr(o_sflag, 3) = nobs_arr(o_sflag, 4);

        % 3rd set of models - Weaver/Rotach/TodaSugita type (only d)
        [z_arr(o_sflag, 5), z_arr(o_sflag, 6), nobs_arr(o_sflag, 5), _] = ...
            z0d_varsim([stdw_(i_sflag) ustar_(i_sflag) L_(i_sflag) u_(i_sflag) stdT_(i_sflag) Tstar_(i_sflag)], ...
                zmax, zv, min_nobs);
        nobs_arr(o_sflag, 6) = nobs_arr(o_sflag, 5);
        [z0_arr(o_sflag, 5), _, _, _, _, _] = ...
            roughnesslength_revised3(z_arr(o_sflag, 5), L_(o_sflag), u_(o_sflag), ustar_(o_sflag), 0.4);
        [z0_arr(o_sflag, 6), _, _, _, _, _] = ...
            roughnesslength_revised3(z_arr(o_sflag, 6), L_(o_sflag), u_(o_sflag), ustar_(o_sflag), 0.4);

        % estimate z0 from variance(w) (after estimating d above)
        [z0_arr(o_sflag, 7), nobs_arr(o_sflag, 7), _] = ...
            z0_varsim_panofsky([stdw_(i_sflag) ustar_(i_sflag) L_(i_sflag) u_(i_sflag)], z_arr(5));
    end
end
