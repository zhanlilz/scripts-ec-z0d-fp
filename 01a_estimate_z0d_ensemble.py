#!/usr/bin/env python
#
# Estimate d, zero displacement height (or equivalently z, height above zero
# displacement plane) and z0, roughness length, using the six approaches in
# Graf et al., 2014.
# 
# Zhan Li, zhanli@gfz-potsdam.de

import os
import argparse
import configparser
import textwrap
import warnings
import tempfile
import subprocess

import numpy as np
import pandas as pd
from scipy import io as spio

def getCmdArgs():
    p = argparse.ArgumentParser(description='Simple script to read a CSV file of required variables and estimate roughness length and height above zero displacement plane using multiple approaches in Graf et al., 2014.')

    p.add_argument('-E', '--example_ini', dest='example_ini', action='store_true', help='Write an example INI file using the name given by PROCESS_CONTROL_FILE.')
    p.add_argument(dest='in_pcf', metavar='PROCESS_CONTROL_FILE', help='Process-Control File (PCF) in INI format that provides parameter values for estimating z0 (roughness length) and z (height above zeros displacement plane) and writing output CSV file.')

    p.add_argument('--graf_dir', dest='graf_dir', required=False, default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'z0d-graf-revision', 'code'), help='Directory to the Matlab/Octave .m files for z0d estimates by Graf et al., 2014. Default: [folder-of-this-script]/z0d-graf-revision/code')
    p.add_argument('--temporary_directory', dest='temporary_directory', required=False, default=None, help='If this option is given, the value in the PROCESS_CONTROL_FILE (PCF) in INI format will be overwritten. Default: None, read value from the PCF.')
    cmdargs = p.parse_args()
    return cmdargs

def getSampleIni():
    """Write a sample INI file
    
    Parameters
    ----------

    Returns
    -------
    
    """
    ini = """
    [meta_variables]
    
    [input_files]
    ; Input CSV file with the first row being header, each column being an
    ; input variable defined in the section of [input_variables].
    input_csv = /this/is/a/path/to/my/input/csv/file

    [input_variables]
    ; Column name of each required variable to estimate z0 and d.
    ; Column name of time stamp of one observation in a row
    col_timestamp = Timestamp

    ; A constant value for receptor/measurement height, meter
    receptor_height = heigh.value
    ; Column name of height of receptor/measurement, meter; If this column is
    ; given, the `receptor_height` value will be ignored. Leave the column name
    ; of receptor/measurement height empty, if you want to use a constant
    ; receptor/measurement height.
    col_receptor_height = zm

    ; Column name of mean wind speed, m*s^-1
    col_alongwind_speed = wind_speed
    
    ; Column name of friction velocity, m*s^-1
    col_friction_velocity = u_
    
    ; Column name of Monin-Obukhov length, meter
    col_obukhov_length = L

    ; Column name of standard deviation of vertical wind, m*s^-1
    col_verticalwind_sd = st_dev_w_

    ; Column name of standard deviation of sonic temperature, kelvin
    col_sonic_temperature_sd = st_dev_ts_

    ; Column name of Tstar
    col_tstar = T_
    
    ; Column name of wind direction, degree
    col_wind_direction = wind_dir

    ; Column name of one or multiple QC, unitless
    col_qc = qc_Tau
             qc_H

    ; One or multiple QC maximum, only process QC <= the given maximum for ALL
    ; the given QC columns. The QC maximums are in the order corresponding to
    ; QC columns by col_qc
    qc_max = 1 
             0

    ; Upper bound to the expected z0 values (e.g. 0.1h), in meters.
    maximum_z0 = 0.2

    ; Upper bound to the expected z values, in meters.
    maximum_z = 3

    ; Lower and upper bounds and steps to generate possible z values, in
    ; meters, for numerical search of optimal z and z0 in the FP-IT approaches
    ; and FV-IT approaches.
    z_lower = 0.1
    z_upper = 4
    z_step = 0.1
    
    ; Bin size of wind direction, in degrees. The z and z0 will be estimated
    ; per each bin to account for possible spatial variation along wind
    ; direction. If wd_binsize <=0, no data binning based on wind direction.
    ; All the input observations will be used in the estimation together.
    wind_direction_binsize = 1

    ; Half window size of wind directions, in degrees, to increase observation
    ; counts per each bin of wind direction and also smooth estimates along wind
    ; direction.  Observations in each bin of wind direction plus and mius this
    ; half window size will be used in the estimation of z and z0 for this bin.
    wind_direction_half_window = 21

    ; Half temporal window size, days. If >= 1 day, we group observations
    ; within a period of +/- temporal_half_window together in the estimation for
    ; each target date.
    ; 
    ; If empty, no temporal window is used and all the observations in the
    ; input CSV will be used together for estimation. All the available dates in
    ; the input CSV will be output. In this case, target dates will also be
    ; ignored. 
    temporal_half_window = 15

    ; Minimum number of valid observations after flitering by the criteria of
    ; an estimation approach for a reliable estimate. 
    minimum_num_obs = 10

    ; Target date or dates to process and output. You are responsible to make
    ; sure rows in the input CSV cover entire temporal window centered at all 
    ; the target dates.
    ; 
    ; If target dates are empty, estimates for every date in the input CSV file
    ; will be output.
    target_dates = 2018-01-01
                   2018-01-02

    [output_files]
    ; Output CSV of estimated z and z0
    output_csv = /this/is/the/path/to/my/output/csv/file 

    [user_runtime_parameters]
    ; Temporary directory
    temporary_directory = /this/is/the/path/to/my/temporary/directory

    """
    return textwrap.dedent(ini)

def main(cmdargs):
    in_pcf = cmdargs.in_pcf
    out_example = cmdargs.example_ini
    graf_dir = cmdargs.graf_dir
    temp_dir = cmdargs.temporary_directory

    if out_example:
        if os.path.isfile(in_pcf):
            msg_str = '''
            {0:s} exists! To use -E/--example_ini option to write an example
            INI to this file, delete it first! if you still need this file, did
            you accidentally type the option or file name wrong? 
            '''.format(in_pcf)
            msg_str = textwrap.dedent(msg_str)
            raise RuntimeError(msg_str)
        else:
            with open(in_pcf, 'w') as fobj:
                fobj.write(getSampleIni())
            return 0

    pcf = configparser.ConfigParser(allow_no_value=True)
    pcf.read(in_pcf)
    
    in_csv = pcf.get('input_files', 'input_csv')
    cols = dict()
    cols['ts'] = pcf.get('input_variables', 'col_timestamp')
    cols['zm'] = pcf.get('input_variables', 'col_receptor_height')
    if 'receptor_height' not in pcf.options('input_variables'):
        zm_const = ''
    else:
        zm_const = pcf.get('input_variables', 'receptor_height')
    if len(cols['zm'])>0 and len(zm_const)>0:
        msg_str = '''
        Both col_receptor_height and receptor_height are given! Ignore
        receptor_height and will read receptor height values from the column
        col_receptor_height.
        '''
        warnings.warn(textwrap.dedent(msg_str))
        zm_const = None
    elif len(cols['zm']) >0 and len(zm_const) <= 0:
        zm_const = None
    elif len(cols['zm']) <=0 and len(zm_const) >0:
        zm_const = float(zm_const)
        cols['zm'] = 'zm'
    else:
        msg_str = ''''
        Either col_receptor_height or receptor_height has to be given. Both are
        empty now. Correct your INI file and try again! 
        '''
        raise RuntimeError(textwrap.dedent(msg_str))
    cols['ws'] = pcf.get('input_variables', 'col_alongwind_speed')
    cols['ustar'] = pcf.get('input_variables', 'col_friction_velocity')
    cols['mo_len'] = pcf.get('input_variables', 'col_obukhov_length')
    cols['sd_w'] = pcf.get('input_variables', 'col_verticalwind_sd')
    cols['sd_ts'] = pcf.get('input_variables', 'col_sonic_temperature_sd')
    cols['tstar'] = pcf.get('input_variables', 'col_tstar')
    cols['wd'] = pcf.get('input_variables', 'col_wind_direction')
    cols['qc'] = pcf.get('input_variables', 'col_qc').split('\n')
    qc_max = pcf.get('input_variables', 'qc_max').split('\n')
    qc_max = [int(val) for val in qc_max]
    z0max = pcf.getfloat('input_variables', 'maximum_z0')
    zmax = pcf.getfloat('input_variables', 'maximum_z')
    zlower = pcf.getfloat('input_variables', 'z_lower')
    zupper = pcf.getfloat('input_variables', 'z_upper')
    zstep = pcf.getfloat('input_variables', 'z_step')
    wd_binsize = pcf.getfloat('input_variables', 'wind_direction_binsize')
    half_wd_win = pcf.getfloat('input_variables', 'wind_direction_half_window')
    half_time_win = None
    if 'temporal_half_window' in pcf.options('input_variables') \
            and len(pcf.get('input_variables', 'temporal_half_window'))>0:
        half_time_win = pcf.getfloat('input_variables', 'temporal_half_window')
    min_nobs = pcf.getint('input_variables', 'minimum_num_obs')
    target_dates = None
    if 'target_dates' in pcf.options('input_variables') \
            and len(pcf.get('input_variables', 'target_dates'))>0:
        target_dates = pcf.get('input_variables', 'target_dates').split('\n')
    out_csv = pcf.get('output_files', 'output_csv')
    if temp_dir is None:
        temp_dir = pcf.get('user_runtime_parameters', 'temporary_directory')
    try:
        # If the temporary directory does not exist, create it first.
        if not os.path.isdir(temp_dir):
            os.mkdir(temp_dir)
    except:
        raise RuntimeError('Failed attempting to create the temporary folder at' + temp_dir)

    if zm_const is None:
        usecols = sum([val if type(val) is list else [val] 
            for val in cols.values()], [])
    else:
        usecols = sum([val if type(val) is list else [val] 
            for val in [cols[k] for k in cols.keys() if k!='zm']], [])
    df = pd.read_csv(in_csv, index_col=cols['ts'], usecols=usecols)
    df.index = pd.DatetimeIndex(df.index)

    if target_dates is not None:
        target_dates = [pd.Timestamp(td) for td in target_dates]
    else:
        target_dates = pd.DatetimeIndex(np.unique((df.index-pd.Timedelta(1, unit='ms')).date))
    tdelta = pd.Timedelta(1, unit='D') - pd.Timedelta(1, unit='ms')
    with open(out_csv, 'w') as fobj:
        fobj.write('{0:s}\n'.format(','.join([ \
                df.index.name, cols['wd'], cols['zm'], \
                'z0-FP-RE-1', 'z0-FP-RE-2', \
                'z0-FP-IT-1', 'z0-FP-IT-2', \
                'z0-FV-IT-1-Eq6', 'z0-FV-IT-2-Eq6', 'z0-FV-IT-1-Eq14', \
                'z-FP-RE-1', 'z-FP-RE-2', \
                'z-FP-IT-1', 'z-FP-IT-2', \
                'z-FV-IT-1', 'z-FV-IT-2', \
                'd-FP-RE-1', 'd-FP-RE-2', \
                'd-FP-IT-1', 'd-FP-IT-2', \
                'd-FV-IT-1', 'd-FV-IT-2', \
                'nobs-FP-RE-1', 'nobs-FP-RE-2', \
                'nobs-FP-IT-1', 'nobs-FP-IT-2', \
                'nobs-FV-IT-1', 'nobs-FV-IT-2', 'nobs-FV-IT-1-Eq14'])))

        if half_time_win is None:
            warningstr = 'Temporal half window is empty. ' \
                    + 'All input observations are used together\n' \
                    + 'Target dates are ignored!'
            warnings.warn(warningstr)
        for td in target_dates:
            if half_time_win is None:
                beg_dt = df.index.min()
                end_dt = df.index.max()
            else:
                if (len(df.loc[td:td+tdelta, :])==0):
                    warnings.warn('Target date {0:s} is not among the dates in the input CSV file'.format(td.strftime('%Y-%m-%d')))
                    continue
                beg_dt = pd.Timestamp(td) - pd.Timedelta(half_time_win, unit='D')
                end_dt = pd.Timestamp(td) + pd.Timedelta(half_time_win+1, unit='D') - pd.Timedelta(1, unit='ms')
            cur_df = df.loc[beg_dt:end_dt, :]

            sflag = np.ones(cur_df.shape[0], dtype=bool)
            for qcval, qmval in zip(cols['qc'], qc_max):
                sflag = np.logical_and(sflag, cur_df[qcval] <= qmval)
            var_list = ['ws', 'wd', 'ustar', 'mo_len', 'sd_w', 'sd_ts', 'tstar']
            if zm_const is None:
                var_list += ['zm']
            for var in var_list:
                sflag = np.logical_and(sflag, np.logical_not(np.isnan(cur_df[cols[var]].values)))

            if np.sum(sflag) == 0:
                warnstr = 'From {0:s} to {1:s}, no observations meet all the following QC criteria'.format(beg_dt.strftime('%Y-%m-%d'), end_dt.strftime('%Y-%m-%d'))
                for qcval, qmval in zip(cols['qc'], qc_max):
                    warnstr = warnstr + '\n  {0:s} <= {1:d}'.format(qcval, qmval)
                warnings.warn(warnstr)
                continue
            colnames = [cols[val] for val in ['ws', 'ustar', 'mo_len', 'sd_w', 'sd_ts', 'tstar', 'wd']]
            i_mat_dict = dict( \
                    data=cur_df[colnames].values[sflag, :], \
                    z0max=z0max, zmax=zmax, zv=np.arange(zlower, zupper, zstep), \
                    wd_binsize=wd_binsize, half_wd_win=half_wd_win, \
                    min_nobs=min_nobs)

            with tempfile.NamedTemporaryFile(suffix='.mat', prefix='tmp.', dir=temp_dir) as i_mat_fobj, \
                tempfile.NamedTemporaryFile(suffix='.mat', prefix='tmp.', dir=temp_dir) as o_mat_fobj:
                i_mat_file = i_mat_fobj.name
                o_mat_file = o_mat_fobj.name
                # i_mat_file = tempfile.mkstemp(suffix='.mat', prefix='tmp.', dir=temp_dir)[1]
                # o_mat_file = tempfile.mkstemp(suffix='.mat', prefix='tmp.', dir=temp_dir)[1]
                spio.savemat(i_mat_file, i_mat_dict, format='4')
                octave_script_str = '''
                warning("off");
                addpath("{2:s}");
                load("{0:s}");
                [z_arr, z0_arr, nobs_arr] = ...
                    z0d_ensemble( ...
                        data, ...
                        z0max, ...
                        zmax, ...
                        zv, ...
                        wd_binsize, ...
                        half_wd_win, ...
                        min_nobs);
                save("{1:s}", "z_arr", "z0_arr", "nobs_arr", "-v4");
                '''.format(i_mat_file, o_mat_file, graf_dir)
                subprocess.run(['octave', '-q', '--eval', octave_script_str])
                o_mat_dict = spio.loadmat(o_mat_file)

                z0_arr = np.zeros((len(cur_df), 7)) + np.nan
                z_arr = np.zeros((len(cur_df), 6)) + np.nan
                nobs_arr = np.zeros((len(cur_df), 7))
                
                z0_arr[sflag] = o_mat_dict['z0_arr']
                z_arr[sflag] = o_mat_dict['z_arr']
                nobs_arr[sflag] = o_mat_dict['nobs_arr']
                d_arr = cur_df[cols['zm']].values[:, np.newaxis] if zm_const is None else zm_const
                d_arr = d_arr - z_arr
                
                out_df = pd.DataFrame(np.hstack([\
                        cur_df[[cols['wd']]].values, \
                        cur_df[[cols['zm']]].values if zm_const is None else np.ones((len(cur_df), 1))*zm_const, \
                        z0_arr, z_arr, d_arr, nobs_arr]), \
                        index=cur_df.index)
                if half_time_win is None:
                    out_df.to_csv(fobj, mode='a', header=False, 
                            na_rep="NaN", date_format="%Y-%m-%d %H:%M")
                    break
                else:
                    out_df.loc[td:td+tdelta, :].to_csv(fobj, mode='a', header=False, 
                            na_rep="NaN", date_format="%Y-%m-%d %H:%M")

                # os.remove(i_mat_file)
                # os.remove(o_mat_file)

if __name__ == "__main__":
    cmdargs = getCmdArgs()
    main(cmdargs)
