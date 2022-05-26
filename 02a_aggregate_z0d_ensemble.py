#!/usr/bin/env python
# 
# Aggregate the estimates of z0 and d from the multiple approaches given in
# Graf et al., 2014.
# 
# Zhan Li, zhanli@gfz-potsdam.de
# Created: Fri Sep 18 15:08:19 CEST 2020

import argparse
import re
import textwrap

import numpy as np
import pandas as pd

def getCmdArgs():
    desc_str = '''
    Aggregate estimates of z0 and d from the multiple approaches given by Graf
    et al., 2014 after passing the following filtering criteria in the given
    order for plausibility check.
    1. z0 > 0, zm >= d > 0
    2. If user-defined bounds of z0 and d are given, z0_upper_bound >= z0 >=
       z0_lower_bound, d_upper_bound >= d >= d_lower_bound
    3. z0 within mean +/- 3*std, mean and std of all the input rows
    4. d within mean +/- 3*std, mean and std of all the input rows
    '''
    p = argparse.ArgumentParser(description=textwrap.dedent(desc_str), formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--cols', metavar='COLUMN_NAME', nargs=3, required=False, default=['Timestamp', 'wind_dir', 'zm'], help='Column names (separated by space) in this order for timestamps, wind direction, and receptor/measurement height. Default: Timestamp wind_dir zm.')
    p.add_argument('in_csv', metavar='INPUT_CSV', help='Path to the input CSV file of multiple estimates of z0 and d using different approaches.')
    p.add_argument('out_csv', metavar='OUTPUT_CSV', help='Path to the output CSV file of aggregated ensemble estimates of z0 and d after plausability filtering.')

    p.add_argument('--z0_bound', metavar='Z0_BOUND', nargs=2, type=float, 
        required=False, default=[None, None], 
        help='Lower and upper bounds of plausible z0 values (inclusive). Default: [None, None], no checking of plausible bounds.')
    p.add_argument('--d_bound', metavar='D_BOUND', nargs=2, type=float, 
        required=False, default=[None, None], 
        help='Lower and upper bounds of plausible d values (inclusive). Default: [None, None], no checking of plausible bounds.')
    cmdargs = p.parse_args()
    return cmdargs

def main(cmdargs):
    in_csv = cmdargs.in_csv
    out_csv = cmdargs.out_csv
    colnames = cmdargs.cols

    z0_bound = cmdargs.z0_bound
    d_bound = cmdargs.d_bound

    z0_re = re.compile('z0-')
    z_re = re.compile('z-')
    d_re = re.compile('d-')
    nobs_re = re.compile('nobs-')
    with open(in_csv, 'r') as fobj:
        linestr = fobj.readline().strip()
    cols = linestr.split(',')
    z0_cols = [val for val in cols if z0_re.match(val)]
    z_cols = [val for val in cols if z_re.match(val)]
    d_cols = [val for val in cols if d_re.match(val)]
    nobs_cols = [val for val in cols if nobs_re.match(val)]
    ts_col, wd_col, zm_col = colnames[0], colnames[1], colnames[2]

    # For the filter using d/z0 ratio in the range (1, 10)
    d_over_z0_pairs = [
            ('d-FP-RE-1', 'z0-FP-RE-1'),
            ('d-FP-RE-2', 'z0-FP-RE-2'),
            ('d-FP-IT-1', 'z0-FP-IT-1'),
            ('d-FP-IT-2', 'z0-FP-IT-2'),
            ('d-FV-IT-1', 'z0-FV-IT-1-Eq6'),
            ('d-FV-IT-2', 'z0-FV-IT-2-Eq6'),
            ('d-FV-IT-1', 'z0-FV-IT-1-Eq14')]

    with open(out_csv, 'w') as fobj:
        fobj.write('{0:s}\n'.format('# QA = 0: at least two valid estimates from the approaches in Graf et al. 2014'))
        fobj.write('{0:s}\n'.format('# QA = 1: at least one valid estimates from the approaches in Graf et al. 2014'))
        fobj.write('{0:s}\n'.format('# QA = 2: no valid estimates from the approaches in Graf et al. 2014'))
        fobj.write('{0:s}\n'.format('# A valid estimate is a value successfully estimated by an approach and this value is within mean +/- 3*std of all the estimates (rows) by this approach.'))
        fobj.write('{0:s}\n'.format(','.join([
                ts_col, wd_col, zm_col, 
                'z0_ensemble', 'z0_ensemble_qa', 
                'z_ensemble', 'z_ensemble_qa', 
                'd_ensemble', 'd_ensemble_qa', 
                'nobs4z0_ensemble', 'nobs4z_ensemble', 'nobs4d_ensemble',
                'nest4z0', 'nest4z', 'nest4d'])))
        
        chunk_df = pd.read_csv(in_csv, index_col=ts_col)

#         for dcol, z0col in d_over_z0_pairs:
#             tmp = chunk_df[dcol] / chunk_df[z0col]
#             sflag = np.logical_and(tmp>0, tmp<40)
#             chunk_df[dcol] = chunk_df[dcol].where(sflag, np.nan)
#             chunk_df[z0col] = chunk_df[z0col].where(sflag, np.nan)

        chunk_df[z0_cols] = chunk_df[z0_cols].where(chunk_df[z0_cols]>0, np.nan)
        # User-defined plausible bounds for z0
        if z0_bound[0] is not None:
            chunk_df[z0_cols] = chunk_df[z0_cols].where(chunk_df[z0_cols]>=z0_bound[0], np.nan)
        if z0_bound[1] is not None:
            chunk_df[z0_cols] = chunk_df[z0_cols].where(chunk_df[z0_cols]<=z0_bound[1], np.nan)
        # Further filtering outliers in the z0 estimates 
        # Outliers are outside annual mean +/- 3*std
        tmp_mean = chunk_df[z0_cols].mean(axis=0).values
        tmp_std = chunk_df[z0_cols].std(axis=0).values
        chunk_df[z0_cols] = chunk_df[z0_cols].where( \
                chunk_df[z0_cols].values > tmp_mean-3*tmp_std, np.nan)
        chunk_df[z0_cols] = chunk_df[z0_cols].where( \
                chunk_df[z0_cols].values < tmp_mean+3*tmp_std, np.nan)

        z0_esb = chunk_df[z0_cols].median(axis=1, skipna=True)
        nest4z0 = chunk_df[z0_cols].count(axis=1)
        nobs4z0_esb = chunk_df[nobs_cols]\
                .where(chunk_df[z0_cols].notna().values, np.nan)\
                .median(axis=1, skipna=True)
        z0_esb_qa = np.zeros_like(nest4z0) + 2
        z0_esb_qa[nest4z0>0] = 1
        z0_esb_qa[nest4z0>1] = 0
        z0_esb_qa = pd.Series(z0_esb_qa, index=nest4z0.index)

        chunk_df[d_cols] = chunk_df[d_cols].where(chunk_df[d_cols]>0, np.nan)
        chunk_df[d_cols] = chunk_df[d_cols].where(chunk_df[d_cols].values<chunk_df[[zm_col]].values, np.nan)
        # User-defined plausible bounds for d
        if d_bound[0] is not None:
            chunk_df[d_cols] = chunk_df[d_cols].where(chunk_df[d_cols]>=d_bound[0], np.nan)
        if d_bound[1] is not None:
            chunk_df[d_cols] = chunk_df[d_cols].where(chunk_df[d_cols]<=d_bound[1], np.nan)
        # Further filtering outliers in the d estimates 
        # Outliers are outside annual mean +/- 3*std
        tmp_mean = chunk_df[d_cols].mean(axis=0).values
        tmp_std = chunk_df[d_cols].std(axis=0).values
        chunk_df[d_cols] = chunk_df[d_cols].where( \
                chunk_df[d_cols].values > tmp_mean-3*tmp_std, np.nan)
        chunk_df[d_cols] = chunk_df[d_cols].where( \
                chunk_df[d_cols].values < tmp_mean+3*tmp_std, np.nan)

        sflag = chunk_df[d_cols].notna().values

        chunk_df[z_cols] = chunk_df[z_cols].where(sflag, np.nan)
        z_esb = chunk_df[z_cols].median(axis=1, skipna=True)
        
        d_esb = chunk_df[d_cols].median(axis=1, skipna=True)
        nest4d = chunk_df[d_cols].count(axis=1)
        nobs4d_esb = chunk_df[nobs_cols[:-1]]\
                .where(sflag, np.nan)\
                .median(axis=1, skipna=True)
        d_esb_qa = np.zeros_like(nest4d) + 2
        d_esb_qa[nest4d>0] = 1
        d_esb_qa[nest4d>1] = 0
        d_esb_qa = pd.Series(d_esb_qa, index=nest4d.index)

        pd.concat([
                chunk_df[[wd_col, zm_col]], 
                z0_esb, z0_esb_qa, 
                z_esb, d_esb_qa, 
                d_esb, d_esb_qa, 
                nobs4z0_esb.round().astype('Int64'), 
                nobs4d_esb.round().astype('Int64'),
                nobs4d_esb.round().astype('Int64'), 
                nest4z0, nest4d, nest4d], axis=1).to_csv(
                        fobj, mode='a', header=False, 
                        na_rep="NaN", date_format="%Y-%m-%d %H:%M")
            
    return

if __name__ == '__main__':
    cmdargs = getCmdArgs()
    main(cmdargs)
