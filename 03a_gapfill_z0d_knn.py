#!/usr/bin/env python

import os
import argparse
import configparser
import textwrap

import numpy as np
import pandas as pd

from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit

def getCmdArgs():
    p = argparse.ArgumentParser(description='Fill gaps in the estimates of z0 and d using kNN in the spatiotemporal space of wind directions (spatial) and dates (temporal).')

    p.add_argument('-E', '--example_ini', dest='example_ini', action='store_true', help='Write an example INI file using the name given by PROCESS_CONTROL_FILE.')
    p.add_argument(dest='in_pcf', metavar='PROCESS_CONTROL_FILE', help='Process-Control File (PCF) in INI format that provides parameter values for filling gaps in the estimates of z0 and, and writing output CSV file.')

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
    ; *** Required ***
    input_csv = /this/is/a/path/to/my/input/csv/file

    [input_variables]
    ; Column name of time stamp of one observation in a row
    ; *** Required ***
    col_timestamp = Timestamp
    
    ; Column name of wind direction, degree
    ; *** Required ***
    col_wind_direction = wind_dir

    ; Column names (one or multiple) of target values to be interpolated and
    ; filled. Each column/target variable is interpolated separately but using
    ; the same parameterized kNN and written into the same output CSV.
    ; *** Required ***
    col_y = z0
            d

    ; Column names (one or multiple) of the QA of target values to be
    ; interpolated and filled. Must be the same size as col_y and the same order
    ; of target variables in col_y 
    ; *** Required ***
    col_y_qa = z0_qa
               d_qa

    ; Names of additional columns in the input CSV to be copied into the output
    ; CSV.
    ; *** Optional ***
    col_y_copy = zm_above_water
                 another_column_to_copy

    ; Weights on time and wind direciton to calculate weighted L2-norm
    ; Minkowski distance in the search of nearest neighbors.
    ; *** NOT IMPLEMENTED ***
    weight_timestamp = 0.5
    weight_wind_direction = 0.5

    ; Number of neighbors for interpolation
    ; *** Required ***
    n_neighbors = 30

    ; Folds for reporting cross-validation scores. If not given, no
    ; cross-validation is performed and no report is printed on screen.
    ; *** Optional ***
    cv_folds = 5

    [output_files]
    ; Output CSV of gap-filled z0 and d 
    ; *** Required ***
    output_csv = /this/is/the/path/to/my/output/csv/file 

    [user_runtime_parameters]

    """
    return textwrap.dedent(ini)


def main(cmdargs):
    in_pcf = cmdargs.in_pcf
    out_example = cmdargs.example_ini
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
    cols['wd'] = pcf.get('input_variables', 'col_wind_direction')
    cols['y'] = pcf.get('input_variables', 'col_y').split('\n')
    cols['y_qa'] = pcf.get('input_variables', 'col_y_qa').split('\n')
    if 'col_y_copy' in pcf.options('input_variables') \
            and len(pcf.get('input_variables', 'col_y_copy'))>0:
        cols['y_copy'] = pcf.get('input_variables', 'col_y_copy').split('\n')
    
    wt = dict()
    if 'ts' in pcf.options('input_variables') \
        and len(pcf.get('input_variables', 'ts'))>0:
        wt['ts'] = pcf.getfloat('input_variables', 'weight_timestamp')
    if 'wd' in pcf.options('input_variables') \
        and len(pcf.get('input_variables', 'ts'))>0:
        wt['wd'] = pcf.getfloat('input_variables', 'weight_wind_direction')

    n_neighbors = pcf.getint('input_variables', 'n_neighbors')
    cv_folds = None
    if 'cv_folds' in pcf.options('input_variables') \
            and len(pcf.get('input_variables', 'cv_folds'))>0:
        cv_folds = pcf.getint('input_variables', 'cv_folds')
    out_csv = pcf.get('output_files', 'output_csv')

    cols2use = [cols[k] for k in cols.keys() if k!='y' and k!='y_qa' and k!='y_copy'] \
            + cols['y'] + cols['y_qa'] + cols['y_copy']
    df = pd.read_csv(in_csv, comment='#', usecols=cols2use, parse_dates=[cols['ts']], 
            index_col=cols['ts'])
    # Use both dates and time of the day to measure temporal distance between
    # two observations. 
    tod = df.index.time
    tod = [val.hour*3600+val.minute*60+val.second for val in tod]
    tod_sin, tod_cos = zip(*[(np.cos(val/3600.*2*np.pi), np.sin(val/3600.*2*np.pi)) for val in tod])
    jd = pd.to_datetime(df.index.date).to_julian_date().values / 183.
    X = [jd, tod_cos, tod_sin, \
            np.cos(np.deg2rad(df[cols['wd']].values)), \
            np.sin(np.deg2rad(df[cols['wd']].values))]

    X = np.array(X).T
#     wt_sum = np.sum(tuple(wt.values()))
#     wt_X = np.array([0.5*wt['ts']/wt_sum, \
#             0.25*wt['ts']/wt_sum, 0.25*wt['ts']/wt_sum, \
#             0.5*wt['wd']/wt_sum, 0.5*wt['wd']/wt_sum])
    wt_X = np.ones(X.shape[1]) 

    for ycol, qacol in zip(cols['y'], cols['y_qa']):
        y = df[ycol].values
        sflag2fill = np.isnan(y)
        sflag2train = np.logical_not(sflag2fill)
        sflag_X = np.logical_not(np.any(np.isnan(X), axis=1))
        sflag2train = np.logical_and(sflag2train, sflag_X)
        sflag2fill = np.logical_and(sflag2fill, sflag_X)

        estimator = KNeighborsRegressor( \
                n_neighbors=n_neighbors, \
                weights='distance', \
                metric='wminkowski', \
                p = 2, metric_params=dict(w=wt_X))
        if cv_folds is not None:
            scores = cross_val_score(estimator, 
                    X[sflag2train], y[sflag2train],
                    cv=ShuffleSplit(n_splits=cv_folds, 
                        test_size=1./cv_folds, 
                        random_state=0))
            msg_str = "{4:d}-fold cross-validation R2 for interpolating {0:s}" \
                    + " using {1:d}-NN on wind direction and dates: " \
                    + "{2:.2f} +/- {3:.2f} (95% CI)."
            msg_str = msg_str.format( \
                    ycol, n_neighbors, scores.mean(), scores.std() * 2, cv_folds)
            print(msg_str)
        estimator.fit(X[sflag2train, :], y[sflag2train])
        y[sflag2fill] = estimator.predict(X[sflag2fill, :])
        df[ycol] = y
        # If any gap cannot be filled by this KNN procedure and still remains a
        # gap (NaN), increase its QA value by 1.
        df.loc[np.isnan(y), qacol] += 1

    df.to_csv(out_csv, header=True, na_rep="NaN", date_format="%Y-%m-%d %H:%M")

if __name__ == '__main__':
    cmdargs = getCmdArgs()
    main(cmdargs)
