- [Scripts to generate flux footprint rasters given single-level micrometeorological measurements.](#scripts-to-generate-flux-footprint-rasters-given-single-level-micrometeorological-measurements)
  - [ONLY TESTED ON LINUX. EXTENSIVE DEBUGGINNG IS NEEDED TO RUN ON WINDOWS](#only-tested-on-linux-extensive-debugginng-is-needed-to-run-on-windows)
  - [References](#references)
  - [Overview](#overview)
  - [Quick start](#quick-start)
    - [Preparation](#preparation)
    - [Step 1: estimate z0d using multiple approaches in Graf et al. 2014](#step-1-estimate-z0d-using-multiple-approaches-in-graf-et-al-2014)
    - [Step 2: aggregate multiple estimates of z0 and d per half-hourly record](#step-2-aggregate-multiple-estimates-of-z0-and-d-per-half-hourly-record)
    - [Step 3: fill gaps in z0d using wind direction and time in a K-Nearest-Neighbor algorithm](#step-3-fill-gaps-in-z0d-using-wind-direction-and-time-in-a-k-nearest-neighbor-algorithm)
    - [Step 4: estimate flux footprints using gap-filled z0d and micrometeorological data and save them in a NetCDF file](#step-4-estimate-flux-footprints-using-gap-filled-z0d-and-micrometeorological-data-and-save-them-in-a-netcdf-file)
    - [Deal with long-running commands](#deal-with-long-running-commands)

# Scripts to generate flux footprint rasters given single-level micrometeorological measurements.

Zhan Li, zhan.li@gfz-potsdam.de

## ONLY TESTED ON LINUX. EXTENSIVE DEBUGGINNG IS NEEDED TO RUN ON WINDOWS

## References
* Graf, A., van de Boer, A., Moene, A., Vereecken, H., 2014. Intercomparison of
Methods for the Simultaneous Estimation of Zero-Plane Displacement and
Aerodynamic Roughness Length from Single-Level Eddy-Covariance Data.
Boundary-Layer Meteorol 151, 373–387. https://doi.org/10.1007/s10546-013-9905-z
* Kormann, R., Meixner, F.X., 2001. An Analytical Footprint Model For Non-Neutral
Stratification. Boundary-Layer Meteorology 99, 207–224.
https://doi.org/10.1023/A:1018991015119

## Overview
The procedure has two componenets:
1. Estimate gap-filled roughness lengths (z0) and zero-displacement heights (d)
using single-level micrometeorological measurements, achieved by the scripts,
    * 01a_estimate_z0d_ensemble.py
    * 02a_aggregate_z0d_ensemble.py
    * 03a_gapfill_z0d_knn.py
2. Generate flux footprint rasters in NetCDF format using gap-filled z0&d and
also micrometeorological measurements, achieved by the script, 
    * 04a_gen_footprints_netcdf.sh

## Quick start

### Preparation
* The scripts are developed and only tested on Linux system and under Bash
  shell. Please make sure you are on a Linux system **and** also are using Bash as
  your terminal shell. If you are on one of the mefe servers at GFZ, you need to
  change your Unix Shell to `/bin/bash`. To make this change, you need to log
  into your GFZ profile, then under the tab *Accounts*, you will see the option
  of *Unix Shell* and please change it to `/bin/bash`.
* Copy this folder to your desired location. Then change directory to this folder. 
    ```bash
    $ cd scripts-ec-z0d-fp
    ```
* Unpack the zipped file of dependent libraries for our scripts and set up the
  environmentl to run the scripts.
    ```bash
    $ mkdir ec-z0d-fp-env
    $ tar -C ec-z0d-fp-env -xzf ec-z0d-fp-env.tar.gz
    $ source ec-z0d-fp-env/bin/activate
    $ conda-unpack
    ```
* If you exit the terminal and come back to the same folder of the scripts, you
  need to **re-run the following two commands again to set up your environment**, 
    ```bash
    $ source ec-z0d-fp-env/bin/activate
    $ conda-unpack
    ```

* Run the following commands to double check to see if your environment is ready
  to run the scripts. If you see "yay!" four times without any error messages,
  you are ready to go. 
    ```bash
    $ python 01a_estimate_z0d_ensemble.py -h > /dev/null && echo "yay!"
    $ python 02a_aggregate_z0d_ensemble.py -h > /dev/null && echo "yay!"
    $ python 03a_gapfill_z0d_knn.py -h > /dev/null && echo "yay!"
    $ export PYTHONPATH="./fluxfm"
    $ python fluxfm/fluxfm/cli/estimate_footprint_km.py -h > /dev/null && echo "yay!"
    ```

* Now you are ready. We will use the example input files in the folder of
  `./tests/data` to walk through the quick start. All the following steps are
  also lumped into a single bash shell script `000_example_ec_z0d_fp.sh` that
  briefs the entire workflow. 

### Step 1: estimate z0d using multiple approaches in Graf et al. 2014
We need to run the script `01a_estimate_z0d_ensemble.py` in Step-1. To get help
about this script, 
```bash
$ python 01a_estimate_z0d_ensemble.py -h
```

This script pools half-hourly records from a prescribed temporal window (days)
and a prescribed wind-direction window (degrees) to estimate z0 and d using
multiple approaches. 

*You are responsible to make sure records (rows) in the input CSV cover entire
temporal window centered at all the target dates on which you would like to
estimate z0 and d*.

**Tasks in Step-1**
* Prepare a CSV file like `./tests/data/test_01a_estimate_z0d_ensemble_input.csv`
  where the first row contains the names of columns and each row is a
  half-hourly record that has at least the following columns, 
  * Time stamp of the half-hourly record in this row
  * Height of receptor/measurement, meter
  * Mean wind speed, m*s^-1
  * Friction velocity, m*s^-1
  * Monin-Obukhov length, meter
  * Standard deviation of vertical wind, m*s^-1
  * Standard deviation of sonic temperature, kelvin
  * Tstar
  * Wind direction, degree
  * QC, unitless

* Prepare an INI file of configuration for the estimation of z0 and d. See
  `./tests/data/test_01a_estimate_z0d_ensemble_config.ini` for explanations and
  as an example to follow. Or you can run the following command to get a sample
  INI configuration file that explains the required contents of the INI file, 
    ```bash
    $ python 01a_estimate_z0d_ensemble.py --example_ini path/to/my_sample_configuration_file.ini
    ```
* Run the following command to estimate z0 and d using multiple approaches and
  generate an output CSV file that should look like
  `./tests/data/example-outputs/test_01a_estimate_z0d_ensemble_output.csv` 
    ```bash
    $ python 01a_estimate_z0d_ensemble.py ./tests/data/test_01a_estimate_z0d_ensemble_config.ini
    ```

### Step 2: aggregate multiple estimates of z0 and d per half-hourly record
We need to run the script `02a_aggregate_z0d_ensemble.py` in Step-2. To get help
about this script, 
```bash
$ python 02a_aggregate_z0d_ensemble.py -h
```

This script filters the multi-approach estimates per record and aggregate
(average) valid estimates after filtering. The filtering first check whether a
value successfully estimated by an approach (i.e., not a NaN) is within mean +/-
3*std of all the estimates (rows) by this approach. Only values within this
range are valid estimates for this particular approach (represented by a
column). Then remaining valid estimates per record (row) are averaged with the
following QA flags assigned, 
* QA = 0: at least two valid estimates from the approaches in Graf et al. 2014
* QA = 1: at least one valid estimates from the approaches in Graf et al. 2014
* QA = 2: no valid estimates from the approaches in Graf et al. 2014

**Tasks in Step-2**
* Run the following command to aggregate multiple estimates of z0 and d by
  different approaches per half-hourly record. This should generate a CSV file
  like `tests/data/example-outputs/test_02a_aggregate_z0d_ensemble_output.csv`
  that contains multi-approach ensemble estimates of z0 and d. 
    ```bash
    # The character '\' is just a line continuatioin character in the command
    # line that allows you to simply type Enter key to break a long command line
    # into muliple lines. Make sure no space after you type '\' and before type
    # Enter key to continue to the next line.
    $ python 02a_aggregate_z0d_ensemble.py \
        --cols Timestamp wind_dir zm \
        -- \
        ./tests/data/test_01a_estimate_z0d_ensemble_output_by_user_run.csv \
        ./tests/data/test_02a_aggregate_z0d_ensemble_output_by_user_run.csv
    ```
> NOTE: Do **NOT** forget to put the characters `--` before the last two
positional arguments.

### Step 3: fill gaps in z0d using wind direction and time in a K-Nearest-Neighbor algorithm
We need to run the script `03a_gapfill_z0d_knn.py` in Step-3. To get help about
this script, 
```bash
$ python 03a_gapfill_z0d_knn.py -h
```

This script fills gaps in the ensemble estimates of z0 and d using the estimates
from the neighboring records of a gap record. The neighors are defined according
to space (wind direction) and time (time of a day and date). 

*You prescribe how many neighors should be used to fill a gap record. You are
responsible to have enough half-hourly records in the CSV file to ensur enough
temporally-close neighors around the records of interests for a good gap
filling.*

Currently it does **NOT** limit the temporal distance (days) between neighors
and a gap to fill. It also does **NOT** limit the spatial distance
(wind-direction differences) between neighors and a gap to fill.

**Tasks in Step-3**
* Prepare an INI file of configuration to run the gap filling. See
  `./tests/data/test_03a_gapfill_z0d_knn_config.ini` for explanations and as an
  example to follow. Or you can run the following command to get a sample INI
  configuration file that explains the required contents of the INI file, 
    ```bash
    $ python 03a_gapfill_z0d_knn.py --example_ini path/to/my_sample_configuration_file.ini
    ```
* Run the following command for gap filling. It should generate an output CSV
  file like
  `./tests/data/example-outputs/test_02a_aggregate_z0d_ensemble_output_gapfilled.csv`.
    ```bash
    $ python 03a_gapfill_z0d_knn.py ./tests/data/test_03a_gapfill_z0d_knn_config.ini
    ```

### Step 4: estimate flux footprints using gap-filled z0d and micrometeorological data and save them in a NetCDF file
We need to run the script `` in Step-4. To get help about this script, 
```bash
$ bash 04a_gen_footprints_netcdf.sh -h
```

This script is a helper script in bash that estimate footprints of multiple
half-hourly records in an input CSV file in parallel and save the footprint
rasters of all the half-hourly records in the input CSV file in a NetCDF file.
The actual estimation of the footprint per a record is done by the Python script
`./fluxfm/fluxfm/cli/estimate_footprint_km.py`. The Python script estimates
footprints using Korman-Meixner-2001 model. To get help about this Python
script, 
```bash
$ export PYTHONPATH="./fluxfm"
$ python fluxfm/fluxfm/cli/estimate_footprint_km.py -h
```

**Tasks in Step-4**
* Prepare a CSV file like
  `./tests/data/test_04a_gen_footprints_netcdf_input.csv` that combines the
  gap-filled ensemble estimates of z0 and d and the other required
  micrometeorological data. This CSV file will be the inputs for the footprint
  estimates. In this CSV file,  each row is a half-hourly record and at minimum
  contains the following columns with the
  **EXACT names**: `Timestamp, qc_Tau, z_ensemble, z0_ensemble, wind_speed, u_,
  L, st_dev_v_, wind_dir`.  where, 
    * `Timestamp`: timestamp of the half-hourly record.  
    * `qc_Tau`: QC flag by turbulance condition of this half-hourly record. If
    qc_Tau >= 2, no footprint is modeled.
    * `z_ensemble`: Aerodynamic/effective height of receptor/measurement (i.e.,
    height above zero-displacement plane), meter
    * `z0_ensemble`: Roughness length, meter
    * `wind_speed`: Mean wind speed, m*s^-1
    * `u_`: Friction velocity, m*s^-1
    * `L`: Monin-Obukhov length, meter
    * `st_dev_v_`: Standard deviation of cross-wind speed, m*s^-1
    * `wind_dir`: wind direction, degree, measured with regard to true north.

* (Optional but recommended) Figure out the settings of the following three
  options for the script `04a_gen_footprints_netcdf.sh`. 
    * --grid_srs: spatial reference system (SRS, that is a coordinate system) into which you want to save your footprint raster grids. 
    * --grid_bounds: the bounds (xmin, xmax, ymin, ymax) of the raster extent that should cover the full extents of most footprints. Coordinates are in the chosen SRS.
    * --receptor_xy: the coordinates of the receptor in the chosen SRS.

    *For Zarnekow, we use `--grid_res="epsg:32633"`, `--grid_bounds="360725.0
    361725.0 5971285.0 5972285.0"`, and `--receptor_xy="361224.023952403
    5971784.23062857"`* 

* Run the following command. It should generate an NetCDF file like
  `./tests/data/example-outputs/test_04a_gen_footprints_netcdf_output.nc` that
  contains the footprint rasters of all the half-hourly records in the input CSV
  file. 
    ```bash
    $ bash 04a_gen_footprints_netcdf.sh \
        --grid_srs="epsg:32633" \
        --grid_bounds="360725.0 361725.0 5971285.0 5972285.0" \
        --receptor_xy="361224.023952403 5971784.23062857" \
        --run_res=0.1 \
        --out_res=1 \
        --njobs=0 \
        ./tests/data/test_04a_gen_footprints_netcdf_input.csv \
        ./tests/data/test_04a_gen_footprints_netcdf_output_by_user_run.nc
    ```

    *The script `04a_gen_footprints_netcdf.sh` usually produces large amounts of
    temporary data. Therefore, it is a good idea to specify your own location
    where the script can dump temporary data (automatically deleted afterwards
    unless some error occurs) using the option
    `--tmp_dir="path/to/my/own/temporary/directory"`*

> NOTE: Remember to **double quote** the values to the options `--grid_srs`,
`--grid_bounds`, and `--receptor_xy`

### Deal with long-running commands
Some of the above commands may take long time to finish, esp. the scripts of
step-1 and step-4. When running such scripts on the Linux server, you may use
the command `nohup` to start the processing and keep the process running on
the server even after you log off from the linux server. 
```bash
$ nohup python 01a_estimate_z0d_ensemble.py ./tests/data/test_01a_estimate_z0d_ensemble_config.ini > my_program_log.txt &
```
Between `nohup` and `>`, you put whatever command and your input arguments there
like you would run on a command line in your terminal. `> my_program_log.txt`
tells `nohup` to save all the ouputs by the command to the terminal in this text
file called `my_program_log.txt` or whatever you like to name. **Do NOT forget
the charater `&` in the end**. Then you can log off the server and process will
keep running on the server.

To check if you process is running after you log off the server, log back into
the same linux server again. Then use the following command to check your
running processes. 
```bash
$ htop -u your_gfz_user_name
```
> NOTE: Press `q` to exit the interface of the command `htop`