#!/usr/bin/env bash

# An example workflow to estimate gap-filled z0 and d and then estimate flux
# footprints.
# 
# Zhan Li, zhan.li@gfz-potsdam.de

# The character '\' in the bash shell script is just a line continuatioin
# character, make sure no space behind '\'

# Step 1 
echo "Step 1: estimate z0d using multiple approaches in Graf et al. 2014"
python 01a_estimate_z0d_ensemble.py ./tests/data/test_01a_estimate_z0d_ensemble_config.ini

# Step 2
echo "Step 2: aggregate multiple estimates of z0d per half-hourly record"
python 02a_aggregate_z0d_ensemble.py \
    --cols Timestamp wind_dir zm \
    -- \
    ./tests/data/test_01a_estimate_z0d_ensemble_output_by_user_run.csv \
    ./tests/data/test_02a_aggregate_z0d_ensemble_output_by_user_run.csv

# Step 3
echo "Step 3: fill gaps in z0d using wind direction and time in a K-Nearest-Neighbor algorithm"
python 03a_gapfill_z0d_knn.py ./tests/data/test_03a_gapfill_z0d_knn_config.ini

# Step 4
echo """Before starting Step 4, please prepare a CSV file like
  ./tests/data/test_04a_gen_footprints_netcdf_input.csv that combines the
  gap-filled ensemble estimates of z0 and d and the other required
  micrometeorological data."""

echo "Do you wish to continue to Step 4? Type 1 for Yes and 2 for No"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;;
    esac
done

echo "Step 4: estimate flux footprints and save into a NetCDF using gap-filled z0d and micrometeorological data"
bash 04a_gen_footprints_netcdf.sh \
    --grid_srs="epsg:32633" \
    --grid_bounds="360725.0 361725.0 5971285.0 5972285.0" \
    --receptor_xy="361224.023952403 5971784.23062857" \
    --run_res=0.1 \
    --out_res=1 \
    --njobs=8 \
    ./tests/data/test_04a_gen_footprints_netcdf_input.csv \
    ./tests/data/test_04a_gen_footprints_netcdf_output_by_user_run.nc