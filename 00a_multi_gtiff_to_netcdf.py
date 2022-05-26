#!/usr/bin/env python

# Generate a NetCDF from a list of Geotiff files and their timestamps.

import argparse

import numpy as np
import pandas as pd

import xarray as xr
import rioxarray as rxr

def getCmdArgs():
    p = argparse.ArgumentParser(description='Convert multiple GeoTiff files with their bands named after timestamps into a NetCDF.')
    p.add_argument('gtiff_list')
    p.add_argument('output_nc')

    cmdargs = p.parse_args()
    return cmdargs

def clean_my_xda(my_xda):
    my_xda = my_xda.loc[{'band':1}].drop_vars('band')
    my_xda.name = 'footprint_' + my_xda.attrs.pop('long_name')
    return my_xda

def main(cmdargs):
    gtiff_list_file = cmdargs.gtiff_list
    output_nc = cmdargs.output_nc

    file_list = pd.read_csv(gtiff_list_file, header=None).iloc[:, 0].to_list()

    xda_list = [rxr.open_rasterio(val, chunks=1024) for val in file_list]

    dt_list = [pd.Timestamp(val.attrs['long_name']) for val in xda_list]
    sort_idx = np.argsort(dt_list)

    xda_list = [xda_list[val] for val in sort_idx]
    dt_list = [dt_list[val] for val in sort_idx]
    xda_list = [clean_my_xda(val) for val in xda_list]

#     # Combine all the DataArrays along the timestamp dimension into one DataArray
#     dt_idx = pd.DatetimeIndex(dt_list)
#     dt_idx.name = 'Timestamp'
#     xda = xr.concat(xda_list, dim=dt_idx)
#     xds = xda.to_dataset(name='footprints')

    # Combine all the DataArrays into separate variables in a Dataset
    xds = xr.merge(xda_list, join='exact')

    # xds = xds.rio.write_crs(xds.rio.crs)
    # xds = xds.rio.write_transform()
    xds = xds.rio.write_coordinate_system()
    # Not sure why rio.write_grid_mapping does not add the attribute
    # 'grid_mapping'. It worked in rioxarray 0.3.1 but not in rioxarray 0.9.0
    xds = xds.rio.write_grid_mapping()
    comp = dict(zlib=True, complevel=9)
    encoding = {var: comp for var in xds.data_vars}
    xds.to_netcdf(output_nc, format='netcdf4',  
        encoding=encoding)

if __name__ == '__main__':
    cmdargs = getCmdArgs()
    main(cmdargs)