#!/usr/bin/env python


import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
        {
         'format': 'netcdf',
         'variable': 'geopotential',
         'product_type': 'monthly_averaged_reanalysis',
         'year': '2022',
         'month': '01',
         'time': '00:00',
        },
     'data/era5/era5_geopotential.nc')
