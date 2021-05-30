import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': 'orography',
        'year': '2020',
        'month': '06',
        'day': '30',
        'time': '00:00',
    },
    'data/era5/era5_orography.nc')