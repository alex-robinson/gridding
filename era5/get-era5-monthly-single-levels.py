import cdsapi
import os, sys, json  

##### User options ###########

# Load from config file
info      = json.load(open("era5_config_monthly.json"))
name      = info["name"]
area      = info["area"]
years_num = info["years"]

years = list(map(str, years_num))
print(years)

year0 = years[0]    # First year
year1 = years[-1]   # Last  year

# Make output folder with the name the same as the region
fldr = "data/era5/"+name 
if not os.path.exists(fldr):
    os.makedirs(fldr)

# Define output filename
if year0 == year1: 
    # Only Downloading one year 
    filename_download = "{fldr}/era5_{name}_{year0}.nc".format(fldr=fldr,name=name,year0=year0)
else:
    # Downloading a range of years
    filename_download = "{fldr}/era5_{name}_{year0}-{year1}.nc".format(fldr=fldr,name=name,year0=year0,year1=year1)


print("Summary:")
print("name: {name}".format(name=name))
print("filename_download: {fnm}".format(fnm=filename_download))
#sys.exit()


vars = [
            '2m_temperature', 'sea_surface_temperature', 'total_precipitation',
        ]
print(vars)

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': vars,
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        'area': area,
        'year': years,
    },
    filename_download)