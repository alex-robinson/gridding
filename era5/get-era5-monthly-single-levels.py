#!/usr/bin/env python

import cdsapi
import os, sys, json  

##### User options ###########

# Load from config file
info      = json.load(open("era5_config_monthly.json"))
name      = info["name"]
area      = None
year0     = info["year0"]
year1     = info["year1"]

years = list(range(year0,year1+1))

#years = list(map(str, years_num))
#year0 = years[0]    # First year
#year1 = years[-1]   # Last  year

vars      = info["vars"]

# Make output folder with the name the same as the region
fldr = "data/era5/"+name 
if not os.path.exists(fldr):
    os.makedirs(fldr)

print("Summary:")
print("name: {name}".format(name=name))
print("years: {years}".format(years=years))
print(vars)

c = cdsapi.Client()

# Do not loop over years. Instead this script will be called 
# once per variable per year, to get request started on CDS server. 
args        = str(sys.argv) 
var_now     = sys.argv[1]
year_now    = sys.argv[2]

if len(sys.argv)==4:
    plev_now = sys.argv[3]
else:
    plev_now = None 

# If var_now is an integer, select the variable name of interest 
try:
    var_int = int(var_now)
    if var_int > len(vars)-1: raise
    var_now = vars[var_int]
except:
    # Pass, assume the variable name was provided
    print("\nOnly 0-{n} variable names available. Try again.\n".format(n=len(vars)-1))
    sys.exit(2)

if not var_now in vars:
    print("\nVariable not available: {var}\n".format(var=var_now))
    sys.exit(2)

print("var  = {var}".format(var=var_now))
print("plev = {plev}".format(plev=plev_now))
print("year = {year}".format(year=year_now))

#sys.exit() 

if plev is None:
    dataset     = "reanalysis-era5-single-levels-monthly-means"
    plev_str    = ""
else:
    dataset     = "reanalysis-era5-pressure-levels-monthly-means"
    plev_str    = "_{plev}".format(plev=plev)

# ajr: disable loop over years 
# for year_now in years:

year_now_str = [str(year_now)]

# Define output filename for this year
filename_download = "{fldr}/era5_{name}_{var}{plev}_{year}.nc".format(
                                fldr=fldr,name=name,var=var_now,plev=plev_str,year=year_now)


### BEGIN DOWNLOAD REQUEST ###

c.retrieve(
    dataset,
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': var_now,
        'pressure_level': plev,
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        'area': area,
        'year': year_now_str,
    },
    filename_download)



