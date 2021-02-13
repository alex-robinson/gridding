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

# Make output folder with the name the same as the region
fldr = "data/ERA5/"+name 
if not os.path.exists(fldr):
    os.makedirs(fldr)

print("Summary:")
print("name: {name}".format(name=name))
print("years: {years}".format(years=years))

vars = [
            '2m_temperature', 'sea_surface_temperature', 'total_precipitation',
        ]
print(vars)

c = cdsapi.Client()

# Do not loop over years. Instead this script will be called 
# once per each year, to get request started on CDS server. 
args = str(sys.argv) 
year_now = sys.argv[1]
download_rhum = False

print("year = {year}".format(year=year_now))

# for year_now in years:

year_now_str = [str(year_now)]

# Define output filename for this year
filename_download = "{fldr}/era5_{name}_{year}.nc".format(fldr=fldr,name=name,year=year_now)
filename_download_rhum = "{fldr}/era5_{name}_rhum_{year}.nc".format(fldr=fldr,name=name,year=year_now)

if not download_rhum:

    # Single pressure-level variables
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
            'year': year_now_str,
        },
        filename_download)

else:

    # Relative humidity at 1000 hPa pressure level
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': 'relative_humidity',
            'pressure_level': '1000',
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
        filename_download_rhum)
