import cdsapi
import os, sys, json  

##### User options ###########

# Load from config file
info      = json.load(open("era5_config_dubai.json"))
name      = info["name"]
area      = info["area"]
year0     = info["year0"]
year1     = info["year1"]

years = list(range(year0,year1+1))

# Make output folder with the name the same as the region
fldr = "data/era5/"+name  
if not os.path.exists(fldr):
    os.makedirs(fldr)

print("Summary:")
print("name: {name}".format(name=name))
print("area: {area}".format(area=area))
print("years: {years}".format(years=years))

#sys.exit()

############################## 

vars_all = ['2m_temperature',
            'soil_temperature_level_1', 
            'soil_temperature_level_2', 
            'soil_temperature_level_3',
            'soil_temperature_level_4', 
            'toa_incident_solar_radiation',
            'surface_solar_radiation_downwards', 
            'total_sky_direct_solar_radiation_at_surface',
            'clear_sky_direct_solar_radiation_at_surface',
            'total_cloud_cover',
            'surface_pressure']

            #'total_precipitation',
            #'10m_u_component_of_wind', 
            #'10m_v_component_of_wind'

vars = vars_all 
print(vars)

c = cdsapi.Client()


# Do not loop over years. Instead this script will be called 
# once per each year, to get request started on CDS server. 
args = str(sys.argv) 
year_now = sys.argv[1]

print("year = {year}".format(year))

#for year_now in years:

year_now_str = [str(year_now)]

# Define output filename for this year
filename_download = "{fldr}/era5_{name}_{year}.nc".format(fldr=fldr,name=name,year=year_now)

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': vars,
        'year': year_now_str,
        'month': [
            '01', '02', '03', '04', '05', '06', 
            '07', '08', '09', '10', '11', '12',
        ],
        'day': [
            '01', '02', '03','04', '05', '06',
            '07', '08', '09','10', '11', '12',
            '13', '14', '15','16', '17', '18',
            '19', '20', '21','22', '23', '24',
            '25', '26', '27','28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': area,
    },
    filename_download)


