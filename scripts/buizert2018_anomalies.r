###### Buizert 2018 -> modify it to be regridded #######

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(fields)

out.fldr=paste("~/apps/gridding/data/buizert2018")

filename=paste(out.fldr,"/GLand_22ka_recon_Buizert_20161228.nc",sep="")
nc=nc_open(filename)
lat2d=ncvar_get(nc, "LAT")
lon2d=ncvar_get(nc, "LON")
tas=ncvar_get(nc, "TS2")  # x, y, m, t
pr=ncvar_get(nc, "PRECIP") # x, y, m, t
time=ncvar_get(nc, "TIME")
yearCE=ncvar_get(nc, "YEAR_CE")
mon=ncvar_get(nc, "MONTH")
dem=ncvar_get(nc, "DEM")
nc_close(nc)

# calculate temperature anomalies (wrt present day)
dtas=tas
nx =dim(tas)[1]
ny=dim(tas)[2]
nt=length(time)
nm=length(mon)
for(x in 1:nx){
  for(y in 1:ny){
    for(m in 1:nm){
      for(t in 1:nt){
        dtas[x,y,m,t]=tas[x,y,m,t] - tas[x,y,m,nt]
      }
    }
  }
}

# calculate precipitation fraction (wrt present day)
prf=pr
nx =dim(pr)[1]
ny=dim(pr)[2]
nt=length(time)
nm=length(mon)
for(x in 1:nx){
  for(y in 1:ny){
    for(m in 1:nm){
      for(t in 1:nt){
        prf[x,y,m,t]=pr[x,y,m,t]/pr[x,y,m,nt]
      }
    }
  }
}

# Replicate DEM for all times
dem.3d=array(NA, dim=c(nx, ny, nt))
for(i in 1:nt){
  dem.3d[,,i]=dem
}
dem=dem.3d

######################################### CREATE A NEW NetCDF FILE  #########################################

# Dimensions
lon=lon2d[,1]
lat=lat2d[1,]
mon=mon
time=time

# Define dimensions
lon.def=ncdim_def("lon", units="degrees", vals=lon, unlim=FALSE, create_dimvar=TRUE, longname="Longitude")
lat.def=ncdim_def("lat", units="degrees", vals=lat, unlim=FALSE, create_dimvar=TRUE, longname="Latitude")
mon.def=ncdim_def("month", units="", vals=mon, unlim=FALSE, create_dimvar=TRUE, longname="Month")
time.def=ncdim_def("time", units="Calendar years before present (present = 1950 C.E.)", vals=time, unlim=FALSE, create_dimvar=TRUE, longname="Age")

# Define variables
vardem.def=ncvar_def("DEM",units="m above sea level", dim=list(lon.def,lat.def,time.def), missval=-9999, longname="Ice elevation", prec="float")
vardtas.def=ncvar_def("tas",units="K", dim=list(lon.def,lat.def,mon.def,time.def), missval=-9999, longname="Surface air temperature at 2m height (anomaly)", prec="float")
varprf.def=ncvar_def("pr",units="m/s water equivalent", dim=list(lon.def,lat.def,mon.def,time.def), missval=-9999, longname="Precipitation fraction", prec="float")
varyearCE.def=ncvar_def("year_CE",units="years Common Era", dim=list(time.def), missval=-9999, longname="Age wrt Common Era", prec="float")

vars <- list(vardem.def,vardtas.def,varprf.def,varyearCE.def)

# Create the file
outputfile = paste(out.fldr,"/GLand_22ka_recon_Buizert_20161228_anom.nc",sep="")
ncout <- nc_create(outputfile, vars,force_v4=TRUE)

# Put variables
ncvar_put(ncout, vardem.def, dem)
ncvar_put(ncout, vardtas.def, dtas)
ncvar_put(ncout, varprf.def, prf)
ncvar_put(ncout, varyearCE.def, yearCE)

# Put global attributes
ncatt_put(ncout,0,"title","Buizert et al., 2018, Greenland 22,000 Year Seasonal Temperature Reconstructions")
ncatt_put(ncout,0,"source","https://www.ncdc.noaa.gov/paleo/study/23430")
ncatt_put(ncout,0,"references","Buizert et al., 2018, Greenland-Wide Seasonal Temperatures During the Last Deglaciation, GRL")

nc_close(ncout)

