#!/home/robinson/apps/R/R/bin/Rscript
# #--vanilla
#args <- commandArgs(TRUE)

## Load functions
source("functions.r")
source("interpolator.r")

require(akima)

## First, if we want a new file, create it using ncdf
## and grid creation functions
create_regional_file <- function(filename,dx=100,missing=9999.0,nt=540,
                                 month0=9,tunits="months since 1957-08-01") {

  x0 = seq( -900.0, 700.0,by=dx)
  y0 = seq(-3500.0,-500.0,by=dx)
  
  # grid = creategrid(-800.0,nx,-3400.0,ny,dx=dx)
  grid = creategrid(x=x0,y=y0)

  # Define dimensions
  x = dim.def.ncdf( "xc", "km", grid$x)
  y = dim.def.ncdf( "yc", "km", grid$y)
  t = dim.def.ncdf( "time", tunits, c(1:nt), unlim=TRUE)
  
  # Make a month variable to ensure the "time" dimension is saved in file
  months0 = c(month0:12)
  nm0 = length(months0)
  nm1 = floor((nt - nm0)/12)
  months = c(months0, rep(c(1:12),nm1))
  nm2 = length(months)
  nm3 = nt - nm2
  if (nm3 > 0) {
    months = c(months, c(1:nm3))
  }

  mon = var.def.ncdf( name="mon",units="",dim=list(t),missval=missing,longname="Months" )

  # Define 2D fields
  xx  = var.def.ncdf( name="xx",units="km",dim=list(x,y),missval=missing,longname="Stereographic x" )
  yy  = var.def.ncdf( name="yy",units="km",dim=list(x,y),missval=missing,longname="Stereographic y" )
  lon = var.def.ncdf( name="lon",units="degrees East",dim=list(x,y),missval=missing,longname="Longitude" )
  lat = var.def.ncdf( name="lat",units="degrees North",dim=list(x,y),missval=missing,longname="Latitude" )
  
  # Now actually create the netCDF file with these variables
  ncvars = list(mon,xx,yy,lon,lat)
  ncnew <- create.ncdf( filename, ncvars )
  
  put.var.ncdf( ncnew, mon, months )
  put.var.ncdf( ncnew, xx,  grid$xx )
  put.var.ncdf( ncnew, yy,  grid$yy )
  put.var.ncdf( ncnew, lon, grid$lon )
  put.var.ncdf( ncnew, lat, grid$lat )

  close.ncdf(ncnew)
  
  return
}

var_to_grid <- function(infile,invarname,outfile,outvarname=invarname,dimnames=c("xc","yc","time"),
                        units="",longname=outvarname,missing=9999.0,transform=NULL)
{
  cat("var_to_grid:",invarname,"=>",outvarname,"\n")

  # First load grid information from output file
  out = list()
  nc = open.ncdf(outfile)
  out$x    = get.var.ncdf(nc,"xc")
  out$y    = get.var.ncdf(nc,"yc")
  out$time = get.var.ncdf(nc,"time")
  out$lon = get.var.ncdf(nc,"lon")
  out$lat = get.var.ncdf(nc,"lat")
  close.ncdf(nc)

  # Determine the range of the output grid
  lon.range = range(out$lon)
  lat.range = range(out$lat)
  
  # Length of time series expected
  nt = length(out$time)

  # Now load variable of interest from input file
  indat = list()
  nc = open.ncdf(infile)
  tmp = names(nc$dim)
  lonnm = "lon"
  if ( !"lon" %in% tmp ) lonnm = "longitude"
  latnm = "lat"
  if ( !"lat" %in% tmp ) latnm = "latitude"

  indat$lon  = get.var.ncdf(nc,lonnm)
  indat$lat  = get.var.ncdf(nc,latnm)
  close.ncdf(nc)

  ## Corrections to dimensions and input variable ##
  
  # Make sure longitude is -180:180
  if ( range(indat$lon)[2] > 180 ) {
    ii = which(indat$lon > 180)
    indat$lon[ii] = indat$lon[ii] - 360.0
  }

  # Load variable
  nc = open.ncdf(infile)
  indat$var = get.var.ncdf(nc,invarname)
  close.ncdf(nc)

  # Make sure var is 3d (including time dimension!)
  if (length(dim(indat$var))==2) dim(indat$var) = c(dim(indat$var),1)
  
  ## End corrections

  points0 = expand.grid(lon=indat$lon,lat=indat$lat)
  points0$x = points0$y = NA 
  
  ebs = 1.0  # 1 degree buffer range for points of interest
  ij = which(points0$lon >= lon.range[1]-ebs & points0$lon <= lon.range[2]+ebs &
             points0$lat >= lat.range[1]-ebs & points0$lat <= lat.range[2]+ebs )
  
  # For points of interest get x,y of stereographic projection
  cat("Generating stereographic projection points.\n")
  # for (q in ij ) {
  #   xy = plane_ellipsoid(coord=c(points0$lon[q],points0$lat[q]))
  #   points0$x[q] = xy$x 
  #   points0$y[q] = xy$y   
  # }
  
  xy = proj_plane(points0$lon[ij],points0$lat[ij])
  points0$x[ij] = xy$x
  points0$y[ij] = xy$y

  # Now perform interpolation using akima package
  cat("Performing interpolation via akima package 'interp' function.\n")
  nt = min(dim(indat$var)[3],nt)
  out$var = array(NA,dim=c(length(out$x),length(out$y),nt))
  for (k in 1:nt) {
    varnow = indat$var[,,k]
    tmp = interp(x=points0$x[ij],y=points0$y[ij],z=varnow[ij],xo=out$x,yo=out$y)
    out$var[,,k] = tmp$z
  }

  ## Apply the transformation function if exists
  if (!is.null(transform)) out$var = transform(out$var)
  
  ## If data is a mask, clean it up
  if (outvarname == "mask") {
    out$var = round(out$var)
  }

  ## Now that new variable is ready, save it to netcdf
  nc = open.ncdf(outfile)
  ncdims = nc$dim[dimnames]
  varnames = names(nc$var)
  
  # Define new variable or extract var info if exists
  if ( ! outvarname %in% varnames ) {
    varout = var.def.ncdf( name=outvarname,units=units,dim=ncdims,
                         missval=missing,longname=longname )
    nc = open.ncdf(outfile,write=TRUE)
    nc = var.add.ncdf(nc,varout)
    close.ncdf(nc)
  } else {
    varout = nc$var[[outvarname]]
    close.ncdf(nc)
  }
  

  # Put the variable
  nc = open.ncdf(outfile,write=TRUE)
  att.put.ncdf( nc, outvarname,"units",units)
  put.var.ncdf( nc, varout, out$var )
  close.ncdf(nc)
  ##

  return(out)
}

######################

#### START SCRIPT ####

######################

## GET ERA40 data onto REMBO 100km grid
if (FALSE) {
  
  outfile = "output/ERA-40-GRL100km_historical_mon_195709-200108.nc"
  create_rembo_input(outfile,nx=16,ny=29,dx=100,missing=9999.0,nt=540,
                     month0=9,tunits="months since 1957-08-31")

  file1 = "data/ECMWF/era.40.invariant.nc"
  tmp = var_to_grid(infile=file1,invarname="z",outfile=outfile,outvarname="zs",
                    units="m",longname="Surface elevation",dimnames=c("x","y"),transform=geo2zs)
  tmp = var_to_grid(infile=file1,invarname="lsm",outfile=outfile,outvarname="mask",
                    dimnames=c("x","y") )
  
  file2 = "data/ECMWF/era.40.monthly.nc"
  tmp = var_to_grid(infile=file2,invarname="r",outfile=outfile,outvarname="rhum",
                    longname="Relative humidity",dimnames=c("x","y","time") )
  
  file3 = "data/ECMWF/era.40.monthly.surface.nc"
  tmp = var_to_grid(infile=file3,invarname="t2m",outfile=outfile,outvarname="t2m",
                    units="K",longname="2m air temperature",dimnames=c("x","y","time") )
  
  file4 = "data/ECMWF/era.40.monthly.850.winduv.nc"
  tmp = var_to_grid(infile=file4,invarname="u",outfile=outfile,outvarname="u850",
                    longname="Wind speed - x",dimnames=c("x","y","time") )
  tmp = var_to_grid(infile=file4,invarname="v",outfile=outfile,outvarname="v850",
                    longname="Wind speed - y",dimnames=c("x","y","time") )

  file5 = "data/ECMWF/era.40.monthly.775.nc"
  tmp = var_to_grid(infile=file5,invarname="u",outfile=outfile,outvarname="u750",
                    longname="Wind speed - x",dimnames=c("x","y","time") )
  tmp = var_to_grid(infile=file5,invarname="v",outfile=outfile,outvarname="v750",
                    longname="Wind speed - y",dimnames=c("x","y","time") )
  tmp = var_to_grid(infile=file5,invarname="w",outfile=outfile,outvarname="w750",
                    longname="Wind speed - w",dimnames=c("x","y","time") )
  tmp = var_to_grid(infile=file5,invarname="t",outfile=outfile,outvarname="t750",
                    longname="Temperature 750 Mb",dimnames=c("x","y","time") )
  tmp = var_to_grid(infile=file5,invarname="z",outfile=outfile,outvarname="z750",
                    longname="Geopotential 750 Mb",dimnames=c("x","y","time") )
  tmp = var_to_grid(infile=file5,invarname="r",outfile=outfile,outvarname="rhum750",
                    longname="Relative humidity 750 Mb",dimnames=c("x","y","time") )

}

## GET ERAINTERIM data onto REMBO 100km grid
if (TRUE) {
  
  outfile = "output/ERA-INTERIM-GRL100KM_historical_mon_197901-201212.nc"
  create_regional_file(outfile,dx=50,missing=9999.0,nt=408,
                       month0=1,tunits="months since 1978-12-31")

  file1 = "data/ECMWF/NEW/ERA-INTERIM-invariant_historical_mon_197901-201212.nc"
  tmp = var_to_grid(infile=file1,invarname="z",outfile=outfile,outvarname="zs",
                    units="m",longname="Surface elevation",dimnames=c("xc","yc"),transform=geo2zs)
  tmp = var_to_grid(infile=file1,invarname="lsm",outfile=outfile,outvarname="mask",
                    dimnames=c("xc","yc") )
  
  file2 = "data/ECMWF/NEW/ERA-INTERIM-surface_historical_mon_197901-201212.nc"
  tmp = var_to_grid(infile=file2,invarname="sp",outfile=outfile,outvarname="sp",
                    units="Pa",longname="Surface pressure",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file2,invarname="tcw",outfile=outfile,outvarname="tcw",
                    units="kg m**-2",longname="Total column water",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file2,invarname="msl",outfile=outfile,outvarname="msl",
                    units="Pa",longname="Mean sea level pressure",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file2,invarname="tcc",outfile=outfile,outvarname="cc",
                    units="1",longname="Total cloud cover",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file2,invarname="v10u",outfile=outfile,outvarname="v10u",
                    units="m s**-1",longname="10-meter velocity, u-component",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file2,invarname="v10v",outfile=outfile,outvarname="v10v",
                    units="m s**-1",longname="10-meter velocity, v-component",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file2,invarname="v2t",outfile=outfile,outvarname="v2t",
                    units="K",longname="2-meter temperature",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file2,invarname="al",outfile=outfile,outvarname="al",
                    units="1",longname="Surface albedo",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file2,invarname="sst",outfile=outfile,outvarname="sst",
                    units="K",longname="Sea surface temperature",dimnames=c("xc","yc","time") )
  

  file3 = "data/ECMWF/NEW/ERA-INTERIM-750Mb_historical_mon_197901-201212.nc"
  tmp = var_to_grid(infile=file3,invarname="z",outfile=outfile,outvarname="v750z",
                    units="m**2 s**-2",longname="750Mb geopotential",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file3,invarname="t",outfile=outfile,outvarname="v750t",
                    units="K",longname="750Mb temperature",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file3,invarname="q",outfile=outfile,outvarname="v750q",
                    units="kg kg**-1",longname="750Mb specific humidity",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file3,invarname="w",outfile=outfile,outvarname="v750w",
                    units="Pa s**-1",longname="750Mb vertical velocity",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file3,invarname="r",outfile=outfile,outvarname="v750r",
                    units="%",longname="750Mb Relative humidity",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file3,invarname="cc",outfile=outfile,outvarname="v750cc",
                    units="1",longname="750Mb cloud cover",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file3,invarname="u",outfile=outfile,outvarname="v750u",
                    units="m s**-1",longname="750Mb velocity, u-component",dimnames=c("xc","yc","time") )
  tmp = var_to_grid(infile=file3,invarname="v",outfile=outfile,outvarname="v750v",
                    units="m s**-1",longname="750Mb velocity, v-component",dimnames=c("xc","yc","time") )
  
}
