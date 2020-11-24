###### Badgeley et al., 2020 -> add ICE-5G elevations to the output #######

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(fields)

domain="Greenland"
res="GRL-16KM"

## B20 data output to be modified
out.fldr=paste("~/apps/gridding_new/output/",domain,"/",res,"/B20_sig50km",sep="")
filename=paste(out.fldr,"/",res,"_CLIM-RECON-B20.nc",sep="")
nc=nc_open(filename)
time.b20=ncvar_get(nc, "time")
xc.b20=ncvar_get(nc, "xc")
yc.b20=ncvar_get(nc, "yc")
nc_close(nc)

## data from ICE-5G
file.topo="~/data/ice_data_yelmo/Greenland/GRL-16KM/GRL-16KM_TOPO-ICE-5G.nc"
nc=nc_open(file.topo)
zs.topo=ncvar_get(nc, "zs")
time.topo=ncvar_get(nc, "time")
nc_close(nc)

# additional modifications
zs.topo[zs.topo<0]=0
time.topo=time.topo*1000

# take only zs between -20000 and 0
zs.topo=zs.topo[,,2:39]
time.topo=time.topo[2:39]


### Interpolate surface elevations over time to get surfaces every 50 years

# define the final array which will contain surface elevations interpolated every 50 years between -20 kyr and 0 kyr.
zs.b20=array(NA,dim=c(dim(zs.topo)[1],dim(zs.topo)[2], length(time.b20)))
# index to increase the count of B20 snapshots given the B20 temporal resolution
k=0

for(j in 1:(length(time.topo)-1)){
  
  # interval of time at which we want to interpolate 
  t1=time.topo[j]
  t2=time.topo[j+1]
  
  # zs at time j and j+1
  z1=zs.topo[,,j]
  z2=zs.topo[,,j+1]
  
  # temporary timeseries (1000 yr, at 50 yr res) at which we want to interpolate
  dt=50
  t=seq(t1,t2,dt)
  
  w = (t - t1)/(t2-t1)
  
  # surface interpolation at each time t
  z=array(NA,dim=c(dim(z1)[1],dim(z1)[2], length(t)))
  for(i in 1:length(t)){
    z[,,i]=z1*(1-w[i]) + z2*w[i]
    zs.b20[,,i+k]=z[,,i]  
  }
  
  # increase index k for final B20 snapshots
  k=k+length(t)-1
}


### Include zs.b20 in B20 climatology file
nc=nc_open(paste(out.fldr,"/",res,"_CLIM-RECON-B20.nc",sep=""), write=T)

# Get dimensions
xc <- nc$dim[["xc"]]
yc <- nc$dim[["yc"]]
time <- nc$dim[["time"]]

# Define variables
varzs.def=ncvar_def("zs",units="m", dim=list(xc,yc,time), missval=-9999, longname="Ice elevation", prec="float")

# Add the new variable
nc = ncvar_add(nc, varzs.def)

# # Put the new variable
ncvar_put(nc, varzs.def, zs.b20)

# Close the nc file
nc_close(nc)