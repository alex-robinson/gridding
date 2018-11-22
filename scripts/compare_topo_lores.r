library(myr)
library(RNetCDF)

# Load data 

topo_con = my.read.nc("../output/Greenland/GRL-20KM/GRL-20KM_TOPO-RTOPO-2.0.1.nc")
topo_gau = my.read.nc("../output/Greenland/GRL-20KM/GRL-20KM_TOPO-RTOPO-2.0.1__gaussian.nc")

image.plot(topo_con$xc,topo_con$yc,topo_con$z_bed-topo_gau$z_bed)
