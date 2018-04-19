# joughin
library(raster)
library(fields)
library(ncdf4)

# *********************************************************************************************
# vx
vx.file="~/Documentos/data/Joughin2017_vel/greenland_vel_mosaic250_vx_v1.tif"
vx=raster(vx.file)
vx = aggregate(vx, factor=4)

# create a raster with lower resolution (2 km)
a=vx
res(a)=c(2000,2000)
vx=resample(vx, a)
vx[is.na(vx)==T]=0

# Convert raster to SpatialPointsDataFrame
r.pts <- rasterToPoints(vx, spatial=TRUE)

# reproject sp object
geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
r.pts <- spTransform(r.pts, CRS(geo.prj))

# Assign coordinates to @data slot, display first 6 rows of data.frame
r.pts@data <- data.frame(r.pts@data, long=coordinates(r.pts)[,1],
                         lat=coordinates(r.pts)[,2])
# *********************************************************************************************
# vy
vy.file="~/Documentos/data/Joughin2017_vel/greenland_vel_mosaic250_vy_v1.tif"
vy=raster(vy.file)
vy = aggregate(vy, factor=4)

# create a raster with lower resolution
a=vy
res(a)=c(2000,2000)
vy=resample(vy, a)
vy[is.na(vy==T)]=0

# Convert raster to SpatialPointsDataFrame
r.ptsy <- rasterToPoints(vy, spatial=TRUE)

# reproject sp object
geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
r.ptsy <- spTransform(r.ptsy, CRS(geo.prj))

# Assign coordinates to @data slot, display first 6 rows of data.frame
r.ptsy@data <- data.frame(r.ptsy@data, long=coordinates(r.ptsy)[,1],
                         lat=coordinates(r.ptsy)[,2])

r.pts@data = data.frame(vx = r.pts@data[,1], vy = r.ptsy@data[,1], long = r.pts@data[,2], lat = r.pts@data[,3])
r.pts@data = round(r.pts@data, digit=6)

# *********************************************************************************************************************+
# write the final data in a .txt file
write.table(r.pts@data,"~/Documentos/GRISLI-UCM/v0.32/skill-score/lhs/LEs_stone/250.sims/greenland_vel_mosaic250_vx_v1.txt",
            row.names=F,col.names=FALSE,quote=FALSE)

