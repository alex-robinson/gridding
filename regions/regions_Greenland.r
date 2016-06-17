
library(RNetCDF)
library(myr)

a = open.nc("../data/Greenland/Greenland_bedrock_topography_V3.nc")
topo = read.nc(a)
close.nc(a)

image.plot(topo$projection_x_coordinate,
            topo$projection_y_coordinate,
            topo$BedrockElevation,zlim=c(-2000,0))

if (FALSE) {
    coords <- locator(type="l")
    poly_grl = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_grl,"poly_grl.txt")

    coords <- locator(type="l")
    poly_grl_inner = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_grl_inner,"poly_grl_inner.txt")

    coords <- locator(type="l")
    poly_iceland = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_iceland,"poly_iceland.txt")

    coords <- locator(type="l")
    poly_ellesmere = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_ellesmere,"poly_ellesmere.txt")

    coords <- locator(type="l")
    poly_svalbaard = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_svalbaard,"poly_svalbaard.txt")
} else {
    poly_grl       = read.table("poly_grl.txt",header=TRUE)
    poly_grl_inner = read.table("poly_grl_inner.txt",header=TRUE)
    poly_iceland   = read.table("poly_iceland.txt",header=TRUE)
    poly_ellesmere = read.table("poly_ellesmere.txt",header=TRUE)
    poly_svalbaard = read.table("poly_svalbaard.txt",header=TRUE)
       
}

points(poly_grl,col=1,pch=20,cex=1.2)
points(poly_grl_inner,col="grey40",pch=20,cex=1.2)

points(poly_iceland,col=1,pch=20,cex=1.2)
points(poly_ellesmere,col=1,pch=20,cex=1.2)
points(poly_svalbaard,col=1,pch=20,cex=1.2)

