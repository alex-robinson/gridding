
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
    my.write.table(poly_grl,"polygon_grl_Bamber2013.txt")

    coords <- locator(type="l")
    poly_grl_inner = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_grl_inner,"polygon_grl-inner_Bamber2013.txt")

    coords <- locator(type="l")
    poly_iceland = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_iceland,"polygon_iceland_Bamber2013.txt")

    coords <- locator(type="l")
    poly_ellesmere = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_ellesmere,"polygon_ellesmere_Bamber2013.txt")

    coords <- locator(type="l")
    poly_svalbard = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_svalbard,"polygon_svalbard_Bamber2013.txt")
} else {
    poly_grl       = read.table("polygon_grl_Bamber2013.txt",header=TRUE)
    poly_grl_inner = read.table("polygon_grl-inner_Bamber2013.txt",header=TRUE)
    poly_iceland   = read.table("polygon_iceland_Bamber2013.txt",header=TRUE)
    poly_ellesmere = read.table("polygon_ellesmere_Bamber2013.txt",header=TRUE)
    poly_svalbard = read.table("polygon_svalbard_Bamber2013.txt",header=TRUE)
       
}

points(poly_grl,col=1,pch=20,cex=1.2)
points(poly_grl_inner,col="grey40",pch=20,cex=1.2)

points(poly_iceland,col=1,pch=20,cex=1.2)
points(poly_ellesmere,col=1,pch=20,cex=1.2)
points(poly_svalbaard,col=1,pch=20,cex=1.2)

