
library(RNetCDF)
library(myr)

a = open.nc("../data/Greenland/Greenland_bedrock_topography_V3.nc")
topo = read.nc(a)
close.nc(a)

image.plot(topo$projection_x_coordinate,topo$projection_y_coordinate,
            topo$BedrockElevation,zlim=c(-2500,0))

contour(topo$projection_x_coordinate,topo$projection_y_coordinate,
        topo$BedrockElevation,levels=c(-2000,-500,-400),lwd=c(2,2,1),col=c("grey40","black","grey80"),add=TRUE)

if (FALSE) {
    coords <- locator(type="l")
    poly_grl = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_grl,"Bamber2013/polygon_grl_Bamber2013.txt")

    coords <- locator(type="l")
    poly_grl_inner = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_grl_inner,"Bamber2013/polygon_grl-inner_Bamber2013.txt")

    coords <- locator(type="l")
    poly_iceland = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_iceland,"Bamber2013/polygon_iceland_Bamber2013.txt")

    coords <- locator(type="l")
    poly_ellesmere = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_ellesmere,"Bamber2013/polygon_ellesmere_Bamber2013.txt")

    coords <- locator(type="l")
    poly_svalbard = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_svalbard,"Bamber2013/polygon_svalbard_Bamber2013.txt")

} else {
    poly_grl       = read.table("Bamber2013/polygon_grl_Bamber2013.txt",header=TRUE)
    poly_grl_inner = read.table("Bamber2013/polygon_grl-inner_Bamber2013.txt",header=TRUE)
    poly_iceland   = read.table("Bamber2013/polygon_iceland_Bamber2013.txt",header=TRUE)
    poly_ellesmere = read.table("Bamber2013/polygon_ellesmere_Bamber2013.txt",header=TRUE)
    poly_svalbard  = read.table("Bamber2013/polygon_svalbard_Bamber2013.txt",header=TRUE)
       
}

if (FALSE) {
    points(poly_grl,col=1,pch=20,cex=1.2)
    points(poly_grl_inner,col="grey40",pch=20,cex=1.2)

    points(poly_iceland,col=1,pch=20,cex=1.2)
    points(poly_ellesmere,col=1,pch=20,cex=1.2)
    points(poly_svalbaard,col=1,pch=20,cex=1.2)
}

