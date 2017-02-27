
library(RNetCDF)
library(myr)

a = open.nc("../ice_data/Antarctica/ANT-40KM/ANT-40KM_TOPO-BEDMAP2.nc")
topo = read.nc(a)
close.nc(a)

image.plot(topo$xc,topo$yc,
            topo$zb,zlim=c(-2500,0))

contour(topo$xc,topo$yc,
        topo$zb,levels=c(-2000,-1000,-500),lwd=c(2,2,1),col=c("grey40","black","grey80"),add=TRUE)

if (FALSE) {
    coords <- locator(type="l")
    poly_ant = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_ant,"BEDMAP2/polygon_ant_BEDMAP2.txt")

    coords <- locator(type="l")
    poly_ant_inner = data.frame(x=coords$x,y=coords$y)
    my.write.table(poly_ant_inner,"BEDMAP2/polygon_ant-inner_BEDMAP2.txt")

} else {
    poly_ant       = read.table("BEDMAP2/polygon_ant_BEDMAP2.txt",header=TRUE)
    poly_ant_inner = read.table("BEDMAP2/polygon_ant-inner_BEDMAP2.txt",header=TRUE)
   
}

if (TRUE) {
    points(poly_ant,col=1,pch=20,cex=1.2)
    points(poly_ant_inner,col="grey40",pch=20,cex=1.2)

}

