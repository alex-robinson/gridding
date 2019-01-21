
library(RNetCDF)
library(myr)

plot_to_file = FALSE 

domain    = "Antarctica"
grid_name = "ANT-5KM"

a = open.nc(paste0("~/models/EURICE/ice_data/",domain,"/",grid_name,"/",grid_name,"_TOPO-RTOPO-2.0.1.nc"))
topo = read.nc(a)
close.nc(a)

zlim = c(-3000,3000)
z_bed = topo$z_bed 
z_bed[z_bed < zlim[1]] = zlim[1]
z_bed[z_bed > zlim[2]] = zlim[2]

xlim = range(topo$xc)
ylim = range(topo$yc)
#xlim = c(-1500,4000)
#ylim = c(-4800,2000)

if (plot_to_file) 
    myfigure("./",paste0("regions_",domain),asp=1.1,pointsize=14,type="png")

par(plt=c(0.1,0.95,0.07,0.95),xaxs="i",yaxs="i")
plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
axis(1)
axis(2)

image(topo$xc,topo$yc,z_bed,add=TRUE,zlim=c(-3000,3000),col=jet.colors)
image(topo$xc,topo$yc,z_bed,add=TRUE,zlim=c(0,5000),col=alpha("white",50))
contour(topo$xc,topo$yc,z_bed,add=TRUE,drawlabels=FALSE,
            levels=c(-2000,-1000,-500,0),lwd=c(2,2,1,1),col=c("grey40","black","grey80","black"))

points(0,0,pch=3,col="white",lwd=2)

# Box containing ANT domain 
x0 = -3120.0
nx = 157 
dx = 40.0 
x1 = x0 + (nx-1)*dx 
y0 = -2920.0
ny = 147 
y1 = y0 + (ny-1)*dx 
antbox = data.frame(x0=x0,x1=x1,y0=y0,y1=y1)
rect(antbox$x0,antbox$y0,antbox$x1,antbox$y1,border=2,col=NA)

if (plot_to_file) {
    box()
    graphics.off() 
}


get_latlon = function(pts,dat)
{   # Given projected x/y values, find nearest lon/lat values in projection

    out     = data.frame(x=pts$x,y=pts$y)

    out$lon = NA 
    out$lat = NA 

    np = length(out$x) 

    for (q in 1:np) {
        x = out$x[q]
        y = out$y[q]

        i = which.min(abs(out$x[q]-dat$xc))
        j = which.min(abs(out$y[q]-dat$yc))
        
        out$lon[q] = dat$lon2D[i,j]
        out$lat[q] = dat$lat2D[i,j]          
    } 

    return(out)
}

if (FALSE) {

    fldr = "polygons"

    ### Continental domains ###

    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_antarctica.txt"),scientific=TRUE,digits=4)

    ### Sub-domains ###

    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_antarctica_inner.txt"),scientific=TRUE,digits=4)
    
    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_wais.txt"),scientific=TRUE,digits=4)
    
    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_eais.txt"),scientific=TRUE,digits=4)
    
}

if (FALSE) {
    points(poly1,col=1,pch=20,cex=1.2)

}

