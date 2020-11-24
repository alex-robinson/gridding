
library(RNetCDF)
library(myr)

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

get_xy = function(pts,dat)
{   # Given lon/lat values, find nearest projected x/y values in projection

    out   = data.frame(lon=pts$lon,lat=pts$lat)
    out$x = NA 
    out$y = NA 

    np = length(out$x) 

    for (q in 1:np) {
        lon = out$lon[q]
        lat = out$lat[q]

        dist =  sqrt((lon-dat$lon2D)^2 + (lat-dat$lat2D)^2)
        ij   = which(dist == min(dist),arr.ind=TRUE)
        i    = ij[1]
        j    = ij[2]
        #i = which.min(abs(out$lon[q]-dat$lon2D))
        #j = which.min(abs(out$lat[q]-dat$lat2D))
        
        out$x[q] = dat$xc[i]
        out$y[q] = dat$yc[j]          
    } 

    return(out)
}

plot_to_file = FALSE 

domain    = "North"
grid_name = "NH-5KM"

if (TRUE) {
    # Plot initial figure 

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

    # Zoom on Northern Greenland/Ellesmere Island
    xlim = c(-3000,2000)
    ylim = c(-3000,0) 

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

    # Box containing GRL domain 
    x0 = -720.0
    nx = 43 
    dx = 40.0 
    x1 = x0 + (nx-1)*dx 
    y0 = -3450.0
    ny = 73 
    y1 = y0 + (ny-1)*dx 
    grlbox = data.frame(x0=x0,x1=x1,y0=y0,y1=y1)
    rect(grlbox$x0,grlbox$y0,grlbox$x1,grlbox$y1,border=2,col=NA)

    # Box containing EIS domain 
    x0 = 380.0
    nx = 89 
    dx = 40.0 
    x1 = x0 + (nx-1)*dx 
    y0 = -5000.0
    ny = 161 
    y1 = y0 + (ny-1)*dx 
    eisbox = data.frame(x0=x0,x1=x1,y0=y0,y1=y1)
    rect(eisbox$x0,eisbox$y0,eisbox$x1,eisbox$y1,border=2,col=NA)

    if (plot_to_file) {
        box()
        graphics.off() 
    }

}

if (TRUE) {
    # Read regions that exist and add them to figure 
    
    fldr = "polygons"


    poly1   = read.table(file.path(fldr,"polygon_laurentide.txt"),header=TRUE)
    coords1 = get_xy(poly1,topo)
    polygon(coords1$x,coords1$y,border='orange1')

    poly2   = read.table(file.path(fldr,"polygon_ellesmere.txt"),header=TRUE)
    coords2 = get_xy(poly2,topo)
    polygon(coords2$x,coords2$y,border='darkmagenta')

    poly3   = read.table(file.path(fldr,"polygon_grl.txt"),header=TRUE)
    coords3 = get_xy(poly3,topo)
    polygon(coords3$x,coords3$y,border='red')

    poly4   = read.table(file.path(fldr,"polygon_grl_and_ellesmere.txt"),header=TRUE)
    coords4 = get_xy(poly4,topo)
    polygon(coords4$x,coords4$y,border='white',lwd=1.5)

    # poly5   = read.table(file.path(fldr,"test.txt"),header=TRUE)
    # coords5 = get_xy(poly5,topo)
    # polygon(coords5$x,coords5$y,border='black',lwd=2)


    # pnow = poly1[1,]
    # i1 = which.min((pnow$x-coords4$x)^2+(pnow$y-coords4$y)^2)
    # pnow = poly1[2,]
    # i2 = which.min((pnow$x-coords4$x)^2+(pnow$y-coords4$y)^2)

    # pnow = poly1[1,]
    # i3 = which.min((pnow$x-coords4$x)^2+(pnow$y-coords4$y)^2)
    # pnow = poly1[2,]
    # i4 = which.min((pnow$x-coords4$x)^2+(pnow$y-coords4$y)^2)
    

}   


if (FALSE) {

    fldr = "new"

    ### Continental domains ###

    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_laurentide.txt"),scientific=TRUE,digits=4)

    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_eis.txt"),scientific=TRUE,digits=4)

    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_asia.txt"),scientific=TRUE,digits=4)

    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_grl.txt"),scientific=TRUE,digits=4)

    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    kk = rev(1:dim(poly1)[1])
    poly1  = poly1[kk,]
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_grl_extension.txt"),scientific=TRUE,digits=4)

    ### Sub-domains ###

    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_ellesmere.txt"),scientific=TRUE,digits=4)
    
    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_iceland.txt"),scientific=TRUE,digits=4)
    
    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_svalbard.txt"),scientific=TRUE,digits=4)
    
    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_britain.txt"),scientific=TRUE,digits=4)
    
    coords = locator(type="l")
    poly1  = get_latlon(coords,topo)
    my.write.table(poly1[,c("lon","lat")],file.path(fldr,"polygon_barents-kara.txt"),scientific=TRUE,digits=4)
    
}

if (FALSE) {
    points(poly_eis,col=1,pch=20,cex=1.2)

}

