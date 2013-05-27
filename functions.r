

grid.interp = function(dat,mask=FALSE,factor=2)
{
  # Get the dimensions
  nx = dim(dat)[1]
  ny = dim(dat)[2]
  x = c(1:nx)
  y = c(1:ny)
  
  obj    = list(x=x,y=y,z=dat)

  # Get new x and y coords
  x1 = seq(1,x[nx],length.out=nx*factor)
  y1 = seq(1,y[ny],length.out=ny*factor)
  newgrid = list(x=x1,y=y1)
  
  obj1   = interp.surface.grid(obj=obj,grid.list=newgrid)
  
  # Add a filter for whether it's a mask or not!!!
  # Temporary - could make land when ice+water touch!
  if (mask) { obj1$z <- round(obj1$z)  }

  return(obj1$z)
}

darker <- function(col,percent=10)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255 + (percent/100)
    cc[cc>1] <- 1; cc[cc<0] <- 0
    c[i] <- rgb(t(cc))
  }
  
  return(c)
}

colshade <- function(frac=0.1,alpha=1)
{ # +1 => lighter
  # -1 => darker 
  
  c = rep(0.5,length(frac)) + frac
  
  # Limit to range 0:1
  c[c<0] = 0
  c[c>1] = 1
  
  col = rgb(c,c,c,alpha)
  
  return(col)
}

alpha <- function(col,percent=50)
{
  c <- rep("#000000",length(col))
  if (length(percent)==1) p <- rep(percent,length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255
    c[i] <- rgb(t(cc),alpha=p[i]/100)
  }
  
  return(c)
}

my.png2 <- function(fldr=".",file="Rplot",date=FALSE,cairo="cairo",asp=1,win=NA,hin=NA,
                    res=300,pointsize=12,bg="white",ext=".png",units="mm",
                    cex=1,cex.lab=1,cex.axis=1)
{
  
  # Make filename
  dte <- ""
  if (date == TRUE) dte <- paste(today,"_",sep="")
  file <- file.path(fldr,paste(dte,file,ext,sep=""))
  
  host <- system("hostname",intern=TRUE)
  os   <- system("uname",intern=TRUE)
  
  # If running on a mac, make sure engine is quartz!
  #if ( os == "Darwin" ) cairo = "quartz"
  
  # Determine width/heights in inches
  if ( is.na(win) & is.na(hin) ) {  # Use default height, determine win via asp
    win <- 189  # Default width for pointsize 12
    hin <- win/asp
  } else if ( is.na(hin) ) {  # only win specified, determine hin
    hin <- win/asp
  } else if ( is.na(win) ) {  # only hin specified, determine win
    win <- asp*hin
  } else {                    # hin and win specified, determine asp
    asp <- win/hin
  }
  
  # Convert quantities if input was not inches
  cat(ext,":",file,"\n")
  cat("width =",win,", height =",hin,"(",units,".) \n")
  conv <- 1.0
  if ( units == "mm" ) conv <- 0.0393700787
  if ( units == "cm" ) conv <- 0.393700787
  hin <- hin*conv; win <- win*conv 
  cat("width =",win,", height =",hin," (in.)\n")
  
  if ( ext == ".png" ) {
    cat("cairo = ",cairo,"\n")
    png(file,width=win,height=hin,units="in",pointsize=pointsize,res=res,type=cairo)
  
  } else if ( ext == ".jpg" ) {
  
    jpeg(file,width=win,height=hin,units="in",pointsize=pointsize,res=res,type=cairo)
  
  } else if ( ext == ".tiff" ) {
  
    tiff(file,width=win,height=hin,units="in",pointsize=pointsize,res=res,type=cairo)
  
  } else if ( ext == ".pdf" ) {
    
    if (cairo %in% c("cairo","cairo1")) {
      cat("**cairo_pdf","\n")
      cairo_pdf(file,width=win,height=hin,pointsize=pointsize)
    } else {
      pdf(file,width=win,height=hin,pointsize=pointsize)
    }
    
  } else if ( ext == ".svg" ) {
  
    svg(file,width=win,height=hin,pointsize=pointsize)
  
  } else if ( ext == ".fig" ) {
  
    xfig(file,width=win,height=hin,pointsize=pointsize)
    
  } else {  # if ( ext == ".ps" )
    
#    postscript(file,width=win,height=hin,pointsize=pointsize,
#               paper="special")
    cairo_ps(file,width=win,height=hin,pointsize=pointsize)
  }
  
#  if ( ext == ".ps" ) {
#     cex <- cex*0.5; cex.axis <- cex.axis*0.5; cex.lab <- cex.lab*0.5
#  }
  
  par(bg=bg,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,tcl=0.2,mgp=c(2.5,0.3,0),las=1)
  return(win)
}

my.par  <- function(mar=c(3.2,3.3,1,1),xaxs="i",yaxs="i",tcl=0.4,mgp=c(2.5,0.3,0),las=1,...)
{
  par(...,mar=mar,tcl=tcl,mgp=mgp,las=las,xaxs=xaxs,yaxs=yaxs)
}

my.axis <- function(side=1,tcl=0.4,mgp=c(2.5,0.25,0),minticks=2,grid=TRUE,...)
{
  
  if ( side %in% c(2,4)) mgp[2]=mgp[2]*1.3
  par(mgp=mgp,las=1)
  #if (side %in% c(1,3)) par(mgp=mgph)
  #if (side %in% c(2,4)) par(mgp=mgpv)
  
  # Calculate the main axis
  at.main = axis(1,labels=F,tick=F)

  # If desired, add ticks between standard values
  if (minticks > 1) {
    dx1 = diff(at.main)[1]/minticks
    at.min = axis(side=side,at=seq(range(at.main)[1]-10*dx1,range(at.main)[2]+10*dx1,by=dx1),labels=FALSE,tick=FALSE)
  } else {
    at.min = at.main
  }
  
  if (grid & side %in% c(1,3)) abline(v=at.min,lwd=1,lty=3,col=8)
  if (grid & side %in% c(2,4)) abline(h=at.min,lwd=1,lty=3,col=8)
  
  # Now actually draw ticks and labels...
  if (!grid & minticks > 1)
    axis(side=side,at=at.min, tcl=tcl/2,labels=FALSE)
  
  # Now plot main axes with labels, etc...
  axis(side=side,at=at.main,tcl=tcl,...)

  return(at.min)
}

my.plot <- function(...,axes=c(1,2,3,4),box=TRUE,grid=TRUE,
                    xlab="",ylab="",xline=1.7,yline=1.9,
                    minticks=2)
{
  #my.par(las=1)
  plot(...,axes=F,ann=F)
  #if (grid) grid()

  if ( 1 %in% axes) my.axis(side=1,minticks=minticks,grid=grid)
  if ( 2 %in% axes) my.axis(side=2,minticks=minticks,grid=grid)
  if ( 3 %in% axes) my.axis(side=3,labels=F,minticks=minticks)
  if ( 4 %in% axes) my.axis(side=4,labels=F,minticks=minticks)

  if (box) box()
  
  if ( !xlab=="" ) title(xlab=xlab,line=xline)
  if ( !ylab=="" ) title(ylab=ylab,line=yline)

}

fliplat <- function(dat,nc)
{ # Reverse the latitude dimension of all relevant variables of a data set
  
  vnms = names(nc$var)

  nms = NULL
  for (q in 1:length(vnms)) { 
    vnm = vnms[q]
    nd = nc$var[[vnm]]$ndim 
    latnow = FALSE
    for ( d in 1:nd ) { 
      if ("lat" %in% nc$var[[vnm]]$dim[[d]]$name) 
            nms = c(nms,vnm)
    }
  }

  cat("fliplat:",nms,"\n")

  # Latitude should be a vector, reverse it as needed
  if ( dat$lat[1] > dat$lat[2] ) {

    jj = rev(c(1:length(dat$lat)))
    dat$lat = dat$lat[jj]
    
    # Now reverse all variables 
    for ( q in 1:length(nms) ) {
      nm = nms[q]
      dims = dim(dat[[nm]])
      
      if ( nm == "lat_bnds" ) {

        dat[[nm]] = dat[[nm]][,jj]

      } else if ( length(dims)== 1 ) {

        dat[[nm]] = dat[[nm]][jj]
      
      } else if ( length(dims)== 2 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][,jj]

      } else if ( length(dims)==3 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][,jj,]

      } else if ( length(dims)==4 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][,jj,,]

      }

    } 
  
  }
  
  return(dat)

}

shiftlon <- function(dat,nc)
{ # Shift longitude from 0:360 => -180:180 for
  # all related variables 
  
  vnms = names(nc$var)

  nms = NULL
  for (q in 1:length(vnms)) { 
    vnm = vnms[q]
    nd = nc$var[[vnm]]$ndim 
    latnow = FALSE
    for ( d in 1:nd ) { 
      if ("lon" %in% nc$var[[vnm]]$dim[[d]]$name) 
            nms = c(nms,vnm)
    }
  }
  
  cat("shiftlon:",nms,"\n")

  if (max(dat$lon) > 180 ) {    # Shift from 0:360 to -180:180

    # Longitude should be a vector
    ii = c(1:length(dat$lon))

    # Longitude
    lonX = dat$lon
    i = which(dat$lon > 180)
    lonX[i] = lonX[i] - 360
    i1 = which(lonX<0)
    i0 = which(lonX>0)
    ii = c(i1,i0)
    
    dat$lon = lonX[ii]
  
    # Now reverse all related variables 
    for ( q in 1:length(nms) ) {
      nm = nms[q]
      dims = dim(dat[[nm]])

      if ( nm == "lon_bnds" ) {

        dat[[nm]] = dat[[nm]][,ii]

      } else if ( length(dims)== 1 ) {

        dat[[nm]] = dat[[nm]][ii]
      
      } else if ( length(dims)== 2 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][ii,]

      } else if ( length(dims)==3 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][ii,,]

      } else if ( length(dims)==4 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][ii,,,]

      }

    } 
    
  }

  return(dat)

}

gen_time = function(info)
{ # Given info for one file, return years and months corresponding to data points
  year0  = info$year0
  yearf  = info$yearf
  month0 = info$month0 
  monthf = info$monthf
  
  # Get start, middle (excluding 1st and last, since they may have less than 12 months)
  # and end years for data set
  nyr0 = 12-month0+1
  yearsstart = rep(year0,nyr0)

  yearsmid = NULL
  if ( yearf-year0 > 1 ) {
    yearsmid = rep(c((year0+1):(yearf-1)),12)
    dim(yearsmid) = c(length(yearsmid)/12,12)
    yearsmid  = as.vector(t(yearsmid))
  }
  
  yearsend = NULL
  if ( yearf-year0 > 0 ) {
    nyrf = monthf
    yearsend = rep(yearf,nyrf)
  }
  # Make a vector for all years now
  nyr0 = 12-month0+1
  nyrf = monthf
  years0 = c(yearsstart,yearsmid,yearsend)


  # Line up time data with months starting with month0 and ending with monthf
  monthsstart = c(month0:12)
  monthsmid = NULL
  if (yearf-year0 > 1 ) monthsmid = rep(c(1:12),(yearf-year0-1))
  monthsend = NULL
  if (yearf-year0 > 0 ) monthsend = c(1:monthf)
  months0 = c( monthsstart,monthsmid,monthsend )

  return(data.frame(years0=years0,months0=months0))
}

arrange_3Dto4D <- function(dat,d3,d4)
{
  
  return(dat4D)
}

get.nc <- function(file,convert=1,vnms=NA,month=-1,missing=1e10,fliplat=FALSE,shiftlon=FALSE)
{
  nc <- open.ncdf(file)
  if (is.na(vnms[1])) vnms <- names(nc$var)
  dims = nc$dim 

  # Definet the initial list for output
  # including all the nc information
  out <- list(nc=nc)
  
  # Get dimensional variables, if they exist
  for ( q in 1:length(dims) ) {
    nm = dims[[q]]$name 
    out[[nm]]   <- get.var.ncdf(nc,nm)
  }

  # Load the variables
  for ( j in 1:length(vnms) ) {
    vnm <- vnms[j]
    
    v3      <- nc$var[[vnm]]
    varsize <- v3$varsize
    ndims   <- v3$ndims

    # Initialize start and count to read one timestep of the variable.
    start <- rep(1,ndims)   # begin with start=(1,1,1,...,1)
    start[ndims] <- 1       # change to start=(1,1,1,...,i) to read timestep i
    count <- varsize        # begin w/count=(nx,ny,nz,...,nt), reads entire var
    count[ndims] <- -1      # change to count=(nx,ny,nz,...,1) to read 1 tstep
    
    if ( month != -1 ) {

      if ( v3$dim[[1]]$name == "month" ) {
        start[1] <- month
        count[1] <- 1
      } else if ( ndims > 1 ) {
        if( v3$dim[[2]]$name == "month"  ) {
          start[2] <- month
          count[2] <- 1
        }
      }

    }
    
    if ( ! vnm %in% c("empty1","empty2","empty3","empty4","empty5") ) {
      #cat("vnm:",vnm,"\n")
      out[[vnm]] <- get.var.ncdf( nc, v3, start=start, count=count )
      #out[[vnm]] <- get.var.ncdf(nc,vnm,start=c(0,13),count=c(-1,13))
      
      if (is.numeric(out[[vnm]][1])) {
        out[[vnm]][abs(out[[vnm]]) >= missing] <- NA
        out[[vnm]][abs(out[[vnm]]) < 1e-12] = 0.0
      }
    }

  }

  # Make a conversion to ka, if needed
  # newer output files have output in years, so they need to be consistent!!
  if (!is.null(out$time)) {
    
    # Check time vector units, if not in kiloyears, convert as appropriate
    if ( pmatch("kiloyears",nc$dim$time$units,nomatch=0) == 0 ) {
      out$time <- out$time*1e-3
    }
    
    # Round to the nearest year, and then apply conversion desired
    out$time <- round(out$time*1e3,0)*1e-3
    out$time <- out$time *convert
  }
  
  # Make sure latest Vtot variable was obtained
  if (!is.null(out$Vtot)) {
  
    n <- length(out$Vtot)
    if ( is.na(out$Vtot[n]) ) out$Vtot[n] <- out$Vtot[n-1]
  
  }
  
  # Check if we are flipping a latitude dimension
  # (up to a three dimensional array assuming lat in second position!)
  if ( fliplat & "lat" %in% names(dims) )
    out = fliplat(dat=out,nc=nc)

  if ( shiftlon & "lon" %in% names(dims) )
    out = shiftlon(dat=out,nc=nc)

  cat("get.nc :",file,"\n")
  
  # Close files
  close.ncdf(nc)
  
  return(out)
}

get.nc0 <- function(file,convert=1,vnms=NA,month=-1,missing=1e10)
{
  nc <- open.ncdf(file)

  if (is.na(vnms[1])) vnms <- names(nc$var)

  # Definet the initial list for output
  out <- list()
  
  # Save meta info
  out$file <- file
  
  # Get dimensional variables, if they exist
  try( out$time <- get.var.ncdf(nc,"time"), silent=TRUE ) 
  try( out$lat  <- get.var.ncdf(nc,"lat"), silent=TRUE ) 
  try( out$lon  <- get.var.ncdf(nc,"lon"), silent=TRUE ) 
  
  # Load the variables
  for ( j in 1:length(vnms) ) {
    vnm <- vnms[j]
    
    v3      <- nc$var[[vnm]]
    varsize <- v3$varsize
    ndims   <- v3$ndims

    # Initialize start and count to read one timestep of the variable.
    start <- rep(1,ndims)   # begin with start=(1,1,1,...,1)
    start[ndims] <- 1       # change to start=(1,1,1,...,i) to read timestep i
    count <- varsize        # begin w/count=(nx,ny,nz,...,nt), reads entire var
    count[ndims] <- -1      # change to count=(nx,ny,nz,...,1) to read 1 tstep
    
    if ( month != -1 & v3$size[1] == 13 ) {
      start[1] <- month
      count[1] <- 1
    }
    
    #cat("vnm:",vnm,"\n")
    out[[vnm]] <- get.var.ncdf( nc, v3, start=start, count=count )
    #out[[vnm]] <- get.var.ncdf(nc,vnm,start=c(0,13),count=c(-1,13))

    out[[vnm]][abs(out[[vnm]]) >= missing] <- NA
  }

  # Make a conversion to ka, if needed
  # newer output files have output in years, so they need to be consistent!!
  if (!is.null(out$time)) {
    
    # Check time vector units, if not in kiloyears, convert as appropriate
    if ( pmatch("kiloyears",nc$dim$time$units,nomatch=0) == 0 ) {
      out$time <- out$time*1e-3
    }
    
    # Round to the nearest year, and then apply conversion desired
    out$time <- round(out$time*1e3,0)*1e-3
    out$time <- out$time *convert
  }
  
  # Make sure latest Vtot variable was obtained
  if (!is.null(out$Vtot)) {
  
    n <- length(out$Vtot)
    if ( is.na(out$Vtot[n]) ) out$Vtot[n] <- out$Vtot[n-1]
  
  }
  
  #out$x    <- try( get.var.ncdf(nc,"x") )
  #out$y    <- try( get.var.ncdf(nc,"y") )

  #   ii_ice <- which(out$mask == 0)
  #   ii_land <- which(out$mask == 1)
  #   ii_water <- which(out$mask == 2)

  # If it's 1-d, convert it into a data frame
  #if ( length(dim(out[[1]])) ) out <- as.data.frame(out)
  
  cat("get.nc :",file,"\n")
  
  # Close files
  close.ncdf(nc)
  
  return(out)
}

sim.filter <- function(data,x1,x0=NA,fill=TRUE,rule=1)
{
  
  # Extract time vector, delete it from data 
  if (is.na(x0[1])) {x0 <- data[[1]]; data[[1]] <- NULL}

  # Set up a data frame with nrows of x1 vector, but ncol of actual data names
  data2 <- data.frame(x1=x1*NA)
  for ( q in 2:ncol(data) ) data2[[q]] <- x1*NA
  names(data2) <- names(data)
  
  # First fill in values that match (rounding to ensure a match!!)
  ii <- round(x0*1e3) %in% round(x1*1e3)
  jj <- round(x1*1e3) %in% round(x0*1e3)
  data2[jj,] <- data[ii,]
  
  # Loop over each variable and interpolate if desired
  if ( fill ) {
    
    # Check which times dont have data (for first column only - should be identical for all)
    # Fill in first missing value manually if needed, then check remaining values
    if ( is.na(data2[1,1]) ) data2[1,] <- data2[2,]    
    ii <- which(is.na(data2[,1]))
    
    if ( length(ii) > 0 ) {
      for ( q in 1:ncol(data2) ) data2[ii,q] <- approx(x0,data[[q]],xout=x1[ii],rule=rule)$y
    }
  }
  
  # Add time vector to beginning of data 
  data2 <- cbind(data.frame(time=x1),data2)
  
  return(data2)
}

get_sector_contribs <- function(dat,sectors,nms=c("pp","snow","melt","runoff","smb"),dx=20)
{ ## Figure out sector contributions to SMB
  
  # Get conversion factor (mm=>Gt)
  conv <- (dx*1e3)^2*1e-12
  
  ns <- max(sectors) 
  ii0 <- which(dat$mask==0)
  
  sect <- as.data.frame(array(NA,dim=c(ns,length(nms)+2)))
  names(sect) <- c("sect","area",nms)
  
  for ( i in 1:ns ) {
    q <- which(sectors==i & dat$mask==0)
    
    sect$sect[i]   <- i
    sect$area[i]   <- length(q)/length(ii0)
    sect$pp[i]     <- sum(dat$pp[q])*conv
    sect$snow[i]   <- sum(dat$snow[q])*conv
    sect$melt[i]   <- sum(dat$melt[q])*conv
    sect$runoff[i] <- sum(dat$runoff[q])*conv
    sect$smb[i]    <- sum(dat$smb[q])*conv
  }
  
  return(sect)
}

get_sector <- function(dat,latmid=72)
{ ## Return an array specifying the sector of each index
  
  zs  <- dat$zs
  lat <- dat$lat; lon <- dat$lon
  x   <- dat$xx;  y <- dat$yy
  
  ridge <- get_ridge(zs,x,y,lat,lon,ymin=-2900,ymax=-1700)
  
  ny <- dim(zs)[2]
  
  sector <- zs*NA
  
  for (q in 1:ny) {
    sector[,q] <- ifelse(lon[,q] > ridge$lon[q],1,2)
  }
  
  sector[lat > latmid & sector==2] <- 3
  sector[lat > latmid & sector==1] <- 4
  
  return(sector)
}

get_ridge <- function(zs,x,y,lat,lon,ymin=-2900,ymax=-1700,spar=0.6)
{ ## Figure out where the central ridge line is based on maximum elevations
  ny <- dim(zs)[2]

  # Make some output vectors
  mid <- data.frame(j=numeric(ny),i=numeric(ny))
  mid$x <- mid$y <- mid$x2 <- mid$i2 <- NA
  mid$lat <- mid$lon <- NA
  
  for ( q in 1:ny ) {
    mid$y[q] <- y[1,q]; mid$j[q] <- q
    
    i <- which.max(zs[,q])
    mid$x[q] <- x[i,q]; mid$i[q] <- i
    
    i <- round(mean(order(zs[,q],decreasing=TRUE)[1:25]))
    mid$x2[q] <- x[i,q]; mid$i2[q] <- i
    
    mid$lat[q] <- lat[i,q]
    mid$lon[q] <- lon[i,q]
  }
  
  #r2 <- smooth.spline(mid$y,mid$x,spar=0.6)
  #mid$xsm <- r2$y; mid$ysm <- r2$x
  
  # Limit the ridge to certain y-values
  j0 <- which( round(mid$y,1) == ymin); j1 <- which( round(mid$y,1) == ymax)
  
  mid$x3 <- mid$x2; mid$y3 <- mid$y
  jj <- which(mid$y < ymin)
  mid$x3[jj] <- mid$x[j0];  mid$y3[jj] <- mid$y[j0]
  #mid$lon[jj] <- mid$lon[j0]; mid$lat[jj] <- mid$lat[j0]
  
  jj <- which(mid$y > ymax)
  mid$x3[jj] <- mid$x2[j1]; mid$y3[jj] <- mid$y[j1]
  #mid$lon[jj] <- mid$lon[j1]; mid$lat[jj] <- mid$lat[j1]
  
  return(mid)
}

get_profile.dist <- function(dat,ii=which(dat$mask==0))
{ ## Figure out profiles of Greenland topography

  ii.water <- which(dat$mask==2)
  dist <- dat$mask*NA
  
  for ( i in ii ) {
    dists <- calcdist(data.frame(x=as.vector(dat$x[ii.water]),y=as.vector(dat$y[ii.water])),
                      data.frame(x=dat$x[i],y=dat$y[i]))
    dist[i] <- min(dists,na.rm=TRUE)
  }
  
  return(dist)
}

get.slice <- function(dat,index=NA,time=NA,nms=c("mask","zs"),ebs=1e-3)
{ ## Get a specific time slice from a 2d dataset
  
  # Determine index of correct time
  if ( !is.na(time) ) {
    t <- which( abs(dat$time-time) < ebs)[1]
  } else {
    t <- index
    if (is.na(t)) t <- which.max(dat$time)
  }
  
  cat("get.slice: time=",dat$time[t],"\n")
  
  # First load generic grid variables (if they exist)
  out <- list()
  if ( "lat" %in% names(dat) ) out <- list(lat=dat$lat,lon=dat$lon,xx=dat$xx,yy=dat$yy)

  # Then load desired variables
  for ( q in 1:length(nms) ) {
    nm <- nms[q]
    if ( is.na(t) | nm %in% c("mask_hydro") ) {
      out[[nm]] <- dat[[nm]]
    } else {
      out[[nm]] <- dat[[nm]][,,t]
    }
  }
  
#   if ( "zs" %in% nms ) {
#     
#     i <- which.max(out$zs)
#     out$x.max <- grid$xx[i]; out$y.max <- grid$yy[i]
#     out$zs.max <- out$zs[i]
#     
#   }
  
  return(out)
}


get_declineTiming <- function(time,Vol,V0=Vol[1],percent=c(20,50,80,85,90,100))
{ # Percent is the percent melted
  
  Vpp <- 100 * (Vol / V0)

  outV <- outt <- numeric(length(percent))*NA
  
  for ( q in 1:length(percent) ) {
    
    ii <- which( Vpp <= (100-percent[q]) )
    
#     if ( length(ii) == 0 ) ii <- c(1:length(Vpp))
    
    if ( length(ii) > 0 ) {
      i <- ii[which.max(Vpp[ii])]
      outV[q] <- Vpp[i]
      outt[q] <- time[i]
    }
    
  }
  
  out <- as.data.frame(as.list(c(outV,outt)))
  names(out) <- c(paste("V",percent,sep=""),paste("t",percent,sep=""))
  
  return(out)

}

get_declineTiming2 <- function(time,Vol,V0=NULL,percent=20)
{ # Percent is the percent melted

  # Make sure Vol is a 2d array, where dim 1 is the simulation number
  # (for a vector, we only have one simulation)
  if (is.null(dim(Vol))) dim(Vol) = c(1,length(Vol))
  
  # How many simulations
  nsim = dim(Vol)[1]

  # Make sure the reference volume is defined  
  if (is.null(V0)) V0 = Vol[,1]
  
  # Get the volume percentage relative to reference volume
  Vpp <- 100 * (Vol / V0)
  
  out = list()
  out$timing  = array(NA,dim=c(nsim,length(percent)))
  out$volume  = array(NA,dim=c(nsim,length(percent)))
  out$percent = t(array(percent,dim=rev(dim(out$timing))))

  for (q in 1:nsim) {
    for (j in 1:length(percent)) {
      
      ii = which( Vpp[q,] <= (100-percent[j]) )

      if ( length(ii) > 0 ) {
        i = ii[which.max(Vpp[q,ii])]
        out$timing[q,j] = time[i]
        out$volume[q,j] = Vpp[q,i]
      }
    }
  }

  return(out)

}

get.param.values <- function(fldrs,file=c("rembo.params","sico.params"),nms=NA,comment="#",skip=4)
{
  
  for ( q in 1:length(fldrs) ) {
    
    fldr <- fldrs[q]
    parameters1 <- matrix(scan(file.path(fldr,file[1]),what="character",comment.char=comment,skip=skip),byrow=TRUE,ncol=3)  
    
    parameters <- parameters1
    
    if ( length(file) > 1 ) {
      if ( file.exists(file.path(fldr,file[2])) ) {
        parameters2 <- matrix(scan(file.path(fldr,file[2]),what="character",comment.char=comment,skip=skip),byrow=TRUE,ncol=3)
        parameters  <- rbind(parameters1,parameters2)
      }
    }

    ## new ##
#     params <- as.list(parameters[,3])
#     names(params) <- parameters[,1]
#     for (i in 1:length(params)) {
#       p <- as.character(as.numeric(params[[i]]))
#       cat(names(params)[i]," : ","p =",p,"  : ","params[[i]] =",params[[i]],"\n")
#       if (!is.na(p) & p == params[[i]]) params[[i]] <- as.numeric(params[[i]])
#     }
    #######
    
#     params <- as.list((as.numeric(parameters[,3])))
    params <- as.list(parameters[,3])
    names(params) <- parameters[,1]
    
    nmstrings <- c("rcpname","ppsname","temper_file")
    
    ii <- which(! parameters[,1] %in% nmstrings )
    for (i in ii) {
      params[[i]] <- as.numeric(params[[i]])
    }
    
    if ( length( grep("PPS",params$temper_file) ) > 0 ) {
      params$ppsname <- substring(params$temper_file,first=14,last=18)
    }
    
    
    if (q == 1) { 
      all <- as.data.frame(params,stringsAsFactors=FALSE)
    } else {
      all <- rbind(all,as.data.frame(params,stringsAsFactors=FALSE))
    }
  
  }
  
  # Get indices of matching names
  pii <- c(1:dim(all)[2])
  if (!is.na(nms[1])) pii <- match(nms,names(all))
  
  bad <- which(is.na(pii))
  if (length(bad) > 0) {
    cat("Invalid parameter name(s):",nms[bad],"\n")
    pii <- pii[which(!is.na(pii))]
  }
  
  # Filter down to just the names I want
  all <- all[,pii]
  
  # Hack to account for changes to parameter names in rembo and sico
  if ( "margin_value" %in% names(all) ) all$margin_value <- NULL
  if ( "anf_frac" %in% names(all) )     all$anf_frac <- NULL
  if ( "rembo_start_file"  %in% names(all) ) { all$restart_file <- all$rembo_start_file; all$rembo_start_file <- NULL }
  
  # For rembo parameters only
  if ( "itm_c" %in% names(all) ) {
    if ( (!"itm_b" %in% names(all)) ) all$itm_b <- 0.0
    if ( (!"year_offset" %in% names(all)) ) all$year_offset <- 0.0
    if ( (!"T_diff" %in% names(all)) ) all$T_diff <- 0.0
  }
  
  # Modify some parameter names (backwards compatibility), if both rembo and sico used
  if ( "itm_c" %in% names(all) & "C_SLIDE_0" %in% names(all) ) {
    if ( "slide0" %in% names(all) ) { all$C_SLIDE_0 <- all$slide0; all$slide0 <- NULL }
    if ( "slide1" %in% names(all) ) { all$C_SLIDE_SEDI <- all$slide1; all$slide1 <- NULL }
    if ( "sea_level" %in% names(all) ) { all$SEA_LEVEL <- all$sea_level; all$sea_level <- NULL }
    if ( "ice_stream" %in% names(all) ) { all$ICE_STREAM <- all$ice_stream; all$ice_stream <- NULL }
    if ( "q_geo"  %in% names(all) ) { all$Q_GEO_0 <- all$q_geo; all$q_geo <- NULL }
    if ( "stream_cut.1" %in% names(all) ) { all$stream.cut.1 <- NULL }
  }
  
  #cat("fldr: ",fldrs,"\n",names(all),"\n")
  
  return(all)
  
}

get.var <- function(file,ny=141,nx=76,what="double",comment="#")
{
  if ( what == "integer" ) {
    dat <- matrix(scan(file,comment=comment,what="character"),byrow=TRUE,nrow=ny,ncol=1)
    out <- matrix(NA,nrow=ny,ncol=nx)
    for ( q in 1:ny ) {
      out[q,] <- as.numeric(strsplit(dat[q,],"")[[1]])
    }
    out <- t(out)[,ny:1]
    
  } else {
    out <- t(matrix(scan(file,comment=comment),byrow=TRUE,nrow=ny,ncol=nx))[,ny:1]
  }
  
  return(out)
}


topo.uplifted <- function(zs,zb,rho=910.0,rho_a=3300.0)
{
  return( zb + (rho/rho_a)*(zs-zb) )
}

topo.stats <- function(zs,zb,mask=NA,dx=20,print=FALSE)
{
  
  # Define a mask
  if (is.na(mask[1])) {
    mask <- zs*0+2; mask[zs>=0] <- 1; mask[zs-zb>0] <- 0
  }
  
  conv.area <- 1e-6 * 1e-6  # m2 => 1e6 km2
  conv.vol  <- 1e-9 * 1e-6  # m3 => 1e6 km3
  
  # Get area in m2
  dx2 <- dx*dx*(1e3*1e3)
  
  H <- (zs - zb)
  ice.vol <- sum(H*dx2) * conv.vol
  
  tmp <- mask*0; tmp[mask==0] <- dx2
  ice.area <- sum(tmp) * conv.area
  
  tmp <- mask*0; tmp[mask<=1] <- dx2
  land.area <- sum(tmp) * conv.area
  
  if (print) {
    cat("\n","    ","  Land","  Ice","\n")
    cat("Area",round(land.area,3),round(ice.area,3),"\n")
    cat("Vol "," --- ",round(ice.vol,3),"\n")
    cat("\n","Elev. range:",round(min(zs),1),round(max(zs),1),"\n")
    cat("\n")
  }
  
  return(data.frame(ice.vol=ice.vol,ice.area=ice.area,land.area=land.area))
}

surface.normal = function(zs)
{
  nx <- dim(zs)[1]; ny <- dim(zs)[2]
  dxx = dyy = dzz = matrix(0,nrow=nx,ncol=ny)

  # Set all vector-z values to 1, since it is pointing away from surface
  dzz[] = 1

  # Get mean elevation of neighbors, then 
  # calculate the xy-gradient
  for ( i in 2:(nx-1) ) {
    for ( j in 2:(ny-1) ) {
      
        zsmean = mean(zs[(i-1):(i+1),(j-1):(j+1)],na.rm=TRUE)

        dx <- c( zs[i,j] - zs[i-1,j  ], zs[i+1,j  ] - zs[i,j] )
        dy <- c( zs[i,j] - zs[i  ,j-1], zs[i  ,j+1] - zs[i,j] )
        
        dx[is.na(dx)] = 0
        dy[is.na(dy)] = 0
        
        dxx[i,j] = mean(dx)
        dyy[i,j] = mean(dy)

    }
  }

  # Normalize
  mag = sqrt(dxx^2+dyy^2+dzz^2)
  dxx = dxx/mag 
  dyy = dyy/mag 
  dzz = dzz/mag
  
  sn = array(0,dim=c(nx,ny,3))
  sn[,,1] = dxx
  sn[,,2] = dyy
  sn[,,3] = dzz

  #return(list(dxx=dxx,dyy=dyy,dzz=dzz))
  return(sn)
}

surface.normal2 = function(zs)
{
  nx <- dim(zs)[1]; ny <- dim(zs)[2]
  
  dxx = dyy = dzz = matrix(0,nrow=nx,ncol=ny)
  
  dx1 = zs[3:nx,1:ny] - zs[2:(nx-1),1:ny]
  dx2 = zs[2:(nx-1),1:ny] - zs[1:(nx-2),1:ny]
  dxx[2:(nx-1),] = (dx1+dx2) / 2

  dy1 = zs[1:nx,3:ny] - zs[1:nx,2:(ny-1)]
  dy2 = zs[1:nx,2:(ny-1)] - zs[1:nx,1:(ny-2)]
  dyy[,2:(ny-1)] = (dy1+dy2) / 2
  
  # Normalize
  mag = sqrt(dxx^2+dyy^2+dzz^2)
  dxx = dxx/mag 
  dyy = dyy/mag 
  dzz = dzz/mag
  
  sn = array(0,dim=c(nx,ny,3))
  sn[,,1] = dxx
  sn[,,2] = dyy
  sn[,,3] = dzz
  
  return(sn)

}


dotprod = function(x,y) return( sum(x*y) )

vec.angle = function(x,y)
{
  xm = sqrt(sum(x^2))
  ym = sqrt(sum(y^2))
  xy = dotprod(x,y)

  costheta = xy / (xm*ym)

  theta = acos(costheta)

  return(theta*180/pi)
}

hgrad <- function(xy,dx=20e3)
{ # xy is 2D field, calculate horizontal gradient vector, dx is grid resolution
  
  nx = dim(xy)[1]
  ny = dim(xy)[2]
  
  dudx = xy*NA
  dudy = xy*NA
  
  inv_2dx = 1/(2*dx)

  for (i in 2:(nx-1)) {
    for (j in 2:(ny-1)) {

      dudx[i,j] = (xy[i+1,j] - xy[i-1,j]) * inv_2dx
      dudy[i,j] = (xy[i,j+1] - xy[i,j-1]) * inv_2dx
        
    }
  }
  
  # Get magnitude too
  mag = sqrt(dudx^2 + dudy^2)

  return(list(x=dudx,y=dudy,xy=mag))
}

col.grad <- function(zs,ltheta=-135,lphi=20)
{ # Generate a scalar to increase or decrease shading
  # -1: darkest color
  #  0: no change
  #  1: lightest color 
  
  # Calculate vector of sunshine (from sun to point)
  # negative, bc points in towards origin...
  torads = pi/180
  dx = sin(lphi*torads)*cos(ltheta*torads)
  dy = sin(lphi*torads)*sin(ltheta*torads)
  dz = -cos(lphi*torads)
  vsun = c(dx,dy,dz)

  # First get surface normals
  vsn = surface.normal2(zs)
  
  # Make into a handy vector for calculating angles
  #vsn = cbind(as.vector(sn$dxx),as.vector(sn$dyy),as.vector(sn$dzz))
  
  # Calculate angle of surface from sunshine for each point
  angles = apply(vsn,FUN=vec.angle,MARGIN=c(1,2),y=vsun)

  # angles = vsn[1,]*0
  # for ( k in 1:dim(vsn)[1] ) {
  #   angles[k] = vec.angle(vsun,vsn[k,])
  # }
  #dim(angles) = dim(zs)
  
  # The larger the angle, the brighter the color should be,
  # 180: directly pointing at surface
  #   0: directly pointing away from surface
  #colshift = (angles/180 - 0.5)
  
  # Get the partially transparent colormask
  # color is black, full transparency or none for brightness
  colshift = 1-(angles/180)
  #colshift = rgb(0,0,0,alpha=colshift)

  return(colshift)
}

add_shade <- function(zs,plot=TRUE)
{
  cat("Calculating contour shading...")
  colshift = col.grad(zs)
  colshade=rgb(1,1,1,seq(0,40,length.out=100),maxColorValue=100)
  cat("done.\n")
  
  if (plot) image(colshift,col=colshade,add=T)

  return(list(level=colshift,col=colshade)) 
}

test.shading <- function(zs,var=NULL)
{

  shade = add_shade(zs,plot=FALSE)

  elev = get.contours("zs")
  image.plot(a$zs,breaks=elev$at,col=elev$cols)
  image(shade$level,col=shade$col,add=T)
  
  breaks = pretty(range(var,na.rm=T),50)
  cols   = colorRampPalette(jet.colors)(length(breaks)-1)
  cols = alpha(cols,65)
  if (!is.null(var)) image(var,add=T,breaks=breaks,col=cols)
  
  contour(a$zs,add=T,levels=seq(0,3500,by=250),drawlabels=FALSE,col="grey50")
  contour(a$zs,add=T,levels=c(2500),drawlabels=FALSE,col="grey50",lwd=3)

}


topo.grad <- function(zs,dx=20,max=TRUE)
{
  nx <- dim(zs)[1]; ny <- dim(zs)[2]

  dzs <- matrix(0,nrow=nx,ncol=ny)

  for ( i in 2:(nx-1) ) {
    for ( j in 2:(ny-1) ) {

      d1 <- zs[i-1,j  ] - zs[i,j]
      d2 <- zs[i+1,j  ] - zs[i,j]
      d3 <- zs[i  ,j-1] - zs[i,j]
      d4 <- zs[i  ,j+1] - zs[i,j]

      dzs[i,j] <- max(abs(c(d1,d2,d3,d4))) / dx
      if (!max) dzs[i,j] <- mean(c(d1,d2,d3,d4) / dx)
    }
  }

  return(dzs)
}

get.asp <- function(zb,dx=20,mean=FALSE)
{
  
  dzb <- zb*0
  
  nx <- dim(zb)[1]; ny <- dim(zb)[2]
  
  for ( i in 2:(nx-1) ) {
    for ( j in 2:(ny-1) ) {
      
      dz <- c( zb[i,j] - zb[i-1,j],
               zb[i,j] - zb[i+1,j],
               zb[i,j] - zb[i,j-1],
               zb[i,j] - zb[i,j+1] )
      
      if (mean) {
        dz <- mean(abs(dz)) / dx
      } else {
        dz <- max(abs(dz)) / dx
      }
      
      dzb[i,j] <- dz
    }
  }
  
  return(dzb)
}

get.curvature <- function(zs,dx=20,n=1)
{ # mean curvature
  # (from: Park et al, 2001)
  # http://casoilresource.lawr.ucdavis.edu/drupal/node/937

  dists = array(1,dim=c(n*2+1,n*2+1))
  for (i in 1:nrow(dists)) {
    for (j in 1:ncol(dists)) {
      dists[i,j] = sqrt( (dx*((n+1)-i))^2 + (dx*((n+1)-j))^2 )
    }
  }
  qq = which(dists != 0)

  czs = zs*NA

  nx = nrow(czs)
  ny = ncol(czs)

  for ( i in (1+n):(nx-n)) {
    for ( j in (1+n):(ny-n)) {

      neighbs = zs[(i-n):(i+n),(j-n):(j+n)]
      
      czs[i,j] = sum((neighbs[qq]-zs[i,j])/dists[qq]) / length(qq)
    }
  }

  return(czs)
}

get.Lfac <- function(asp)
{ # From Surenda's poster
  
  Lfac <- -0.7*asp^2 - 0.18*asp + 1.0
  return(Lfac)
}

stream.gen <- function(zs,fmax=1,zmax=500)
{ ## Generate a sliding mask, where the sliding factor
  ## increases with lower elevations
  
  # First make a mask of zeros (assuming no sliding anywhere)
  mask <- zs*0
  
  # Obtain points that are of interest, ie, less than zmax
  ii <- which(zs <= zmax)
  
  # Determine how strong the sliding factor is at each point
  # (the lower the elevation, the stronger the sliding factor
  mask[ii] <- fmax*(1-zs[ii]/zmax)
  
  # If the elevation is below zero, or the mask
  # erroneously gives a value greater than fmax, reset to fmax
  mask[zs <= 0 | mask > fmax] <- fmax
  
  return(mask)
}

filter.grid <- function(data,dx,dx1)
{
  
  dim0 <- dim(data)
  dim1 <- (dim0-1)*dx/dx1 + 1

  ii <- seq( from=0,by=dx1/dx,length.out=dim1[1]) + 1
  jj <- seq( from=0,by=dx1/dx,length.out=dim1[2]) + 1
  
  # Reduced grid to filtered points
  data1 <- data[ii,jj]

  return(data1)
} 

extend.hydrobasins <- function(grid,basins)
{  ## Calculate the hydrological basins for the whole grid
    
  # Store original data
  bas0  <- bas1 <- basins
  
  nx <- dim(grid$xx)[1]
  ny <- dim(grid$xx)[2]
  
  # Eliminate non-existent basins
  bas0[bas0==0] <- NA
  
  # Loop over all points and fill in basins
  # for points lacking a label (with the nearest neighbor)
  k <- 0
  n <- nx*ny
  
  for ( i in 1:nx ) {
    for (j in 1:ny ) {
      
      if (is.na(bas0[i,j])) {
        
        p0 <- list(x=grid$xx[i,j],y=grid$yy[i,j])
        p1 <- list(x=as.vector(grid$xx),y=as.vector(grid$yy))
        
        dists <- calcdist(p0,p1)
        
        qq <- order(dists)
        q <- qq[which(!is.na(bas0[qq]))][1]
        
        bas1[i,j] <- bas0[q]
      }
     
    }
    
    cat("Row", i,"of",nx,"\n")
  }

  # Return new basins
  return(bas1)

}

rmse <- function(x,ii=c(1:length(x)),round=NA)
{
  # Filter for desired points
  x <- x[ii]
  
  ii1 <- which(!is.na(x))
  rmse <- sqrt(sum(x[ii1]^2)/length(x[ii1]))
  
  if (!is.na(round)) rmse <- round(rmse,round)
  
  return(rmse)
  
}

rsquared <- function(x,xfit,ii=c(1:length(x)),round=NA)
{ # Calculate the R-squared value: R-squared = 1 - MSE/VAR(Y)

  # Filter for desired points
  ii0 <- c(1:length(x))
  ii1 <- which(!is.na(x) & !is.na(xfit) & ii0 %in% ii)
  x <- x[ii1]
  xfit <- xfit[ii1]
  
  xm <- mean(x)
  
  varx <- sum((x-xm)^2)
  
  MSE <- sum((xfit-x)^2)
  
  rsq <- 1 - MSE/varx
  
  if (!is.na(round)) rsq <- round(rsq,round)
  
  return(rsq)
  
}

## GEOID
geo2zs <- function(geoz,lat=NA)
{    # Convert from geopotential height to elevation
     # obtained from: http://www.ofcm.gov/fmh3/text/appendd.htm
  
  dd <- dim(geoz)
  
  G    <- 9.80665       # Gravitational constant [m/s2]
  
  if ( !is.na(lat[1]) ) {
    a    <- 6378.137e3    # Equatorial radius       [m]
    b    <- 6356.7523e3   # Polar radius            [m]
    lat1 <- lat * pi/180  # Latitude converted to radians
    
    # Radius of earth at given latitude [m]
    R_e  <- sqrt( ( (a^2*cos(lat1))^2 + (b^2*sin(lat1))^2 )  /
                  ( (a*cos(lat1))^2   + (b*sin(lat1))^2 )  )

    # Gravity at latitude [m/s2]
    g_f  <- G * ( 1 - 0.0026370*cos(2*lat1) + 0.0000059*cos(2*lat1)^2 )
                      
                    
    # Gravity ratio [--]
    Gr   <- g_f * R_e / G

    # Get the geometric elevation
    z    <- (geoz * R_e) / (Gr-geoz)
    
    #cat(lat, R_e/1d3, geoz, z,"\n")
  
  } else {
    
    # Divide geopotential height by gravitational constant
    z    <- geoz / G
    
  }
  
  ##
  ##z <- g_f
  ##
  
  dim(z) <- dd
  
  return(z)

}

## Earth grid weighting ##
## AREA of grid boxes
gridarea <- function(lon,lat,Re=6371,Dim=TRUE)
{
  # Take from: http://map.nasa.gov/GEOS_CHEM_f90toHTML/html_code/src/grid_mod.f.html
  #            2*PI*Re^2    {                                     }
  #    Area = ----------- * { sin( YEDGE[J+1] ) - sin( YEDGE[J] ) }
  #              IIGLOB     {                                     }
  # Re = radius of earth
  # YEDGE = position of latitude grid box edge
  # IIGLOB = number of longitudes
  
  nx <- length(lon)
  ny <- length(lat)
  
  # Convert to radians
  latr <- lat*pi/180
  lonr <- lon*pi/180
  
  dx <- (lonr[2] - lonr[1])/2
  dy <- (latr[1] - latr[2])/2
  
  area <- array(NA,dim=c(nx,ny))
  
  for ( j in 1:ny ) {
    a = (2*pi*Re^2/nx) * ( sin( latr[j]+dy ) - sin( latr[j]-dy ) )
    area[,j] = rep(a,nx)
  }
  
  # Dim's weighting...
  if (Dim) {
    a <- cos(latr)
    area <- matrix(rep(a,nx),byrow=TRUE,ncol=ny,nrow=nx)
  }
  
  return(area)
}

mean.areawt <- function(var,area,...)
{
  ii <- which(!is.na(var) )
  
#   wt <- area*NA
#   wt <- area[ii]/sum(area[ii])
  
  ave <- sum(var[ii]*area[ii]/sum(area[ii]))
  
  return(ave)
}

monthly2daily.smooth <- function(Tm,day=c(1:360),method="linear")
{
  Tm <- as.numeric(Tm[tolower(month.abb)])
  d0 <- c(1:12)*30 - 15
  
  if ( method == "linear" ) {
    Td <- approx(d0,Tm,day,rule=2)
  } else if ( method == "spline" ) {
    Td <- spline(d0,Tm,xout=day,method="periodic")
  }
  
  return(Td$y)
}

monthly2daily <- function(Tm,day=c(1:360))
{ ## Interpolate data in time: monthly => daily
  
  Tm <- as.numeric(Tm)
  
  ndm = 30
  nm  = 12
  
  n <- length(day)
  m0  = m1  = numeric(n)
  wt0 = wt1 = numeric(n)
  Td        = numeric(n)
  
  # Get length of arrays to fill, midpoint of month (in days)
  # and the total weight (total number of days in a month)
  mid = ndm / 2
   
  for (k in 1:n) {
    
    d <- day[k]
    
    for( m in 1:(nm+1)) {
      if ( m*ndm-mid > d ) break
    }
    m1[k] = m; m0[k] = m1[k]-1
    
    if ( m1[k] > 12 ) m1[k] =  1
    if ( m0[k] <  1 ) m0[k] = 12
    
#     wt1[k] = abs( mod(k-mid,ndm) )
    wt1[k] = abs( (d-mid) %% ndm )
    wt0[k] = ndm - wt1[k]
    
    wttot = wt0[k] + wt1[k]
    wt0[k] = wt0[k]/wttot; wt1[k] = wt1[k]/wttot
    
    Td[k] = wt0[k]*Tm[m0[k]] + wt1[k]*Tm[m1[k]]
  }
   
  return(Td)
  
}

effectiveT <- function(T,sigma=5)
{
  inv_sqrt2   = 1.0/sqrt(2.0)
  inv_sqrt2pi = 1.0/sqrt(2.0*pi)

  inv_sigma   = 1.0/sigma

  Teff = sigma*inv_sqrt2pi*exp(-0.5*(T*inv_sigma)^2) +
               T*0.5*erfcc(-T*inv_sigma*inv_sqrt2)

  return(Teff)
}

erfcc <- function(x)
{
  z = abs(x)
  t = 1.0/(1.0+0.5*z)

  erfc = t*exp(-z*z-1.265512230+t*(1.000023680+t*(0.374091960 +
    t*(0.096784180+t*(-0.186288060+t*(0.278868070+t*
    (-1.135203980+t*(1.488515870+t*(-0.822152230+
    t*0.170872770)))))))))

  #if (x < 0.0) { erfcc = 2.0-erfcc }
  erfc[x < 0.0] = 2.0-erfc[x < 0.0]

  return(erfc)
}

### RCM data loading

dma2cma <- function(dma,area=1.74,conv=1e-2*area*1e6*1e-6)
{ # dma given in percent extent
  # convert percent to fraction
  # convert area to km2 from 1e6 km2
  # convert final cma from km2 to 1e6 km2
  return(sum(dma)*conv)
}

get_cma_series <- function(dat,nm="MAR_MELT1B",mnm="MSK",mice=1,t0=1958,dx=25,nd=365,area=1.74)
{

  mask <- dat[[mnm]]; ii <- which(mask!=mice)
  n <- dim(dat[[nm]])[3]
  
  t <- seq(from=t0,by=1,length.out=n)
  
  ma <- array(0,dim=dim(dat[[nm]]))
  cma <- ama <- numeric(n)
  
  # Loop over each year to output the filtered mask and the CMA
  for ( j in 1:n ) {
    
    melt <- dat[[nm]][,,j]; melt[ii] <- 0
    
    # Store the array of melt area for this year
    ma[,,j] <- melt
    
    # Now determine CMA for this year
    cma[j] <- sum(melt)*(dx)^2 * 1e-6
    
    ama[j] <- 100* (cma[j] / nd) / area
  }
  
  return(list(time=t,cma=cma,ama=ama,ma=ma))
}

Gt2sl <- 1/362   # 1mm = 362 Gt

aline <- function(time=NA,var=NA,t0=1979,tf=2050,y0=0.7556875,slope=36e4*1e-6,sd=17e4*1e-6/3)
{
  # Get the fit here
  if (!is.na(time[1])) {
    ii <- which(time >= t0 & time <= tf)
    x00 <- time[ii]; y00 <- var[ii]
    fit <- lm(y00~x00)
    
    sd <- sd(fit$residuals)
    y0 <- fit$coefficients[1]; slope <- fit$coefficients[2]
    y0 <- y0 + t0*slope
  }
  # fit obtained
  
  # Get fit and confidence intervals
  t <- seq(from=t0,to=tf,by=1)
  y <- predict(fit,data.frame(x00=t),interval="confidence")
  y1 <- y[,2]; y2 <- y[,3]
  y <- y[,1]
  
#   y <- y0 + (t-t0)*slope
#   
#   y1 <- y0 + (t-t0)*(slope-sd)
#   y2 <- y0 + (t-t0)*(slope+sd)
  
  return(list(slope=slope,sd=sd,time=t,y=y,y1=y1,y2=y2))
}

load.rcmdata <- function(file="../rcms/RCMGCM_mar.era40.txt",
                             time=NA,trange=c(1900,2500),trange.base=c(1958,2001))
{ ## Load temperature time-series and add to datasets
  
  # Read the file and extract the average of the months of interest
  tmp <- read.table(file,header=TRUE)
  names(tmp)[2:13] <- toupper(month.abb)
  tmp$tjja <- apply(tmp[,c("JUN","JUL","AUG")],FUN=mean,MARGIN=c(1))
  
  n <- length(tmp$DEC)
  tmp$DECm1 <- tmp$DEC
  tmp$DECm1[2:n] <- tmp$DEC[1:(n-1)]
  
  tmp$tdjf  <- apply(tmp[c("DECm1","JAN","FEB")],FUN=mean,MARGIN=1)
  tmp$tdjfm <- apply(tmp[c("DECm1","JAN","FEB","MAR")],FUN=mean,MARGIN=1)
  
  # Filter to the time range of interest
  ii <- which(tmp$time >= trange[1] & tmp$time <= trange[2])
  if (!is.na(time[1])) ii <- which(tmp$time %in% time)
  
  out <- as.list(tmp[ii,])
  
  # Determine the base temperature during the base time period and
  # subtract to get anomalies
  ii <- which(out$time >= trange.base[1] & out$time <= trange.base[2])
  
  ## WRONG - I needed to add ii here (18.05.2011)
  out$tjja0 <- mean(out$tjja[ii],na.rm=TRUE)
  out$dtjja <- out$tjja - out$tjja0
  
  out$tdjf0 <- mean(out$tdjf[ii],na.rm=TRUE)
  out$dtdjf <- out$tdjf - out$tdjf0
  
  if ( "snow" %in% names(out) ) out$precip <- out$snow + out$rain
  
  return(out)
}


## Load and process RCM data
save.rcmdata <- function()
{
  
  # Load RCM daily data
  dmas <- read.table("../rcms/Fettweis_dma_1958-2001.txt",header=TRUE)
  dmas$Day <- dmas$Day*360/365   # Scale days to 360-day year
  
  dmas2 <- read.table("../rcms/Fettweis_dma_1979-2009.txt",header=TRUE)
  dmas2$Day <- dmas2$Day*360/365   # Scale days to 360-day year_offset
  
  ## Get CMA from daily dataset
  cma.thresh <- c(7.75,8.50,9.25)
  cma.mar <- c(dma2cma(dmas$MAR1a),dma2cma(dmas$MAR1b),dma2cma(dmas$MAR1c))
  cma.rac <- c(dma2cma(dmas$RACMO1a),dma2cma(dmas$RACMO1b),dma2cma(dmas$RACMO1c))
  
  nd <- 360; conv <- 100/1.74
  ama.mar <- conv* cma.mar / nd
  ama.rac <- conv* cma.rac / nd
  
  # Load yearly melt data (2D)
  dat <- get.nc("../rcms/MAR-RACMO2-SAT_melt_extent.nc")
  
  # MAR
  mar1b <- load.rcmdata(file="../rcms/RCMGCM_mar.era40.txt",trange.base=c(1958,2001))                         
  tmp   <- get_cma_series(dat,nm="MAR_MELT1B",mnm="MSK",mice=1,t0=1958,dx=25) 
  ii    <- which(tmp$time %in% mar1b$time)
  ii0 <- which(mar1b$time %in% tmp$time[ii])
  mar1b$cma <- mar1b$ama <- numeric(length(mar1b$time))*NA
  mar1b$ma <- array(NA,dim=c(dim(tmp$ma)[1:2],length(mar1b$time)))
  mar1b$cma[ii0]  <- tmp$cma[ii]
  mar1b$ama[ii0]  <- tmp$ama[ii]
  mar1b$ma[,,ii0] <- tmp$ma[,,ii]
  
  # RACMO
  rac1b <- load.rcmdata(file="../rcms/RCMGCM_racmo.era40.txt",trange.base=c(1958,2001))  
  tmp   <- get_cma_series(dat,nm="RAC_MELT1B",mnm="MSK",mice=1,t0=1958,dx=25)
  ii    <- which(tmp$time %in% rac1b$time)
  ii0 <- which(rac1b$time %in% tmp$time[ii])
  rac1b$cma <- rac1b$ama <- numeric(length(rac1b$time))*NA
  rac1b$ma <- array(NA,dim=c(dim(tmp$ma)[1:2],length(rac1b$time)))
  rac1b$cma[ii0]  <- tmp$cma[ii]
  rac1b$ama[ii0]  <- tmp$ama[ii]
  rac1b$ma[,,ii0] <- tmp$ma[,,ii]
  
  # SAT
  sat1b <- get_cma_series(dat,nm="SAT_MELT1B",mnm="MSK",mice=1,t0=1958,dx=25)
  # Make an artificial temperature series for satellite data
  ii1 <- which(mar1b$time %in% sat1b$time)
  ii2 <- which(rac1b$time %in% sat1b$time)
  
  nms <- c(toupper(month.abb),"tjja","tdjf","tdjfm")
  for (q in 1:length(nms)) {
    nm <- nms[q]
    sat1b[[nm]] <- apply( rbind(mar1b[[nm]][ii1],rac1b[[nm]][ii2]),FUN=mean,MARGIN=c(2) )
  }
  sat1b$tjja0 <- mean( sat1b$tjja[ which(sat1b$time >= 1958 & sat1b$time <= 2001) ] )
  sat1b$dtjja <- sat1b$tjja - sat1b$tjja0
  
  # Get base tjja for 1980-1999 like GCMs
  t0 <- 1980; t1 <- 1999
  ii <- which(mar1b$time >= t0 & mar1b$time <= t1)
  mar1b$tjja0gcm <- mean(mar1b$tjja[ii],na.rm=TRUE)
  ii <- which(rac1b$time >= t0 & rac1b$time <= t1)
  rac1b$tjja0gcm <- mean(rac1b$tjja[ii],na.rm=TRUE)
  ii <- which(sat1b$time >= t0 & sat1b$time <= t1)
  sat1b$tjja0gcm <- mean(sat1b$tjja[ii])
  dtjja0 <- mean(sat1b$dtjja[ii],na.rm=TRUE)
  
  # Add dummy field to mimic gcm fields
  mar1b$dtjja0 <- mar1b$dtjja
  rac1b$dtjja0 <- rac1b$dtjja
  sat1b$dtjja0 <- sat1b$dtjja
  
  ## Load GCM cmas and temperatures...
  
  # CMA data
  datcma <- read.table("../rcms/RCMGCM_proj.cma.txt",header=TRUE)
  datcma[,2:dim(datcma)[2]] <- datcma[,2:dim(datcma)[2]]*1e1
  nd.gcm <- 365
  
  # MAR - ECHAM5 - A1B
  MAR.echam5.A1B <- load.rcmdata(file="../rcms/RCMGCM_mar.echam5.A1B.txt",trange.base=c(1980,1999)) 
  ii <- which(datcma$time %in% MAR.echam5.A1B$time)
  MAR.echam5.A1B$cma <- datcma$MAR.echam5.A1B[ii]                      
  MAR.echam5.A1B$ama <- conv* MAR.echam5.A1B$cma / nd.gcm
  
  # Offset the anomaly dtjja to account for bias in GCM
  ii <- which(MAR.echam5.A1B$time >= t0 & MAR.echam5.A1B$time <= t1)
  dtjja1 <- mean(MAR.echam5.A1B$dtjja[ii])
  MAR.echam5.A1B$dtjja.offset <- dtjja1-dtjja0
  MAR.echam5.A1B$dtjja0 <- MAR.echam5.A1B$dtjja
  MAR.echam5.A1B$dtjja <- MAR.echam5.A1B$dtjja0 - MAR.echam5.A1B$dtjja.offset
  
  
  # MAR - ECHAM5 - E1
  MAR.echam5.E1 <- load.rcmdata(file="../rcms/RCMGCM_mar.echam5.E1.txt",trange.base=c(1980,1999)) 
  ii <- which(datcma$time %in% MAR.echam5.E1$time)
  MAR.echam5.E1$cma <- datcma$MAR.echam5.E1[ii]
  MAR.echam5.E1$ama <- conv* MAR.echam5.E1$cma / nd.gcm
  
  # Offset the anomaly dtjja to account for bias in GCM
  ii <- which(MAR.echam5.E1$time >= t0 & MAR.echam5.E1$time <= t1)
  dtjja1 <- mean(MAR.echam5.E1$dtjja[ii])
  MAR.echam5.E1$dtjja.offset <- dtjja1-dtjja0
  MAR.echam5.E1$dtjja0 <- MAR.echam5.E1$dtjja
  MAR.echam5.E1$dtjja <- MAR.echam5.E1$dtjja0 - MAR.echam5.E1$dtjja.offset
  
  
  # MAR - HADCM3 - A1B   
  MAR.hadcm3.A1B <- load.rcmdata(file="../rcms/RCMGCM_mar.hadcm3.A1B.txt",trange.base=c(1980,1999)) 
  ii <- which(datcma$time %in% MAR.hadcm3.A1B$time)
  MAR.hadcm3.A1B$cma <- datcma$MAR.hadcm3.A1B[ii]
  MAR.hadcm3.A1B$ama <- conv* MAR.hadcm3.A1B$cma / nd.gcm
  
  # Offset the anomaly dtjja to account for bias in GCM
  ii <- which(MAR.hadcm3.A1B$time >= t0 & MAR.hadcm3.A1B$time <= t1)
  dtjja1 <- mean(MAR.hadcm3.A1B$dtjja[ii])
  MAR.hadcm3.A1B$dtjja.offset <- dtjja1-dtjja0
  MAR.hadcm3.A1B$dtjja0 <- MAR.hadcm3.A1B$dtjja
  MAR.hadcm3.A1B$dtjja <- MAR.hadcm3.A1B$dtjja0 - MAR.hadcm3.A1B$dtjja.offset
   
  # Calculate projections
  pmar1b <- aline(time=mar1b$time,var=mar1b$cma,y0=NA,t0=1979,tf=2050)
  psat1b <- aline(time=sat1b$time,var=sat1b$cma,y0=NA,t0=1979,tf=2050)
  prac1b <- aline(time=rac1b$time,var=rac1b$cma,y0=NA,t0=1979,tf=2050)
  
  save(dmas, dmas2, cma.thresh, cma.mar, cma.rac, ama.mar, ama.rac,
       mar1b, rac1b, sat1b, MAR.echam5.A1B, MAR.echam5.E1,MAR.hadcm3.A1B,
       file="rcmdata.Rdata")
  
}

## NEW ##
load_MAR0 <- function(fldr="data/rcms/MARv2/MARv2_ERA-40_1958-1999")
{
  
  vnms = c("ice","me","melt","rf","ru","sf","stt","su","tt")
  files = list.files(fldr)
  
  if (length(vnms) != length(files) ) {
    cat("Wrong files. Fix function 'load_MAR'!\n")
    cat(vnms,"\n")
    cat(files,"\n")
  }
  
  dat = list()
  for (q in 1:length(vnms)) {
    vnm = vnms[q]
    q0 = grep(paste("-",vnm,".dat",sep=""),files)

    if ( length(q0) > 0 ) {

      tmp = read.table(file.path(fldr,files[q0]))

      if (q==1) dat$time = tmp[[1]]
      dat[[vnm]] = tmp[,c(2:13)]
      names(dat[[vnm]]) = tolower(month.abb)
    }

  }
  
  cat("Loaded MAR folder:",fldr,"\n")

  return(dat)
}

load_MAR <- function(fldr,time=c(1950:2100),time.base=c(1980:2000))
{
  icearea = 1744980  # km2

  # Now load MAR data set(s)
  fldr1 = NULL
  if (length(grep("rcp",fldr))>0) {
    fldr1 = fldr
    fldr1 = gsub("rcp26_2006-2100","histo_1965-2005",fldr1)
    fldr1 = gsub("rcp45_2006-2100","histo_1965-2005",fldr1)
    fldr1 = gsub("rcp85_2006-2100","histo_1965-2005",fldr1)
  }
  if (length(grep("ERA-40",fldr))>0) {
    fldr1 = fldr
    fldr1 = gsub("40_1958-1999","INTERIM_1979-2011",fldr1)
  }
  
  dat0 = list(load_MAR0(fldr))
  if (!is.null(fldr1)) dat0[[2]] = load_MAR0(fldr1)
  
  # dat = list(time=time,month=c(1:12),tjja=NA,dtjja=NA,cma=NA)
  dat = list(time=time,month=c(1:12))
  nt = length(time)
  dat$tt   = as.data.frame(array(NA,dim=c(nt,12)))
  dat$melt = as.data.frame(array(NA,dim=c(nt,12)))
  dat$smb  = as.data.frame(array(NA,dim=c(nt,12)))

  for (q in 1:length(dat0)) {
    
    tmp = dat0[[q]]
    
    # Filter to the time range of interest
    ii0 <- which(tmp$time %in% dat$time)
    ii1 <- which(dat$time %in% tmp$time)

    dat$tt[ii1,]   = tmp$tt[ii0,]
    dat$melt[ii1,] = tmp$melt[ii0,]
    dat$smb[ii1,]  = tmp$sf[ii0,] + tmp$rf[ii0,] - tmp$ru[ii0,] - tmp$su[ii0,]

  }
  
  dat$tjja = apply(dat$tt[,c(6,7,8)],FUN=mean,MARGIN=c(1))
  dat$cma  = apply(dat$melt*icearea,FUN=sum,MARGIN=c(1))
  dat$cma[dat$cma==0] = NA

  ## Get temp anomalies relative to base period
  ii = which(dat$time %in% time.base)
  T0 = mean(dat$tjja[ii])
  dat$dtjja = dat$tjja - T0
  
  ## Get REMBO base period if ERA-40
  if (length(grep("ERA-40",fldr))>0) {
    ii = which(dat$time %in% c(1958:2001)) 
    T0 = mean(dat$tjja[ii])
    dat$dtjja.rembo = dat$tjja - T0
  }
  
  return(dat)
}

###

### GRIDDING
normalize <- function(x,sd1=NA,ave=NA)
{
  if (is.na(sd1)) sd1 <- sd(x,na.rm=TRUE)
  if (is.na(ave)) ave <- mean(x,na.rm=TRUE)
  
  x1 <- (x-ave)/sd1
  return(x1)
}

standardize <- function(pres,ibase=c(1:length(pres)))
{ # Standardize the data of a given month
  
  # Determine the average pressure during the base time period
  # Then calculate the anomalies from this average pressure
  ave <- mean(pres[ibase],na.rm=TRUE)
  dpres <- pres - ave
  
  # Calculate the standard deviation of the anomalies
  sd1 <- sd(dpres[ibase],na.rm=TRUE)
  
  # Divide by the standard deviation of the anomalies to normalize (standardize) it
  stdpres <- dpres / sd1
  
  return(stdpres)
}

findpoint <- function(dat,p=c(lat=37.82,lon=-25.73),ebs=3)
{ # Find the nearest grid point corresponding to a lat/lon location

  ii.point <- which( dat$lat >= (p[1]-ebs) &
                     dat$lat <= (p[1]+ebs) &
                     dat$lon >= (p[2]-ebs) &
                     dat$lon <= (p[2]+ebs) )
  
  d1 <- 1e6
  for (i in ii.point) {
    d <- calcdistw(p1=data.frame(lat=p[1],lon=p[2]),
                   p2=data.frame(lat=dat$lat[i],lon=dat$lon[i]))
    if ( d < d1 ) { d1 <- d; i1 <- i}
  }
  
  return(i1)
}

## INSOLATION LOADING ##
get_insol <- function(time,lat,day_mean=c(171))
{
  nc = open.ncdf("~/wrk_local/data/paleo/insolation.NH.150ka.nc")

  nc.time = get.var.ncdf(nc,"time")*1e-3
  nc.lat  = get.var.ncdf(nc,"latitude")

  k0 = which.min(abs(nc.time-min(time)))
  k1 = which.min(abs(nc.time-max(time)))
  k0 = max(0,k0-1)
  k1 = min(length(nc.time),k1+1)
  nk = k1-k0+1

  l0 = which.min(abs(nc.lat-min(lat)))
  l1 = which.min(abs(nc.lat-max(lat)))
  l0 = max(0,l0-1)
  l1 = min(length(nc.lat),l1+1)
  nl = l1-l0+1

  lat1  = nc.lat[c(l0:l1)]
  time1 = nc.time[c(k0:k1)]

  ins3D = get.var.ncdf(nc,"ins",start=c(1,l0,k0),count=c(-1,nl,nk))
  close.ncdf(nc)
  
  if (length(day_mean)>1) {
    ins2D = apply(ins3D[day_mean,,],c(2,3),mean)
  } else {
    ins2D = ins3D[day_mean,,]
  }

  ## Now interpolate to desired times and lats 
  obj = list(x=lat1,y=time1,z=ins2D)

  lat0 = lat 
  if (!is.null(dim(lat0))) lat = as.vector(lat0)

  loc = expand.grid(x=lat,y=time)
  new = interp.surface(obj,loc)

  if (!is.null(dim(lat0))) {
    dim(new) = c(dim(lat0),length(time))
  } else {
    dim(new) = c(length(lat0),length(time))
    if (length(lat0)==1) new = as.vector(new)
  }

  return(new)
}

## LOADING BINARY (GRADS) FILES
# ----------------------------------------------------------------------
# Function: get.ctl
# Author  : Alex Robinson
# Purpose : Get info from a ctl file
# ----------------------------------------------------------------------
get.ctl <- function(fnm="clima.ctl")
{

  ctl <- system(paste("/home/robinson/python/sico/embfunctions.py",fnm),intern=TRUE)
  ctl <- strsplit(ctl," ")[[1]]

  gx = ctl[1]
  nx = as.integer(ctl[2])
  x0 = as.double(ctl[3])
  dx = as.double(ctl[4])
  ny = as.integer(ctl[5])
  y0 = as.double(ctl[6])
  dy = as.double(ctl[7])

  nvar = as.integer(ctl[8])
  vnms = ctl[9:length(ctl)]

  x = seq(from=x0,to=x0+(nx-1)*dx,by=dx)
  y = seq(from=y0,to=y0+(ny-1)*dy,by=dy)

  return(list(gx=gx,vnms=vnms,nvar=nvar,nxny=c(nx,ny),x=x,y=y))
}

# ----------------------------------------------------------------------
# Function: get.binary
# Author  : Alex Robinson
# Purpose : Get data from a binary output file
# ----------------------------------------------------------------------
get.binary <- function(fldr=".",gx="cntro.clima.gx",ctl="na",nx=151,ny=281,nk=1)
{

  if (ctl != "na") {
    ctl  <- get.ctl(file.path(fldr,ctl))

    gx   <- ctl$gx
    nvar <- ctl$nvar
    nx   <- ctl$nxny[1]
    ny   <- ctl$nxny[2]

    vnms <- ctl$vnms

  } else {
    nvar = length(vnms)
  }

  ndat = nx*ny
  newdata <- file(file.path(fldr,gx), "rb")
  datavals <- readBin(newdata, real(), size = 4, n = nvar*ndat*nk, endian = "little")
  close(newdata)

  cat("\n Loaded: ",file.path(fldr,gx),"\n")


  # Currently nk is not handled, also output should be in an array
  #var <- array(datvals,c(nx,ny,nk,nvar)

  var <- list()

  if (nk == 1) {

    # Obtain list of 2d variables
    for (q in 1:nvar)
    {
      i0 <- (q-1)*ndat + 1
      i1 <- i0 + ndat - 1
      var[[q]] <- matrix(datavals[i0:i1],nrow=nx,ncol=ny)
    }

    names(var) <- vnms

  } else {

    # Obtain list through time of list of 2d variables
    for (k in 1:nk) {

      var[[k]] <- list()

      for (q in 1:nvar)
      {
        i0 <- (k-1)*ndat*nvar + (q-1)*ndat + 1
        i1 <- i0 + ndat - 1
        var[[k]][[q]] <- matrix(datavals[i0:i1],nrow=nx,ncol=ny)
      }

      names(var[[k]]) <- vnms
    }

  }

  return(var)
}

## Load 1d binary file
get.binary2 <- function(fldr="/scratch/01/andrey/PPSCM",ctl="grads/history.ctl",file="OUT/history.dat",n=3)
{
  ## Load variable names from CTL file
  ctl <- scan(file.path(fldr,ctl),what="character",sep="\n")
  ctl <- strsplit(ctl,split=" ")
  nms <- vector(mode="character")
  vars <- FALSE
  for ( i in 1:length(ctl)) {
    s <- ctl[[i]][1]
    if ( s == "ENDVARS" )vars <- FALSE
    if ( vars ) nms <- c(nms,s)
    if ( s == "VARS" ) vars <- TRUE
  }
  
  # Determine total number of variables
  nvar <- length(nms)
  
  ## Done loading variable names
  
  newdata <- file(file.path(fldr,file), "rb")
  datavals <- readBin(newdata, real(), size = 4, n = nvar*n, endian = "little")
  close(newdata)
  
  dat <- as.data.frame(t(array(datavals,dim=c(nvar,n))))
  names(dat) <- nms
  
  cat("\n Loaded: ",file.path(fldr,file),"\n")
  
  return(dat)
}

## Load 1d/2d/3d binary file
get.binary3 <- function(fldr="/scratch/01/andrey/PPSCM",ctl="grads/history.ctl",file="OUT/history.dat",
                        n=3,nx=NA,ny=NA,mynms=c("ma","zs","zb","b"))
{
  ## Load variable names from CTL file
  ctl <- scan(file.path(fldr,ctl),what="character",sep="\n")
  ctl <- strsplit(ctl,split=" ")
  nms <- vector(mode="character")
  vars <- FALSE
  for ( i in 1:length(ctl)) {
    s <- ctl[[i]][1]
    if ( s == "ENDVARS" )vars <- FALSE
    if ( vars ) nms <- c(nms,s)
    if ( s == "VARS" ) vars <- TRUE
  }
  
  # Determine total number of variables and data fields
  nvar <- length(nms)
  nt   <- n
  ntot <- prod(c(nx,ny,nvar),na.rm=T)

  ## Done loading variable names
  
  newdata <- file(file.path(fldr,file), "rb")
  datavals <- readBin(newdata, real(), size = 4, n = nt*ntot, endian = "little")
  close(newdata)
  
  if ( is.na(nx) ) {  # 0D (time series)
  
    dat <- as.data.frame(t(array(datavals,dim=c(nvar,nt))))
    names(dat) <- nms

  } else {    # Multi-dimensional!!
    
    dat   <- list()
    empty <- array(NA,dim=c(nt,nx,ny))
    
    tmp <- array(datavals,dim=c(nx,ny,nvar,nt))
    
    for ( q in 1:nvar) {
      nm <- nms[q]

      if ( nm %in% mynms ) {
        dat[[nm]] <- empty
      for (k in 1:nt) dat[[nm]][k,,] <- tmp[,,q,k]
      }
      
    }
    
    #names(dat) <- nms

    # Add lat/lon
    dat$lat <- seq(from=  21,by=0.75,length.out=ny)
    dat$lon <- seq(from=-180,by=1.5, length.out=nx)
    
  }

  cat("\n Loaded: ",file.path(fldr,file),"\n")
  
  return(dat)
}

## Load two binary files (from climber-sicopolis)
get.binaries <- function(fldr="/scratch/01/andrey/PPSCM",n,t0=0,accel=2,x1=NA,twoD=FALSE,x2=NA)
{
  
  ndat  <- n*1000/accel
  ndat2 <- n*1000
  n2d   <- n

  time1 <- c(1:ndat)*1e-3*accel+t0
  time2 <- c(1:ndat2)*1e-3 + t0
  time2d <- c(1:n2d) + t0

  if (is.na(x1[1])) x1 <- time1
  
  datclim <- get.binary2(fldr=fldr,ctl="grads/history.ctl",file="OUT/history.dat",n=ndat)
  datice  <- get.binary2(fldr=fldr,ctl="grads/time.ctl",file="sico_out/0d.ser",n=ndat2)
  
  dat2d <- NA
  if (twoD) dat2d   <- get.binary3(fldr=fldr,ctl="grads/2d_ser.ctl",file="sico_out/2d.ser",n=n2d,nx=241,ny=87)

  # Filter
  datclim <- sim.filter(datclim,x1=x1,x0=time1)
  datice  <- sim.filter(datice,x1=x1,x0=time2)
  
  tmp <- dat2d
  x0 <- c(1:dim(dat2d$zs)[1]) + t0
  kk <- which(x0 %in% x2)

  nms <- names(dat2d)
  for ( q in 1:length(nms) ) {
    nm <- nms[q]
    if ( ! nm %in% c("lat","lon") ) {
      tmp[[nm]] <- dat2d[[nm]][kk,,]
    }
  }
  dat2d      <- tmp
  dat2d$time <- x0[kk]

  # Combine series
  datice$time <- NULL
  dat <- cbind(datclim,datice)
  
  return(list(series=dat,dat2d=dat2d))
}



## HYSTERESIS
get.basins <- function(infob,xhyst,yhyst) 
{
  
  # Get number of cases
  n <- dim(infob)[1]
  
  # Get the x values and generate y values to 
  # make up the basins image
  x <- sort( unique(infob$T_warming) )
  y <- sort( unique(infob$Vtot.init) )

  # Determine current max/min volumes of bounding hyst curves
  Vmins <- Vmaxs <- numeric(length(x))
  ebs <- 0.05
  for ( i in 1:length(x) ) {
    Tnow <- x[i]
    q <- which(abs(xhyst[,1]-Tnow) < ebs)
    Vmins[i] <- mean(yhyst[q,1],na.rm=TRUE)
    q <- which(abs(xhyst[,2]-Tnow) < ebs)
    Vmaxs[i] <- mean(yhyst[q,2],na.rm=TRUE)
  }
  
  # Create the basins object, x=temperature,y=volume,z=final volume
  basins <- list(x=x,y=y,z=matrix(0,nrow=length(x),ncol=length(y)))
  basins$time <- basins$z0 <- basins$zf <- basins$z
  
  Vscale <- max(infob$Vtot.eq,na.rm=TRUE)
  
  for ( i in 1:n ) {
    
    # Get current values
    Tnow  <- infob$T_warming[i]
    Veq   <- infob$Vtot.eq[i]
    Vinit <- infob$Vtot.init[i]
    
    # Get x and y indices of current case
    ii <- which( x == Tnow )
    jj <- which( y == Vinit )
    
    # Get bounding volumes for current dT
    Vmin <- Vmins[ii]; Vmax <- Vmaxs[ii]
    Vscale <- Vmax
    
    # Determine basin value (-1,-2: ambigous, 0: upper basin, 1: lower basin)
    Bval <- 0
    if ( !is.na(Veq) ) {
      
      # Check for ambiguous points
      # (points that are in the middle, higher or lower than starting point...)
      if ( abs(Veq-Vmin)/Vscale > 0.15 & abs(Veq-Vmax)/Vscale > 0.15 ) Bval <- -1
      
      # Check for ambiguous points that are lower than starting point
      # (points that have reduced volume, but did not go to Vmin)
      #if ( (Veq-Vinit)/Vscale < -0.05 & (Veq-Vmin)/Vscale > 0.05 ) Bval <- -2
      
      # Check for points that are close to the minimum equil. value
      if ( abs(Veq-Vmin)/Vscale < 0.1 ) Bval <- 1
      
      # If Vmin and Vmax are converging, then basins dont exist
      if ( abs(Vmax-Vmin)/Vscale < 0.5 ) Bval <- 0
    }
        
    # Store values in basins list
    basins$z0[ii,jj]   <- Vinit
    basins$zf[ii,jj]   <- Veq
    basins$z[ii,jj]    <- Bval
    basins$time[ii,jj] <- infob$time.eq[i]
    
    # Add hysteresis curve to facilitate plotting
    # xhyst, yhyst, should contain two curves each, one for anf_dat=2 and anf_dat=3
    
    # Eliminate NA values
    jj <- which(!is.na(apply(cbind(xhyst,yhyst),FUN=sum,MARGIN=1)))
    xhyst <- xhyst[jj,]
    yhyst <- yhyst[jj,]
    
    basins$xhyst <- xhyst
    basins$yhyst <- yhyst
  }

  return(basins)
}

## Add basins to a plot
add.basins <- function(basins,b.pts=FALSE,b.lines=FALSE,tmax=NA,lwd=2,col=1,lwd.box=2)
{
  ccc <- colorRampPalette(c("darkred","red","yellow"))(7)
  c1 <- ccc[4]
  c1 <- rgb(t(col2rgb(c1)),alpha=230,max=255)
  c2 <- rgb(t(col2rgb(c1)),alpha=80,max=255)

  # Colors: c(ambiguous, upper basin, lower basin)
  c1 <- c(c2,"white",c1)
  
  # Plot the basins first
  image(basins,col=c1,add=TRUE)
  
  # Make horizontal lines indicating each 10% initial starting volume
  if (b.lines) {
    jj <- c(1:length(basins$y))
    yy <- basins$y[which(jj %% 2 == 0)]  #; yy<-yy[-1]
    fracs <- c("20%","30%","40%","50%","60%","70%","80%","90%")
    col <- rep("grey20",length(yy)); col[1:3] <- "grey95"
    abline(h=yy,col=col[8],lwd=2)
    x1 <- xlim[2]-(xlim[2]-xlim[1])*0.08
    for ( j in 1:length(yy) ) {
      x1 <- x1 - (xlim[2]-xlim[1])*0.06
      text(x1,yy[j]-0.08,fracs[j],cex=0.5,col=col[j])
    }
  }
  
  # Set color of hyst plot to black
  #cols <- cols*0 + 1
  
  # Get polygon outline of hyst curve, so that I can
  # color over basin values outside of the curve
  xx1 <- c(min(basins$xhyst[,1]),max(basins$xhyst[,1]),basins$xhyst[,1] )
  yy1 <- c(0,0,basins$yhyst[,1] )
  
  xx2 <- c(max(basins$xhyst[,1]),min(basins$xhyst[,1]),basins$xhyst[,2] )
  yy2 <- c(5,5,basins$yhyst[,2] )
  
  polygon(xx1,yy1,col=0,border=0)
  polygon(xx2,yy2,col=0,border=0)
  
  ## Replot hysteresis lines
  for (i in 1:2) {
    lines(basins$xhyst[,i],basins$yhyst[,i],lwd=lwd,col=col)
  }
  
  # Plot ambiguous points of basins
  if (b.pts) {
    x1 <- matrix(basins$x,nrow=length(basins$x),ncol=length(basins$y))
    
    x1 <- as.vector(x1); V1 <- as.vector(basins$zf)
    t1 <- as.vector(basins$time); if (!is.na(tmax)) t1[t1 > tmax] <- tmax
    c2 <- rep("lightgreen",length(t1)); c2[ t1 >= (0.90*max(t1)) ] <- "darkgreen"
    c2[ ! as.vector(basins$z) %in% c(-1,2) ] <- "grey70"
    ii <- which(c2=="grey70"); points(x1[ii],V1[ii],pch=3,col=c2[ii],lwd=2.5)
    ii <- which(c2!="grey70"); points(x1[ii],V1[ii],pch=3,col=c2[ii],lwd=2.5)
  }
  
  box(lwd=lwd.box)
}

### Set up hysteresis plotting function ###
plot.hyst <- function(sets,nms=NA,basins=NA,by="itm_c",v="Vtot.eq",xlab="Temperature anomaly (°C)",
                      ylab=expression(paste("Volume (million ",km^{3},")")),
                      title=NA,panel=NA,legend="topright",ylim=c(0,4),cols=NA,
                      lwd.box=1,col.box=1,col.lab=1,x.by=0.5,y.by=1,tmax=NA,
                      xlim=c(-2,4),lwd=c(1,2.5),lty=c(1,1),pch=c(21,20),cex.pch=c(0.8,1.1),
                      axes=c(TRUE,TRUE),b.pts=FALSE,...)
{
  
  # Determine how many sets there are (based on 'by')
  s <- unique(sets[,by])
  nsets <- length(s)
  
  if ( is.na(cols[1]) & length(cols)==1 ) {
    #cols <- c("slateblue4","mediumpurple1","orangered3","green","magenta")
    #cols <- c("blue","mediumpurple4","orangered2","green","magenta")
    #cols <- c("slateblue1","darkviolet","salmon")
    cols <- c("blue","black","red")
  }
  upcols <- cols  #gen.bad.cols(cols,bad=c(1,1,1),m=180)

  ylim <- range(ylim,sets[,v],na.rm=TRUE)  # Make sure volume range is at least 0,4
  
  # Determine the axes
  xx0 <- xlim[1]
  if ( xlim[1] %% 1 != 0.0 ) {
    offset <- xlim[1] %% 1
    if ( (xlim[1]+offset) %% 1 != 0.0 ) offset <- -offset
    xx0 <- xlim[1] + offset
  }
  xmajor <- seq(from=xx0,to=xlim[2],by=x.by)
  xminor <- seq(from=xx0,to=xlim[2],by=x.by/2)
  xminor[match(xmajor,xminor)] <- NA
  
  ymajor <- seq(from=ylim[1],to=ylim[2],by=y.by)
  yminor <- seq(from=ylim[1],to=ylim[2],by=y.by/2)
  yminor[match(ymajor,yminor)] <- NA
  
  # Make the initial blank plot
  plot(sets$T_warming,sets[,v],type="n",xlim=xlim,ylim=ylim,
       xlab="",ylab="",axes=FALSE)
       
  cex.lab <- 0.8
  if ( !xlab=="" ) mtext(xlab,side=1,line=1.5,cex=cex.lab,col=col.lab)
  ylab0 <- ylab; if ( is.expression(ylab) ) ylab0 <- "title exists"
  if ( !ylab0=="" ) mtext(ylab,side=2,las=0,line=1.5,cex=cex.lab,col=col.lab) 
  
  if (length(basins) >= 3) {
    
    ccc <- colorRampPalette(c("darkred","red","yellow"))(7)
    c1 <- ccc[4]
    #c1  <- colorRampPalette(c("red","black"))(10)[2]
    c1 <- rgb(t(col2rgb(c1)),alpha=230,max=255)
    c2 <- rgb(t(col2rgb(c1)),alpha=80,max=255)
    
#     
#     c1 <- ccc[4]; c2 <- ccc[6]

    # Colors: c(ambiguous, upper basin, lower basin)
    c1 <- c(c2,"white",c1)
    
    # Plot the basins first
    image(basins,col=c1,xlim=xlim,ylim=ylim,add=TRUE)
    
    # Make horizontal lines indicating each 10% initial starting volume
    if (FALSE) {
      jj <- c(1:length(basins$y))
      yy <- basins$y[which(jj %% 2 == 0)]  #; yy<-yy[-1]
      fracs <- c("20%","30%","40%","50%","60%","70%","80%","90%")
      col <- rep("grey20",length(yy)); col[1:3] <- "grey95"
      abline(h=yy,col=col[8],lwd=2)
      x1 <- xlim[2]-(xlim[2]-xlim[1])*0.08
      for ( j in 1:length(yy) ) {
        x1 <- x1 - (xlim[2]-xlim[1])*0.06
        text(x1,yy[j]-0.08,fracs[j],cex=0.5,col=col[j])
      }
    }
    
    # Set color of hyst plot to black
    #cols <- cols*0 + 1
    
    # Get polygon outline of hyst curve, so that I can
    # color over basin values outside of the curve
    #set <- sets[sets$set==1,]
    set <- sets
    i <- order(set$T_warming,decreasing=TRUE); set <- set[i,]
    
    # Get indices of up and dn values
    up <- which( set$anf_dat == 2 )
    dn <- which( set$anf_dat != 2 )
    
    xx1 <- c(seq(from=-10,to=10),   set$T_warming[up] )
    yy1 <- c(seq(from=-10,to=10)*0, set[up,v]         )
    
    xx2 <- c(seq(from=-10,to=10),   set$T_warming[dn] )
    yy2 <- c(seq(from=-10,to=10)*0 + 5, set[dn,v]     )
    
    polygon(xx1,yy1,col=0,border=0)
    polygon(xx2,yy2,col=0,border=0)
    
    # Plot ambiguous points of basins
    if (b.pts) {
      x1 <- matrix(basins$x,nrow=length(basins$x),ncol=length(basins$y))
      
      x1 <- as.vector(x1); V1 <- as.vector(basins$zf)
      t1 <- as.vector(basins$time); if (!is.na(tmax)) t1[t1 > tmax] <- tmax
      c2 <- rep("lightgreen",length(t1)); c2[ t1 >= (0.90*max(t1)) ] <- "darkgreen"
      c2[ ! as.vector(basins$z) %in% c(-1,2) ] <- "grey70"
      ii <- which(c2=="grey70"); points(x1[ii],V1[ii],pch=3,col=c2[ii],lwd=2.5)
      ii <- which(c2!="grey70"); points(x1[ii],V1[ii],pch=3,col=c2[ii],lwd=2.5)
    }
    
  } else {
    
    # Add the background grid
    #grid(lty=1)
    abline(h=ymajor,v=xmajor,lty=3,col=8)
    abline(h=yminor,v=xminor,lty=3,col=8)
    
  }

  leg.lwd=numeric(nsets)+lwd[2]
  leg.lty=numeric(nsets)+lty[2]
  leg.pch=numeric(nsets)+pch[2]
  leg.col=cols
  
  if (!is.na(panel)) { 
    x0 <- xlim[1]+0.5; y0 <- 0.3
    text(x0,y0,panel,cex=0.5,vfont=c("serif","bold"),col="grey30")
  }
  
  if (!is.na(nms[1])) {
    leg.bg <- "grey90"
    legend(legend,nms,     #title=title,
           lwd=leg.lwd,lty=leg.lty,
           #pch=leg.pch,
           col=leg.col,cex=0.6,inset=c(0.02,0.02),
           box.lwd=0,box.col=leg.bg,bg=leg.bg)
#     legend(legend,nms,     #title=title,
#            lwd=leg.lwd,lty=leg.lty,
#            #pch=leg.pch,
#            col=leg.col,cex=0.65,ncol=3,#inset=c(0.01,0.02),
#            box.lwd=0,box.col=leg.bg,bg=leg.bg)
  }
  
  
  ### FOR hyst.supplementary ### 
#     colnow <- cols
#     for (i in 1:nsets) {
#       # Get the current set
#       ii <- which( sets[,by] == s[i] )
#       set <- sets[ii,]
#       
#       lines(set$T_warming,set[,v],lwd=lwd[i],lty=lty[i],col=colnow[i])
#     }
  ###
  
  
  # Loop through ups and downs
  for (q in 2:3) {
    
    colnow <- cols
    if (q==2) colnow <- upcols
    qq <- q-1
    
    # Loop through the sets of hysteresis data
    for (i in 1:nsets) {
      
      # Get the current set
      ii <- which( sets[,by] == s[i] & sets$anf_dat == q )
      set <- sets[ii,]
      
      # Plot ups
      lines(set$T_warming,set[,v],lwd=lwd[qq],lty=lty[qq],col=colnow[i])
      #points(set$T_warming,set[,v],pch=pch[qq],cex=cex.pch[qq],col=colnow[i])

    }
  }
    
  # Add axes!
  if (axes[1]) {
    axis(1,at=xmajor,labels=TRUE,...)
    #axis(1,at=xminor,labels=FALSE)
    #axis(3,at=xmajor,labels=FALSE)
    #axis(3,at=xminor,labels=FALSE)
  }
  if (axes[2]) {
    axis(2,at=ymajor,labels=TRUE,...)
    #axis(4,at=ymajor,labels=TRUE,...)
    #axis(2,at=yminor,labels=FALSE)
    #axis(4,at=ymajor,labels=FALSE)
    #axis(4,at=yminor,labels=FALSE)
  }
  
  # Add the border box
  box(lwd=lwd.box,col=col.box)
  
}

## Adding histogram to a plot
my.hist <- function(hi,freq=TRUE,b=FALSE,col='grey95',border='grey70',lwd=1,lty=1,filled="standard")
{
  
  x  <- hi$mids
  x2 <- hi$breaks
  y  <- hi$counts
  if (!freq) y <- hi$density
  
  if ( b ) { # USE NEW METHOD DATA
    x  <- hi$mids2
    x2 <- hi$breaks2
    y  <- hi$counts2
    if (!freq) y <- hi$density2
  }

  ## Standard histogram plot
  # for ( i in 1:length(x) ) {
  #   dx <- (x2[i+1]-x2[i]) / 2
  #   rect(x[i]-dx,0,x[i]+dx,y[i],col=col,border=border,lwd=lwd,lty=lty)
  # }
  
  ## histogram outline without vertical lines between bins
   n <- length(x2)
   xfill <- c(x2[1],x2[1])
   yfill <- c(0,y[1])
   for ( i in 1:n ) {
     xfill <- c(xfill,x2[i],x2[i+1],x2[i+1],x2[i+1])
     yfill <- c(yfill,y[i],y[i],y[i],y[i+1])
   }
   xfill <- c(xfill,x2[n],x2[n])
   yfill <- c(yfill,y[n],0)
   polygon(xfill,yfill,col=col,border=NA)
   lines(xfill,yfill,col=border,lwd=lwd,lty=lty)
  
  ## Polygon around values
  # polygon(hi$density.call,col=col,border=border,lwd=lwd,lty=lty)
  #polygon(c(x[1],x,x[length(x)]),c(0,y,0),col=col,border=border,lwd=lwd,lty=lty)

}

get.conf <- function(x,weight,interval=95)
{ # Function to determine confidence/credence intervals
  # given a sample x and weights (weight should sum to 1)
  # Returns the indices that fall within the desired interval
  
  # Number of points
  n <- length(weight)

  # Order the weights from smallest to largest
  i1 <- order(weight)
  y1 <- weight[i1]

  # Loop over ordered weights, until interval is reached 
  # (ie, count until the interval of interest begins)
  w <- 0.0
  for ( i in 1:n) {
    w <- w + y1[i]
    if ( w >= (100-interval)*1e-2 ) break
  }
  
  # Collect the indices of points inside the desired interval
  # and the actual values
  ii95 <- i1[i:n]
  x95 <- x[ii95]

  # Return the indices
  return(ii95)
}

dnorm.integrated <- function(x,n=100,mean=0,sd=1)
{ # Get the integrated prob. density of bins (instead of just 1 value)
  # (function requested by Dr Rougier for statistics of hysteresis paper)

  x0 <- sort(unique(x))
  dd1 <- diff(x0)/2
  dd1 <- c(mean(dd1),dd1,mean(dd1))
  pp1 <- numeric(length(x0))
  p <- x*NA
  for ( q in 1:length(x0) ) {
    r0 <- x0[q]-dd1[q]
    r1 <- x0[q]+dd1[q+1]
    pp1[q] <- sum(dnorm(seq(r0,r1,length.out=n),mean=mean,sd=sd))/n
    p[which(x==x0[q])] <- pp1[q]
  }
  
  return(p)
}

samples.syn <- function(x,p,n=1e3,dx=0.01,spar=0.6,from=range(x,na.rm=T)[1],to=range(x,na.rm=T)[2]) 
{ ## Method: Get spline of cdf and sample a large number of
  ## points from it. Return individual points with probabilities attached and
  ## and evenly spaced pdf as derivative of smoothed spline
  
  # Determine the sampled cumulative density
  ii   <- order(x)
  x1   <- x[ii]
  dx1  <- diff(x1)
  p1   <- p[ii]
  cdf1 <- cumsum(p1)
  
  # Add additional points to beginning and end of x to make spline 0:1
  xbuffer <- cumsum(rep(mean(dx1,na.rm=T)*10,10))
  x1 <- c(x1[1]-rev(xbuffer),x1,x1[length(x1)]+xbuffer)
  cdf1 <- c( rep(0.0,length(xbuffer)), cdf1, rep(1.0,length(xbuffer)) )
  
  # Calculate the spline function of the cdf
  cdf1fun <- splinefun(x=x1,y=cdf1,method="monoH.FC")
  
  ### FIRST calculate spline at random x points to get individual probabilities
  min1 = range(x[p>0])[1]; max1 = range(x[p>0])[2]
  xr <- sort(runif(n,min=min1,max=max1))
  cdfr <- cdf1fun(xr,deriv=0)
  cdfr[cdfr>1] <- 1
  cdfr[cdfr<0] <- 0

  # Calculate the probability of each point p[i] = cdf[i]-cdf[i-1]
  pr <- cdfr[2:length(cdfr)] - cdfr[1:(length(cdfr)-1)]
  pr <- c(0,pr)
  pr <- pr/sum(pr)
  
  ## NEXT calculate spline at evenly-spaced intervals to get the pdf and cdf
  mids <- seq(from=from,to=to,by=dx)
  cdf <- cdf1fun(mids,deriv=0)
  cdf[cdf>1] <- 1
  cdf[cdf<0] <- 0
  
  ## SMOOTH THE CDF
  ## Note: running mean doesn't work, because
  ## the cdf flattens out too much

  ## Smooth the spline now using smooth.spline
  ## (This wasn't possible originally, bc smooth.spline 
  ##  doesn't have option for monotonic splines)
  cdf = smooth.spline(x=mids,y=cdf,spar=spar)$y
  
  # Calculate the probability density
  # (Derivative of cdf, using central-differencing)
  pdf <- cdf*0
  for ( i in 3:(length(pdf)-2) ) {
    
    # Central difference (4th order error)
    pdf[i] <- ( -cdf[i+2] + 8*cdf[i+1] - 8*cdf[i-1] + cdf[i-2] ) / (12*diff(mids)[1])
  }
  pdf[pdf<0] = 0
  
  return(list(mids=mids,density=pdf,cdf=cdf,x1=x1,cdf1=cdf1,xr=xr,pr=pr))    

}



my.density <- function(x,weights=rep(1,length(x)),from=NA,to=NA,n=512,dx=NA,dx2=0.1,type="cdf",calchist=TRUE)
{ # Function to calculate a PDF based on a sample of data x,
  # and given weights. Output should be similar to standard hist() function
  # Along with confidence intervals and a stats summary
  
  # First remove all na values
  ii =which(!is.na(x))
  x = x[ii]
  weights = weights[ii]

  lims0 <- range(x)

  # Check range of x values
  goodrange = TRUE
  if (diff(lims0)==0) goodrange = FALSE

  if (is.na(from)) from <- lims0[1]
  if (is.na(to))   to   <- lims0[2]

  # If dx is given (desired output binwidth), determine n
  if (!is.na(dx)) {
    tmp <- seq(from=from,to=to,by=dx)
    n   <- length(tmp)
  } else {
    tmp <- seq(from=from,to=to,length.out=n)
    dx  <- diff(tmp)[1]
  }
  
  # Make sure weights are normalized
  weights = weights / sum(weights)
  
  if (!goodrange) {
    density.call = NA
    mids = tmp
    density = mids*0+1
    density = density / sum(density) / dx 
  } else if ( type=="kde" ) { # Use kernel density estimate
    density.call <- density(x,weights=weights,from=from,to=to,n=n)
    mids    <- density.call$x
    density <- density.call$y
  } else { # Use smoothed empirical CDF
    density.call <- samples.syn(x=x,p=weights,from=from,to=to,dx=dx)
    mids    <- density.call$mids
    density <- density.call$density
  }

  dx <- diff(mids)[1]
  db <- dx/2
  breaks <- seq(from=mids[1]-db,to=max(mids)+db,length.out=length(mids)+1)
  counts  <- density*dx*length(x)
  
  ## NEW - calculate histogram instead of density
  if (calchist) {
    breaks2 <- seq(from=from,to=to,by=dx2)
    mids2   <- breaks2[1:(length(breaks)-1)]+dx2/2
    counts2 <- density2 <- density*0
    for ( i in 1:length(mids2) ) {
      ii <- which( x >= breaks2[i] & x < breaks2[i+1] )
      density2[i] <- sum(weights[ii])
      counts2[i]  <- length(ii)
    }
    density2 <- density2 / diff(breaks2)[1]
  
  } else {
    mids2  <- breaks2 <- counts2 <- density2 <- NA
  }
  ## END NEW
  
  # Get various confidence intervals of interest (from density output)
  x2 <- mids; weights2 <- density/sum(density); weights2[is.na(weights2)]<-0
  i95 <- get.conf(x2,weights2,interval=95)
  i90 <- get.conf(x2,weights2,interval=90)
  i66 <- get.conf(x2,weights2,interval=66)
  i10 <- get.conf(x2,weights2,interval=10)
  i01 <- get.conf(x2,weights2,interval=1)
  
  idens = get.conf(x2,weights2,interval=5)
  best  = mean(mids[idens],na.rm=T)
  # best  = mids[which.max(density)]

  # Generate an output table (from density output)
  table2 <- data.frame(best=best,
                      r66a=range(x2[i66])[1],r66b=range(x2[i66])[2],
                      r90a=range(x2[i90])[1],r90b=range(x2[i90])[2],
                      r95a=range(x2[i95])[1],r95b=range(x2[i95])[2])
                      #r100a=999,r100b=999,
                      #bestsim=x2[which.max(weights2)])
  
  # Get various confidence intervals of interest (from original data)
  i95 <- get.conf(x,weights,interval=95)
  i90 <- get.conf(x,weights,interval=90)
  i66 <- get.conf(x,weights,interval=66)
  i10 <- get.conf(x,weights,interval=10)
  i01 <- get.conf(x,weights,interval=1)
  
  # Generate an output table (from original data)
  # Determine best (from top few points)
  idens = get.conf(mids,density,interval=5)
  best  = mean(mids[idens],na.rm=T)
  # best  = mids[which.max(density)]

  table  <- data.frame(best=best,
                       r66a=range(x[i66])[1],r66b=range(x[i66])[2],
                       #r90a=range(x[i90])[1],r90b=range(x[i90])[2],
                       r95a=range(x[i95])[1],r95b=range(x[i95])[2])
                       #r100a=range(x)[1],r100b=range(x)[2],
                       #bestsim=x[which.max(weights)])
  
  # Store all results
  out <- list(mids=mids,counts=counts,density=density,breaks=breaks,density.call=density.call,
              mids2=mids2,breaks2=breaks2,counts2=counts2,density2=density2,
              x=x,y=weights,i95=i95,i90=i90,i66=i66,i10=i10,i01=i01,summary=table,summary2=table2)
  
  cat("Calculated my.density. Best x =",round(table2$best,3),'\n')

  return(out)
}

combine <- function(x,y,normalize=FALSE) 
{ # Function to combine all combinations of two vectors
  
  tmp <- expand.grid(x,y)
  out <- tmp[[1]]*tmp[[2]]
  
  if (normalize) out <- out/sum(out)
  
  cat("Combined values:",length(out),"\n")
  
  return(out)
}

plot.histweights <- function(dist)
{
  col <- rep(2,length(dist$x))
  col[dist$i95] <- "violet"
  col[dist$i66] <- 1
  col[dist$i10] <- 4
  plot(dist$x,dist$y,col=col,xlab="Temp (°C)",ylab="Weighting",xlim=c(0,6))

}
pointfit2 <- function(x0,y0,col=1,lwd=1,pch=1,lty=1)
{
  fit <- lm(y0~x0+I(x0^2))
  x1 <- sort(x0); y1 <- predict(fit,newdata=data.frame(x0=x1))
  points(x0,y0,pch=pch,col=col)
  lines(x1,y1,col=alpha(col,60),lwd=lwd,lty=lty)

}

# For weighted means, vars
require("Hmisc")
require("nortest")

### NON-LINEAR TRENDS ####
require(Rssa)
require(compiler)

series_fill <- function(y)
{ # Fill in the missing values of the time series 
  # with a linear interpolation
  
  # Generate a generic x vector
  # and a new y vector
  x = c(1:length(y))
  ynew = y

  # Find indices of missing values
  ii = which(is.na(y))
  
  # Approximate the missing values with a linear fit
  if (length(ii)>0) ynew[ii] = approx(x,y,xout=x[ii],rule=2)$y 

  return(ynew)
}

ssatrend <- function(y,L=(length(y)-1) %/% 2,fill=TRUE)
{ # Compute the trend from SSA decomposition
  
  # If fill requested, make sure that NAs are filled in
  if (fill) y = series_fill(y)

  # Get a new ssa object and
  # Reconstruct the first element
  ssa   = new.ssa(y,L=L)
  trend = reconstruct(ssa,groups=list(1))$F1
  
  # That's it! Return the trend...
  return(trend)
}

cmp.ssatrend = cmpfun(ssatrend)

### END NON-LINEAR TRENDS ####

## PLOTTING ##


gen.styles <- function(nm,n=length(nm),col=rep(1,n),lty=rep(1,n),
                       lwd=rep(1,n),pch=rep(23,n),pch.lwd=rep(1,n),cex=rep(1,n),bg=rep(NA,n))
{
  styles <- data.frame(nm=nm,col=col,lty=lty,lwd=lwd,pch=pch,pch.lwd=pch.lwd,bg=bg,cex=cex)
  styles$lty     <- as.numeric(styles$lty)
  styles$lwd     <- as.numeric(styles$lwd)
  styles$pch     <- as.numeric(styles$pch)
  styles$pch.lwd <- as.numeric(styles$pch.lwd)
  styles$bg      <- as.numeric(styles$bg)
  styles$cex     <- as.numeric(styles$cex)
  return(styles)
}

library(gridBase)

land.colors  <- colorRampPalette(c("tan","brown"),bias=5)
water.colors <- colorRampPalette(c("skyblue3","lightblue"))
grey.colors  <- colorRampPalette(c("grey50","grey80"))
abc.labels <- c("a","b","c","d","e","f","g")

# Default jet colors
jet.colors = c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F",
                 "yellow", "#FF7F00", "red", "#7F0000")

panel.contour.grl <- function(x,y,subscripts,...,datasets,
                              extra=extra,extra2=extra2,ptext=NA,plabel=NA,pobs=NULL,latlon=TRUE)
{
  # Determine which data set is being plotted
  pnow <- panel.number()
  
  # Get info for current data set
  land <- extra$land; mask <- extra$mask; elev <- extra$elev; points <- extra$points
  var  <- extra$var
  vnm  <- extra$vnms[ min(pnow,length(extra$vnms)) ]
  
  cat("var: ",vnm,"\n")
  
  # Get current data set
  pdat <- datasets[[ min(pnow,length(datasets)) ]]
  out  <- pdat[[vnm]]  
  ma   <- pdat$mask; zs  <- pdat$zs
  xx   <- pdat$xx; yy  <- pdat$yy
  ii0 <- c(1:length(out))
  
  cat("n=",length(out),"\n")  
  
  # Plot actual variable for subscripts desired 
  #panel.levelplot(x,y,subscripts=subscripts,...)
  ii1 <- which(ma %in% var$mask)
  panel.contourplot(xx, yy, out,subscripts=ii1,
                    at=var$at, col.regions=var$cols,
                    col=var$contour.col,lty=var$contour.lty,lwd=var$contour.lwd,
                    region = TRUE, contour = var$contour.plot, labels = FALSE)
  
  # Plot land if desired
  if (land$plot) {
    ii1 <- which(ma %in% land$mask)
    panel.contourplot(xx, yy, zs,subscripts=ii1,
                      at=land$at, col.regions=land$cols,
                      region = TRUE, contour = FALSE, labels = FALSE)
  }
  
  # Plot the ocean
  if ( ocean$plot ) {
    ii1 <- which(ma %in% ocean$mask)
    panel.contourplot(xx, yy, zs,subscripts=ii1,
                      at=ocean$at, col.regions=ocean$cols,
                      region = TRUE, contour = FALSE, labels = FALSE)
  }
  
  # Plot land/ice outline from mask
  if (mask$plot) {
    # First plot the desired mask contours
    panel.contourplot(xx, yy, ma, subscripts=ii0,lty=1,col=1,lwd=2,
                      at=mask$at,
                      region = FALSE, contour = TRUE, labels = FALSE)
  }
  
  # Add elevation contours (thicker lines for round numbers..)
  if (elev$plot) {
    ii1 <- which(ma %in% elev$mask)
    panel.contourplot(xx, yy, zs,subscripts=ii1,lty=2,col="grey20",lwd=0.5,
                      at=elev$at,region = FALSE, contour = TRUE, labels = FALSE)
  }
  
  # Add additional contour field (eg, data for comparison
  if (!is.null(extra2)) {
    # Plot the desired mask contours
    panel.contourplot(xx, yy, extra2$mask, subscripts=ii0,lty=extra2$lty,col=extra2$col,lwd=extra2$lwd,
                      at=extra2$at,
                      region = FALSE, contour = TRUE, labels = FALSE)
  }
  
  # Add observational or other points of interest to plot
  if (points$plot) {
    #c2 <- rgb(t(1-col2rgb(points$col)))
    panel.points(points$x,points$y,col=points$col,pch=points$pch,fill=points$bg,alpha=0.6,lwd=points$lwd,cex=points$cex)
    ltext(points$x,points$y,points$lab,cex=points$cex.lab,col=points$col.lab) #,adj=c(0.5,-1))
  }
  
  ## Add latitude/longitude lines
  if (latlon) {
  
    latticks <- c(60,70,80); lonticks <- c(-75,-60,-45,-30,-15)
    panel.contourplot(xx, yy, pdat$lat,subscripts=ii0,lty=2,col="grey40",lwd=0.5,
                      at=latticks,
                      region = FALSE, contour = TRUE, labels = FALSE)
    panel.contourplot(xx, yy, pdat$lon,subscripts=ii0,lty=2,col="grey40",lwd=0.5,
                      at=lonticks,
                      region = FALSE, contour = TRUE, labels = FALSE)
    
    ## Latitude ticks
    xs <- numeric(length(latticks))+min(pdat$xx)
    ys <- c(-3305,-2120,-835)
    yslab <- c("60°N","70°","80°")
    ysrt  <- c(-12,-15,-42)
    for (i in 1:length(ys)) ltext(xs[i],ys[i],yslab[i],srt=ysrt[i],adj=c(-0.1,0),col=1,cex=1.1)
    
    ## Longitude ticks
    xs <- c(-455,-240,-55,110,320)
    ys <- numeric(length(lonticks))+max(pdat$yy)-60
    yslab <- c("-75°","-60°","-45°","-30°","-15°E")
    for (i in 1:length(ys)) ltext(xs[i],ys[i],yslab[i],adj=c(0.5,0),col=1,cex=1.1)
  
  }
  
  # Panel label
  if (!is.na(plabel[1])) {
    x1 <- min(pdat$xx)+130; y1 <- max(pdat$yy)-100
    ltext(x1,y1,plabel[pnow],cex=2.5,col="grey40") #,vfont=c("serif","bold"))
  }
  
  # And title
  if (!is.na(ptext[1])) ltext(300,-2750,ptext[ min(pnow,length(ptext)) ],cex=1.5,pos=1)

}

get.col = function(x,col=c("blue","red"),n=20,mid=NA,extend=10,ii=c(1:length(x)))
{
  nx = length(x)
  xlim = range(x,na.rm=T)
  
  if (!is.na(mid)) {
    # Caculate the xlim, so that mid is the midpoint! (eg, mid=0)
  }
  
  dxlim = diff(xlim)
  xlim[1] = xlim[1] - (extend/100)*dxlim
  xlim[2] = xlim[2] + (extend/100)*dxlim

  breaks = pretty(xlim,n)
  db = diff(breaks)[1]
  nb = length(breaks)

  palette = colorRampPalette(col)(nb-1)
  cols = rep(NA,nx)
  for ( i in 1:nx ) {
    j = which.min( abs(x[i] - breaks+db/2) )
    if ( length(j)==1 ) cols[i] = palette[j]
  }
  
  jj = which(! c(1:length(x) %in% ii))
  cols[jj] = NA
  return(list(breaks=breaks,palette=palette,col=cols))
}

get.contours <- function(nm,var=NULL,zrange=NULL,n=15,alpha=100,darker=0,
                col=NULL,bias=NULL,at=NULL,at.mark=NULL,rev=FALSE)
{
  
  bias1 <- 1
  
  if (!is.null(var)) {
    zrange1 <- range(abs(var),na.rm=TRUE); zrange1 <- c(0,zrange1[2])
  }
  
  if ( nm == "tt" ) {
    
    at1       <- c(-40,-30,-20,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10)
    at.mark1  <- as.character(at1)
    bias1     <- sum(at1>0) / sum(at1<0)
    col1   <- c("magenta","blue","cyan","white","yellow","red","darkred") 
  
  } else if ( nm %in% c("pp","snow") ) {
    
    if (is.null(var)) zrange1 <- c(0,2500)

    at1      <- pretty(zrange1,n)
    at.mark1 <- as.character(at1)
    col1     <- c("magenta","blue","cyan","green","yellow","red","darkred")
    #col1     <- c("white","yellow","green","darkgreen","cyan","blue","magenta","darkviolet")  
    #col1     <- c("white","yellow","green","cyan","blue","magenta")
    
    cat(" -used pp/snow contours.\n")
    
  } else if ( nm == "smb" ) {
    
    at1      <- c(-6000,-4000,-2000,-1000,-800,-600,-400,-200,-100,0,100,200,400,600,800,1000,2000,4000)
    at.mark1 <- as.character(at1)
    bias1    <- 8/9
    col1     <- c("magenta","blue","cyan","white","yellow","red","darkred")
    #col1     <- c("darkred","red","yellow","white","cyan","blue","magenta")
    
  } else if ( nm %in% c("melt","runoff") ) {
    
    at1      <- c(0,100,200,400,600,800,1000,2000)
    at.mark1 <- as.character(at1)
    col1     <- c("magenta","blue","cyan","white","yellow","red","darkred")
    #col1     <- c("darkred","red","yellow","white","cyan","blue","magenta")
    
  } else if ( nm == "dsmb" ) {
    
    at1      <- c(-6000,seq(-1000,1000,by=100),6000)
    at.mark1 <- as.character(at1)
    at.mark1[1] <- at.mark1[length(at.mark1)] <- ""
    bias1    <- 8/9
    col1     <- c("magenta","blue","cyan","white","yellow","red","darkred")
    #col1     <- c("darkred","red","yellow","white","cyan","blue","magenta")
    
  } else if ( nm %in% c("zs") ) {
    
    at1      <- seq(from=0,to=3300,by=300)
    at.mark1 <- as.character(at1)
    col1     <- c("grey80","white")
    #col1     <- c("darkslategray","white")
    
  } else if ( nm %in% c("land") ) {
    
    at1      <- seq(from=-100,to=3300,by=100)
    at.mark1 <- as.character(at1)
    #at.mark1[which(!at.mark1 %in% c("0","500","1000","1500","2000","2500","3000","3500"))] <- " "
    bias1    <- 5
    col1     <- c("tan","brown")
  
  } else if ( nm %in% c("ocean") ) {
    
    at1      <- seq(from=-4500,to=100,by=100)
    at.mark1 <- as.character(at1)
    #at.mark1[which(!at.mark1 %in% c("0","500","1000","1500","2000","2500","3000","3500"))] <- " "
    col1     <- c("paleturquoise","paleturquoise")
    #col1     <- c("white","white")
  
  } else if ( nm %in% c("dttp1") ) {  # Generic (positive or negative only)

    zrange1  <- zrange
    at1      <- pretty(zrange1,n)
    at.mark1 <- as.character(at1)
    col1     <- c("magenta","cyan","white","yellow","orange","red","darkred")
    
    cat("dttp colors for ",nm,":",zrange1,"\n")
    
  } else if ( nm %in% c("pos","neg") ) {  # Generic (positive or negative only)

    zrange1  <- zrange
    at1      <- pretty(zrange1,n)
    at.mark1 <- as.character(at1)
    col1     <- c("magenta","blue","cyan","yellow","red","darkred")
    
    cat("Generic colors for ",nm,":",zrange1,"\n")
    
  } else {  # Generic (negative and positive)
    
    #zrange1  <- max(abs(zrange))
    #zrange1  <- c(-zrange1,zrange1)
    zrange1  <- zrange
    at1      <- pretty(zrange1,n)
    at.mark1 <- as.character(at1)
    bias1    <- sum(at1>0) / max(1,sum(at1<0))
    col1     <- c("magenta","blue","cyan","lightgreen","white","yellow","orange","red","darkred")
    
    cat("Generic colors for ",nm,":",zrange1,"\n")
  }
  
  # Assign newly generated values if needed
  if (is.null(col))       col     <- col1
  if (is.null(bias))      bias    <- bias1
  if (is.null(at))        at      <- at1
  if (is.null(at.mark))   at.mark <- at.mark1
  
  if (rev) col <- rev(col)
  
  # Now generate actual colors via colorwheel
  # Shift colors darker or lighter, and add alpha layer (transparency)
  colorwheel <- colorRampPalette(col,bias=bias)
  cols <- colorwheel(length(at)-1)
  cols <- darker(cols,darker) 
  cols <- alpha(cols,alpha)
  
  return(list(at=at,at.mark=at.mark,cols=cols,colorwheel=colorwheel,
              contour.plot=FALSE,contour.lwd=0.5,contour.lty=1,contour.col="white"))
}


fill_mask <- function(mask)
{ # Fill an ice sheet mask so there are no isolated grid types
  
  nx <- dim(mask)[1]
  ny <- dim(mask)[2]
  
  fill <- 0
  dn   <- 1
  
  for ( i in 3:(nx-2) ) {
    for ( j in 3:(ny-2) ) {
      
      # For now only fill land points that should be ice
      if ( mask[i,j] == 1 ) {
        
        neighbs <- mask[(i-dn):(i+dn),(j-dn):(j+dn)]
        tot <- base::sum(neighbs)-1
        
        # If all neighbors == 0, then it is isolated by ice
        if ( tot < 4 ) {
          mask[i,j] <- 0
          fill <- fill+1
        }

      }
    }
  }
  
  cat("Filled mask",fill,"times.\n")
  
  return(mask)
}

## Add latitude and longitude lines and labels
add_latlon <- function(grid,x=grid$xx[,1],y=grid$yy[1,],col="grey30",lty=2,lwd=0.5,cex=0.9,shift.lat=0,shift.lon=0)
{
  lats <- grid$lat
  lons <- grid$lon
  
  # Make sure lat/lon is range [-180:180]
  if (max(lons) > 180 ) {
    ii = which(lons>180)
    lons[ii] = lons[ii] - 360
  } 

  # Lat/lons
  latticks <- c(60,70,80)
  lonticks <- c(-75,-60,-45,-30,-15)
  contour(x=x,y=y,z=lats,add=TRUE,levels=latticks,drawlabels=FALSE,lty=lty,lwd=lwd,col=darker(1,30))
  contour(x=x,y=y,z=lons,add=TRUE,levels=lonticks,drawlabels=FALSE,lty=lty,lwd=lwd,col=darker(1,30))

  ## Latitude ticks
  xs <- numeric(length(latticks))+min(x)
  ys <- c(-3315,-2125,-845) - 10 + shift.lat
  yslab <- c("60°N","70°","80°")
  ysrt  <- c(-12,-15,-38)
  for (i in 1:length(ys)) text(xs[i],ys[i],yslab[i],srt=ysrt[i],adj=c(-0.1,0),col=col,cex=cex)
  
  ## Longitude ticks
  xs <- c(-455,-240,-55,110,320)
  ys <- numeric(length(lonticks))+max(y)-72 + shift.lon
  yslab <- c("-75°","-60°","-45°","-30°","-15°E")
  for (i in 1:length(ys)) text(xs[i],ys[i],yslab[i],adj=c(0.5,0),col=col,cex=cex)

}

  
# My own image/contour plotting function
my.image <- function(dat,nm="zs",ii=c(1:length(dat[[nm]])),col=NULL,breaks=NULL,mask=TRUE)
{
  
  x <- dat$xx[,1]; y <- dat$yy[1,]
  
  ## Plot the elevation: filled contour
  var <- dat[[nm]]*NA
  var[ii] <- dat[[nm]][ii]
  
  image(x=x,y=y,z=var,axes=FALSE,col=col,breaks=breaks,xlab="",ylab="")
  
  ## Add land/ice mask
  if (mask) contour(x=x,y=y,z=dat$mask,add=TRUE,drawlabels=FALSE,nlevels=2,lwd=2,col=alpha(1,80))
  contour(x=x,y=y,z=dat$zs,add=TRUE,drawlabels=FALSE,nlevels=10,lwd=0.5,lty=1,col=alpha(1,50))
  
  ## Add lat/lon
  contour(x=x,y=y,z=dat$lat,add=TRUE,drawlabels=FALSE,nlevels=4,lwd=0.5,lty=2)
  contour(x=x,y=y,z=dat$lon,add=TRUE,drawlabels=FALSE,nlevels=4,lwd=0.5,lty=2)
  
  box()
}

my.image3 <- function(dat,nm="zs",suffix=NA,sim=1,k=1,ii=c(1:length(dat[[nm]][sim,,,k])),col=NULL,breaks=NULL,
                      plt.var=TRUE,plt.var.cont=FALSE,plt.ocean=TRUE,plt.land=TRUE,plt.mask=TRUE,plt.latlon=TRUE,lwd=0.5,
                      col.land=NULL,col.latlon=alpha(1,50),cex.latlon=0.9,lwd.mask=2,fill=FALSE)
{
  
  vnm      <- nm
  vnm.mask <- "mask"
  vnm.zs   <- "zs"
  if (!is.na(suffix)) {
    vnm      <- nm
    vnm.mask <- "mask"
    vnm.zs   <- "zs"
  }

  ## Plot the variable: filled contour
  var0 <- dat[[vnm]][sim,,,k]
  var <- var0*NA
  var[ii] <- var0[ii]

  # Get x and y vectors for plotting
  dims = dim(dat[[vnm]][sim,,,k])
  x = seq(0,1,length.out=dims[1])
  y = seq(0,1,length.out=dims[2])
  if ( "xx" %in% names(dat) ) x = dat$xx[,1]
  if ( "yy" %in% names(dat) ) y = dat$yy[1,]
  
  ## Limit variable value to range of contours
  if (!is.null(breaks)) {
    brange <- range(breaks)
    var[var>brange[2]] <- brange[2]
    var[var<brange[1]] <- brange[1]
  }
  
  # Get generic variables
  mask  <- dat[[vnm.mask]][sim,,,k]
  if (fill) mask <- fill_mask(mask)
  
  land  <- dat[[vnm.zs]][sim,,,k]; land[mask != 1]  <- NA
  ocean <- dat[[vnm.zs]][sim,,,k]; ocean[mask != 2] <- NA
  
  zs    = dat[[vnm.zs]][sim,,,k]

#   if (is.null(breaks)) {
#     conts <- get.contours(nm="pos",zrange=range(var,na.rm=TRUE),n=10)    
#     breaks <- conts$at
#     col    <- conts$cols
#   }
  
  # Start the empty plot
  par(xaxs="i",yaxs="i")
  plot(range(x),range(y),type="n",axes=FALSE,xlab="",ylab="")
  
  # Add ocean
  if ( plt.ocean ) {
    cont.ocean <- get.contours(nm="ocean",darker=20)
    ocean[ocean<min(cont.ocean$at)] = min(cont.ocean$at)
    ocean[ocean>max(cont.ocean$at)] = max(cont.ocean$at)
    image(x=x,y=y,z=ocean,add=TRUE,col=cont.ocean$cols,breaks=cont.ocean$at)
  }
  
  # Add land
  if ( plt.land ) {
    cont.land <- get.contours(nm="land",col=col.land)
    land[land<min(cont.land$at)] = min(cont.land$at)
    land[land>max(cont.land$at)] = max(cont.land$at)
    image(x=x,y=y,z=land,add=TRUE,col=cont.land$cols,breaks=cont.land$at)
  }
  
  # Get contours for ice elevation (for consistency
  cont.zs <- get.contours(nm="zs")
  cont.zs$lwd <- rep(lwd,length(cont.zs$at))
  cont.zs$lwd[cont.zs$at %in% c(2400)] <- lwd*2
  
  # Add variable
  if (plt.var) {
    var[var<min(breaks)] = min(breaks)
    var[var>max(breaks)] = max(breaks)
    image(x=x,y=y,z=var,add=TRUE,col=col,breaks=breaks)
    if (plt.var.cont) contour(x=x,y=y,z=var,add=TRUE,col="grey70",lwd=0.8,levels=breaks,drawlabels=FALSE)
  }
  
  ## Add land/ice mask
  if (plt.mask) contour(x=x,y=y,z=mask,add=TRUE,drawlabels=FALSE,nlevels=2,
                        lwd=lwd.mask,col=alpha(1,70))
  
  ## Add land contours
  if (plt.var) {
    contour(x=x,y=y,z=zs,add=TRUE,drawlabels=FALSE,levels=cont.zs$at,
            lwd=cont.zs$lwd,lty=1,col=alpha(1,50))
  }
  
  if (plt.latlon) {
    ## Add lat/lon
    add_latlon(dat,x=x,y=y,col=col.latlon,cex=cex.latlon,shift.lat=30,shift.lon=10)
  }
  
  
  box()
  
  return(list(col=col,breaks=breaks))
}

my.leg2 <- function(breaks,col,labs=paste(breaks),units="mm",plt=c(0,1,0,1),
                    xlab="",ylab="",xlim=NULL,ylim=NULL,zlim=range(breaks),
                    cex=1,cex.lab=1,new=FALSE,vertical=TRUE,line=1.8,
                    asp=1,mgp=c(3,0.5,0),at=breaks,labels=NULL,extend=F,...)
{
  
  kk = which(breaks >= zlim[1] & breaks <= zlim[2])

  n <- length(breaks[kk])
  x00 <- 0; y00 <- 0; dr <- 1/(n-1)
  
  # Generate points ranging from 0:1 for rectangles
  if ( vertical ) {
    yy <- seq(from=0,to=1,length.out=n)
    y0 <- yy[1:n]; y1 <- yy[1:n]+dr
    x0 <- rep(0,n); x1 <- rep(1,n)
  } else {
    xx <- seq(from=0,to=1,length.out=n)
    x0 <- xx[1:n]; x1 <- xx[1:n]+dr
    y0 <- rep(0,n); y1 <- rep(1,n)
  }

  dx = dy = 1
  drx <- dx/(n-1)
  dry <- dy/(n-1)
  
  #rect(x0-drx,y0-dry,x1+drx,y1+dry,col="grey95",border=NA)
  
  xlim0 = range(x0,x1)
  ylim0 = range(y0[1:(n-1)],y1[1:(n-1)])
  if (!vertical) xlim0 = range(x0[1:(n-1)],x1[1:(n-1)])
  if (!vertical) ylim0 = range(y0,y1)

  par(new=new,plt=plt,xpd=NA,xaxs="i",yaxs="i",cex=cex,...)
  plot( xlim0,ylim0, type="n",axes=F,ann=F)

  if ( vertical ) {
    
    for ( k in 1:n ) {
      # if (k < n) rect(x0[k],y0[k],x1[k],y1[k],col=col[k],border="grey60",lwd=0.4)
      if (k < n) rect(x0[k],y0[k],x1[k],y1[k],col=col[kk[k]],border=NA,lwd=0.4)
      #text(x1[k],y0[k],labs[k],pos=4,cex=cex,col="grey10") 
    }
    text(x0[n],y1[n-1]+dy*0.05,units,pos=4,cex=cex*1.2,col="grey10")
    
    par(new=T,xaxs="i",yaxs="i",...)
    plot(c(0:1),zlim,type="n",axes=F,ann=F,xlim=xlim,ylim=ylim)
    axis(4,at=at,labels=labels,mgp=c(3,0.5,0),tcl=0.0,col="grey10",cex.axis=cex)
    box(col="grey10")

  } else {
    
    for ( k in 1:n ) {
      if (k < n) rect(x0[k],y0[k],x1[k],y1[k],col=col[k],border=NA,lwd=0.4)
      #text(x0[k],y0[k],labs[k],pos=1,cex=cex,col="grey10") 
    }
    #text(x1[n-1]+dx*0.00,0.4,units,pos=4,cex=cex*1.2,col="grey20")
    mtext(side=4,line=1,las=1,units,col="grey20",cex=cex)
    
    # Add triangles to ends
    if (extend) {
      pp = list(x=c(x0[n],x0[n]+drx*0.2,x0[n]),y=c(y0[n],y0[n]+0.5,y1[n]))
      polygon(pp,col=col[n-1],border="grey10")
      pp = list(x=c(x0[1],x0[1]-drx*0.2,x0[1]),y=c(y0[1],y0[1]+0.5,y1[1]))
      polygon(pp,col=col[1],border="grey10")
    }
    
    # Adjust axis distance if horizontal
    mgp[2] = mgp[2]*0.2

    par(new=T,xaxs="i",yaxs="i",...)
    plot(zlim,c(0:1),type="n",axes=F,ann=F,xlim=xlim,ylim=ylim)
    axis(1,at=at,labels=labels,mgp=mgp,tcl=0.0,col="grey10",cex.axis=cex*1.1)
    box(col="grey10")

  }
  
  mtext(side=1,line=line,xlab,cex=cex.lab)
  mtext(side=2,line=line,ylab,cex=cex.lab)

  par(xpd=FALSE)
}

my.leg <- function(breaks,col,units="mm",x00=-300,y00=-3500,dr=70)
{
  n <- length(col)

  x11 <- x00 + (n-1)*dr/2+dr
  y11 <- y00+2*dr
  
  par(xpd=NA)
  
  rect(x00-2*dr,y00-2.5*dr,x11+2*dr,y11+0.5*dr,col="grey95",border=1)
  
  for ( k in 1:n ) {
    x0 <- x00+(k-1)*dr/2; x1 <- x0+dr
    y0 <- y00; y1 <- y0+2*dr
    rect(x0,y0,x1,y1,col=col[k],border=NA)
    if (k==1             ) text(x0,y0,paste(breaks[k]),pos=1,cex=0.9,col="grey10")
    if (k==length(breaks)-1) text(x1,y0,paste(breaks[k+1],units),pos=1,cex=0.9,col="grey10")
  }
  
  par(xpd=FALSE)
}

par.uin <- function()
  # determine scale of inches/userunits in x and y
  # from http://tolstoy.newcastle.edu.au/R/help/01c/2714.html
  # Brian Ripley Tue 20 Nov 2001 - 20:13:52 EST
{
    u <- par("usr")
    p <- par("pin")
    c(p[1]/(u[2] - u[1]), p[2]/(u[4] - u[3]))
}

quiver<- function(xx=NULL,yy=NULL,uu,vv,scale=0.13,length=0.008,add=TRUE,col=1,lwd=1,lty=1,thin=1)
# first stab at matlab's quiver in R
# from http://tolstoy.newcastle.edu.au/R/help/01c/2711.html
# Robin Hankin Tue 20 Nov 2001 - 13:10:28 EST
{
    if (thin > 1) {
      # Thin arrays for vector plotting 
      ii = seq(2,dim(marc$xx)[1],by=as.integer(thin))
      jj = seq(2,dim(marc$xx)[2],by=as.integer(thin))
      xx = xx[ii,jj]
      yy = yy[ii,jj]
      uu = uu[ii,jj]
      vv = vv[ii,jj]
    }

    if (is.null(xx)) {
      xx <- col(uu)
      xx = xx/max(xx)
    }
    if (is.null(yy)) {
      yy <- max(row(uu))-row(uu)
      yy = yy/max(yy)
    }    

    speed <- sqrt(uu*uu+vv*vv)
    maxspeed <- max(speed,na.rm=TRUE)

    ux <- uu*scale/maxspeed
    vy <- vv*scale/maxspeed

    #matplot(xpos,ypos,add=add,type="p",cex=0)
    arrows(xx,yy,xx+ux,yy+vy,length=length*min(par.uin()),col=col,lwd=lwd,lty=lty)
}

plot.series <- function(x,y,xlim=NULL,ylim=NULL,col=NULL,good=NULL,
                        axis=c(1,2),lwd=2,lty=1,alpha=20,title=NULL)
{ # x: x variable (vector)
  # y: y variables (array, n X nx)

  if (is.null(xlim)) xlim = range(x,na.rm=T)
  if (is.null(ylim)) {
    kk = which(x >= xlim[1] & x <= xlim[2])
    ylim = range(y[,kk],na.rm=T)

  }
  
  plot(xlim,ylim,type="n",axes=F,ann=F)
  grid()
  axis(1)
  axis(2)
  
  sims = c(1:dim(y)[1])
  if (is.null(good)) good = sims
  
  # First plot bad sims with transparency
  qq = which(! sims %in% good)
  for (q in qq) lines(x,y[q,],lwd=lwd,lty=lty,col=alpha(col[q],alpha))
  
  # Then plot good sims
  qq = which(sims %in% good)
  for (q in qq) lines(x,y[q,],lwd=lwd,lty=lty,col=col[q])
  
  # Plot label
  if (!is.null(title)) text(xlim[1],ylim[2]-diff(ylim)*0.05,pos=4,title,cex=1)
  
  box()
}

## 
## PLOT FUNCTION TEMP...
plot.paleo <- function(runs,nm="Vtot",xnm="time",ii=NA,xlab="",ylab="",title="",col=NULL,bad=NA,axis=NA,lwd=2,lty=1,
                       hline=NA,h.lwd=2,h.lty=1,h.col="grey60",lwd.box=1,
                       xmaj=6,ymaj=4,xmin=xmaj*2-1,ymin=ymaj*2-1,
                       x.at=NULL,y.at=NULL,
                       convert=1,xlim=range(runs[[xnm]],na.rm=TRUE),cex=1,at=NULL,
                       ylim=range(runs[[nm]],na.rm=TRUE)*convert,label="",add=FALSE,interp=NA )                       
{
  
  # Determine number of runs to plot (and indices of runs)
  nruns <- dim(runs[[nm]])[2]
  
  # Generate x array (in case time is used, make array for number of runs)
  x00 <- runs[[xnm]]
  if (length(x00) != length(runs[[nm]])) x00 <- array(runs[[xnm]],dim=dim(runs[[nm]]))
  
  # Generate a color palette based on the number of runs
  if (is.null(col)) col <- colorRampPalette(c("red","yellow","blue","black","green","magenta","cyan"))(nruns)
  if (length(col)==1) col <- rep(col,nruns)
  if (is.na(bad[1]))  bad <- numeric(nruns)
  if (length(lwd)==1) lwd <- rep(lwd,nruns)
  if (length(lty)==1) lty <- rep(lty,nruns)
  if (length(lty)==2) { lty2 <- lty[2]; lty <- rep(lty[1],nruns); lty[which(bad!=0)] <- lty2 }
  
  # Determine indices to plot (if not given)
  if (is.na(ii[1])) ii <- c(1:nruns)
  
  # Extend ylim a bit to allow for "title"
  yr <- diff(ylim)*0.15
  t0 <- title; if ( is.expression(title) ) t0 <- "title exists"
  if ( t0 != "" ) ylim[2] <- ylim[2] + yr
  xr <- diff(xlim)*0.02
  
  # Determine the axes
  xx0 <- xlim[1]
  if ( xlim[1] %% 1 != 0.0 ) {
    offset <- xlim[1] %% 1
    if ( (xlim[1]+offset) %% 1 != 0.0 ) offset <- -offset
    xx0 <- xlim[1] + offset
  }
  # Get major and minor axes
  xmajor <- pretty(xlim,n=xmaj); xminor <- pretty(xlim,n=xmin)
  ymajor <- pretty(ylim,n=ymaj); yminor <- pretty(ylim,n=ymin)
  
  if (!is.null(x.at)) { 
    xmajor <- x.at
    xminor <- seq(from=min(x.at),to=max(x.at),length.out=1+2*(length(x.at)-1))
  }
  
  if (!is.null(y.at)) { 
    ymajor <- y.at
    yminor <- seq(from=min(y.at),to=max(y.at),length.out=1+2*(length(y.at)-1))
  }
  
  # Generate width of masking rectangle
  tmp <- diff(xlim); xp <- mean(xlim) + c(-tmp,tmp)*2
  tmp <- diff(ylim); yp <- mean(ylim) + c(-tmp,tmp)*2
  
  ## Now make the plots ##
  if (!add) {
    plot(range(x00[,ii],na.rm=TRUE),range(runs[[nm]][,ii],na.rm=TRUE)*convert,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE)
    #grid(col=8)#,lty=1)
    abline(v=xminor,h=yminor,col="lightgrey",lty=3)
    #cat("xmajor:",xmajor,"\n")
    #cat("xminor:",xminor,"\n")
    if ( !is.na(hline[1]) ) abline(h=hline,lwd=h.lwd,lty=h.lty,col=h.col)
  }
  
  # Loop over runs, but plot grey (bad) runs first
  for ( i in ii ) {
    if ( bad[i] > 0 ) {
      x <- x00[,i]; y <- runs[[nm]][,i]*convert
      if ( !is.na(interp) ) { 
        x1 <- seq(from=min(x),to=max(x),length.out=length(x)*interp)
        tmp <- approx(x,y,x1)
        x <- tmp$x; y <- tmp$y
      }
      lines(x,y,lwd=lwd[i],col=col[i],lty=lty[i])
    }
  }
  
  #   # Now plot semi-transparent rectangle to dim the "bad" runs
#   cshade <- col2rgb("white"); cshade <- rgb(t(cshade),alpha=200,max=255)
#   rect(xp[1],yp[1],xp[2],yp[2],col=cshade,border=NA)

  for ( i in ii ) {
    if ( bad[i] == 0 ) {
      x <- x00[,i]; y <- runs[[nm]][,i]*convert
      if ( !is.na(interp) ) { 
        x1 <- seq(from=min(x),to=max(x),length.out=length(x)*interp)
        tmp <- approx(x,y,x1)
        x <- tmp$x; y <- tmp$y
      }
      lines(x,y,lwd=lwd[i],col=col[i],lty=lty[i])
    }
  }
  
  ## Plot the meta text
  if (!add) {
    text(xlim[2]-xr,ylim[1],pos=2,xlab,cex=cex)
    
    xtxt <- xlim[1]+xr; adj <- c(0,1)   #; if ( 4 %in% axis ) { xtxt <- xlim[2]-xr; adj <- c(1,1) }
    text(xtxt,ylim[2]-yr*0.1,title,adj=adj,cex=cex)
    
    # Label (a, b, c...)
    if (!label=="") {
      xtxt <- xlim[2]-xr*1.5; adj <- c(1,1)
      text(xtxt,ylim[2],label,adj=adj,cex=cex,col="grey20")#,vfont=c("serif","bold"))
    }
    
    ## Now the axes and box
    for (i in 1:length(axis)) {
      if (axis[i] %in% c(1,3)) axis(axis[i],at=x.at)
      if (axis[i] %in% c(2,4)) axis(axis[i],at=y.at)
    }
    
#     if ( !is.na(axis[1]) ) axis(axis[1],at=x.at)
#     if ( !is.na(axis[2]) ) axis(axis[2],at=y.at)
#     if ( !is.na(axis[3]) ) axis(axis[3],at=x.at)
#     if ( !is.na(axis[4]) ) axis(axis[4],at=y.at)
  }
  
  box(lwd=lwd.box)
}

plot.panel <- function(x,y,ii=c(1:length(x)),xlab="",ylab="",title="",col=1,col.xlab=1,bad=NA,axis=NA,pch=21,lwd=2,
                       hline=NA,h.lwd=2,h.lty=1,h.col="grey60",label="",grid.h=NA,grid.v=NA,
                       convert=1,xlim=range(x,na.rm=TRUE),cex=1,cex.pch=1.8,at=NULL,
                       ylim=range(y,hline,na.rm=TRUE)*convert )                       
{
  
  nruns <- length(y)
  
  # Determine number of runs to plot (and indices of runs)
  if (is.na(bad[1]))  bad <- rep(0,nruns)
  if (length(pch)==1) pch <- rep(pch,nruns)
  if (length(col)==1) col <- rep(col,nruns)
  
  # Extend ylim a bit to allow for "title"
  yr <- diff(ylim)*0.1
  ylim[2] <- ylim[2] + yr
  xr <- diff(xlim)*0.02

  ## Now make the plots ##
  plot(x,y*convert,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE)
  if (is.na(grid.h[1]) & is.na(grid.v[1])) grid(col="lightgrey")
  abline(h=grid.h,v=grid.v,lty=3,col="lightgrey")
  
  if ( !is.na(hline[1]) ) abline(h=hline,lwd=h.lwd,lty=h.lty,col=h.col)
  
  # Loop over runs, but plot grey (bad) runs first
  ii1 <- ii[which(bad[ii] > 0)]
  points(x[ii1],y[ii1]*convert,pch=pch[ii1],lwd=lwd,col=col[ii1],cex=cex.pch)
  
  # Now plot good points
  ii1 <- ii[which(bad[ii] == 0)]
  points(x[ii1],y[ii1]*convert,pch=pch[ii1],lwd=lwd,col=col[ii1],cex=cex.pch)
  
  ## Now the axes
  if ( !is.na(axis[1]) ) { 
    axis(axis[1],at=at)
    for ( i in 2:length(axis)) axis(axis[i]) 
  }
  
  #text(xlim[2]-xr,ylim[1],pos=2,xlab,cex=cex,col=col.xlab)
  text(xlim[2]-xr,ylim[1],adj=c(1,1),xlab,cex=cex,col=col.xlab)
  
  xtxt <- xlim[1]; adj <- c(0,1)   #; if ( 4 %in% axis ) { xtxt <- xlim[2]-xr; adj <- c(1,1) }
  text(xtxt,ylim[2],adj=adj,title,cex=cex)
  
  # Label (a, b, c...)
  if (!label=="") {
    xtxt <- xlim[2]-xr*1.5; adj <- c(1,1)
    text(xtxt,ylim[2],label,adj=adj,cex=cex,col="grey20")
  }
  
  box(lwd=1)

}

## Function to plot panels corresponding to parameter choices
plot.panels <- function(info,pp=c("itm_c","slide0","pkappa"),fpp="Vtot",ylim=c(2,3.5),
                        title="",pch=20,col=1,ann=FALSE)
{
  
  for ( i in 1:length(pp) ) {
    
    p <- pp[i]
    
    xlim <- range(info[[p]]); r <- diff(xlim)*0.1; xlim[1] <- xlim[1]-r; xlim[2] <- xlim[2]+r
    
    plot(info[[p]],info[[fpp]],type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=ylim)
    grid()
    
    points(info[[p]],info[[fpp]],pch=pch,cex=2,col=col)
    
    axis(1,at=unique(info[[p]])); axis(2); box()
    
    if ( ann ) title(p) 
    if ( i == 1 ) mtext(title,side=2,line=2.5,cex=0.9,las=0)
    
  }

}

## Function to plot borehole temperature profiles
plot.core <- function(core,cnm="grip",nm="temp",xlab="Temperature (°C)",ylab="Depth (km)",
                           col=1,bad=NA,ii=NA,lwd=1,convert.y=1e-3,
                           ylim=NULL,xlim=NULL,title=""  )
{
    
  # Define the names (to include the core name)
  ynm <- paste("grip","zs",sep=".")
  xnm <- paste(cnm,nm,sep=".")
  
  # Determine number of runs to plot (and indices of runs)
  nruns <- dim(core[[xnm]])[2]
  if (is.na(ii[1])) ii <- c(1:nruns)
  if (is.na(bad[1])) bad <- numeric(nruns)
  if (length(col)==1) col <- rep(col,nruns)
  if (length(lwd)==1) lwd <- rep(lwd,nruns)
  
  depth <- core[[ynm]]*convert.y
  if (is.null(ylim)) ylim <- range(depth,na.rm=TRUE)
  if (is.null(xlim)) xlim <- range(core[[xnm]],na.rm=TRUE)
  
  # Make the initial (empty) plot
  plot(xlim,ylim,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="")
  #xx <- c(rep(-50,2),rep(10,2)); yy <- c(-4010,10,10,-4010)
  #g1 <- "grey95"; polygon(xx,yy,border=g1,col=g1)
  
  #xx <- c(rep(vline[1],2),rep(vline[2],2)); yy <- c(hline,hline[2:1])
  #polygon(xx,yy,border=NA,col=0)
  grid(col="grey80")
  
  # Loop over runs, but plot grey (bad) runs first
  for ( i in ii ) {
    if ( bad[i] > 0 ) lines(core[[xnm]][,i],depth,lwd=lwd[i],col=col[i])
  }
  
  for ( i in ii ) {
    if ( bad[i] == 0 ) lines(core[[xnm]][,i],depth,lwd=lwd[i],col=col[i])
  }

  #polygon(xx,yy,border=1,lwd=2)
  #abline(h=hline,v=vline,col=2,lwd=2)
  mtext(side=1,line=2,xlab,cex=1)
  mtext(side=2,line=2.5,ylab,cex=1,las=0)
  
  title(main=title)
  
  box()
}


plot2d <- function(mask,zs,mar=NA,plt=NA,title="",cex.title=1,asp=0.7,title.loc=c(0.9,0.15),lwd=1,
                   col=c("white","grey40","grey60"),add=TRUE,lab=TRUE,labcex=0.8,interp=FALSE,shade=FALSE)
{
  def.par <- par(no.readonly=TRUE)
  if ( is.na(mar[1]) ) mar <- def.par$mar
  if ( is.na(plt[1]) ) plt <- def.par$plt
  
  #par(new=add,mar=mar)
  par(new=add,plt=plt)
  
  if (interp != FALSE) {
    # tmp <- list(mask=mask,zs=zs)
    # new <- interp.clim(tmp,nms=c("mask","zs"),factor=2)
    # zs   <- new$zs
    # mask <- new$mask

    zs   = grid.interp(zs,factor=interp)
    mask = grid.interp(mask,factor=interp,mask=TRUE)

  }

  zs[mask==2] <- NA
  image(mask,axes=FALSE,col=col,asp=1/asp)

  #contour(mask,add=TRUE,levels=c(2),drawlabels=FALSE,col=darker(col[2],-5),lwd=lwd)

  contour(zs*1e-3,add=TRUE,levels=seq(from=0.0,to=3.5,by=0.5),lwd=lwd,
          drawlabels=FALSE,col="grey40" )
  
  contour(zs*1e-3,add=TRUE,levels=c(1.5,2.5),lwd=lwd*2.5,
          drawlabels=FALSE,col="grey40")

  #contour(mask,add=TRUE,levels=c(1),drawlabels=FALSE,col=alpha(1,50),lwd=lwd*0.5)
  
  if (shade) add_shade(zs)

  text(title.loc[1],title.loc[2],title,cex=cex.title,pos=2)
  
  if (add) par(new=add,plt=def.par$plt,mar=def.par$mar)  # return to normal
  
}

plot2d.ani <- function(dat, ...) 
{
  cat("Generating frames:",length(dat$time),"..")
  
  for ( i in 1:length(dat$time) ) {
    
    par(mar=c(0,0,0,0))
    plot2d(dat$mask[,,i],dat$zs[,,i],title=paste(dat$time[i],"ka"),add=FALSE )
    
    Sys.sleep(ani.options("interval"))
  }
  
  cat("done.\n")
  
}

plot.blank <- function(mar=NA)   
{ ## plot an empty space, in which text can be added
  def.par <- par(no.readonly=TRUE)
  if (is.na(mar[1])) mar <- c(1,1,1,1)
  par(mar=mar)
  plot(c(0,1),c(0,1),type="n",axes=FALSE, xlab="", ylab="")
  
  par(def.par$mar) # return to normal
}

