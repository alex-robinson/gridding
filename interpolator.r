# Math constants
torads = pi/180; todegs = 180.0/pi

# Earth: radii
R_earth = 6394800.0
a_earth = 6378137.0; b_earth = 6356752.3142
e_earth = sqrt((a_earth^2-b_earth^2)/(a_earth^2))

# Grid: reference for stereographic projection
lambda0 = -39.0 * torads;  phi0    =  71.0 * torads

sinphi0    = sin(phi0); sinlambda0 = sin(lambda0)
cosphi0    = cos(phi0); coslambda0 = cos(lambda0)
  
#reference x- and y-values for stereographic projection
x0 = -800.0; y0 = -3400.0    # Bottom-left corner
x1 =  700.0; y1 =  -600.0    # Upper-right corner

# Grid resolution
#dx = 20.0   # km
#nx = 76; ny = 141

get.xy <- function(dims=c(76,141),dx=20,xy0=NA)
# Generate x and y vectors for the domain of the grid (in km)
{

  nx <- dims[1];        ny <- dims[2]
  if (nx == 76) {
    dx = 20
  } else if (nx == 151) {
    dx = 10
  } else if (nx == 301) {
    dx = 5
  } else {
    cat("get.xy: No proper nx found.\n")
  }

  xlim <- dx*(nx-1)/2;  ylim <- dx*(ny-1)/2

  if (is.na(xy0[1])) {
    x0 <- -xlim; y0 <- -ylim
  } else {
    x0 <- xy0[1]; y0 <- xy0[2]
  }

  x <- seq(from=x0,by=dx,length.out=dims[1])
  y <- seq(from=y0,by=dx,length.out=dims[2])

  return(list(x=x,y=y))
}

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

### INTERPOLATION FUNCTIONS ###

calcdist <- function(p1,p2)
{   # Calculates the cartesian distance between two points
  
  x1 <- p1$x; y1 <- p1$y
  x2 <- p2$x; y2 <- p2$y
  
  out <- sqrt( (y2-y1)^2 + (x2-x1)^2 )
  
  return(out)
}

calcdistw <- function(p1,p2)
{   # Calculates the great circle distance between two points on the earth (lat/lon)
  
  phi1 <- p1$lat; lambda1 <- p1$lon
  phi2 <- p2$lat; lambda2 <- p2$lon
  
  # Convert the longitude to fit range -180:180
  if ( lambda1 > 180 ) { lambda1 <- lambda1-360 }
  if ( lambda2 > 180 ) { lambda2 <- lambda2-360 }
  
  f = 1 / 298.257223563

  L = abs(lambda2 - lambda1) * torads        # diff in longitude
  
  if (abs(phi1 - phi2) < 1e-12) {
    # Points are at the same latitude -> find arclength of latitude circle
    dist = L * R_earth*cos(phi1*torads)
  } else {
  
    #Vincenty’s formula as it is used in the script:
    # a_earth, b_earth = major & minor semiaxes of the ellipsoid   
    # f = flattening (a_earth−b_earth)/a_earth
    # φ1, φ2 = geodetic latitude   

    U1 = atan((1-f)*tan(phi1*torads))        # U is 'reduced latitude'
    U2 = atan((1-f)*tan(phi2*torads))  
    sinU1 = sin(U1); sinU2 = sin(U2)
    cosU1 = cos(U1); cosU2 = cos(U2)
    
    lmda = L; lmda_prime = 2.0 * pi
    
    i = 0
    while (abs(lmda - lmda_prime) > 1e-12) {
      sinsig = sqrt( (cosU2*sin(lmda))^2 + (cosU1*sinU2 - sinU1*cosU2*cos(lmda))^2 ) # eq 14
                  
      if (abs(sinsig - 0.0) < 1e-12) {
        cat("Error in calcdistw: using co-incident points.\n")
        cat("Lat/Lon 1:",phi1,lambda1,"\n")
        cat("Lat/Lon 2:",phi2,lambda2,"\n")
        cat("U1, U2, L:",U1,U2, L,"\n")
        cat("Using co-incident points. Try again#\n")
        break
      }
      cossig = sinU1*sinU2 + cosU1*cosU2*cos(lmda)             # eq 15
      sig = atan2(sinsig, cossig)                              # eq 16
      sinalph = cosU1*cosU2*sin(lmda) / sinsig
      cos2alph = 1 - sinalph^2                                 # trig identity; paragraph 6
      cos2sigm = cossig - 2.0*sinU1*sinU2 / cos2alph
      C = (f / 16.0) * cos2alph * (4.0 + f*(4-3*cos2alph))     # eq 10
      lmda_prime = lmda
      lmda = L + (1-C)*f*sinalph*(sig + C * sinalph *(cos2sigm + C * cossig*(-1+2.0*cos2sigm^2)))   # eq 11
                
      i = i + 1
      if (i == 1000) { cat("Too many iterations of Vincenty equation#\n"); dist <- NA; break }
    }
    
    uu2 = cos2alph * (a_earth^2 - b_earth^2) / b_earth^2
    AA = 1 + (uu2 / 16384.0) * (4096.0 + uu2*(-768.0 + uu2*(320.0 - 175.0*uu2)))   # eq 3
    BB = (uu2 / 1024.0) * (256.0 + uu2 * (-128.0 + uu2*(74.0 - 47.0*uu2)))         # eq 4
    dsig = BB * sinsig * (cos2sigm + (BB/4.0)*(cossig*(-1.0 + 2.0*cos2sigm^2) -
              (BB/6.0)*cos2sigm*(-3.0+4.0*sinsig^2)*(-3.0+4.0*cos2sigm^2)))        # eq 6
    
    dist = b_earth * AA * (sig-dsig)           # eq 19
    
  }
  
  return(dist)
    
}

interp2d <- function(obj,factor=2,m=FALSE)
{
  # Ranges of x and y values
  xlims <- range(obj$x);  ylims <- range(obj$y)
  
  
  # Length of incoming x,y vectors
  nx0 <- length(obj$x);   ny0 <- length(obj$y)
  
  # Length of new vectors
  nx1 <- (nx0-1)*factor+1; ny1 <- (ny0-1)*factor+1
  
  # Interval of new vectors
  dx <- (obj$x[2]-obj$x[1]) / factor
  
  # Generate the new object of the right dimensions
  obj2 <- list(x = seq(from=xlims[1],to=xlims[2],by=dx),
               y = seq(from=ylims[1],to=ylims[2],by=dx),
               z = matrix(NaN,nrow=nx1,ncol=ny1))
  
  # Perform the interpolation
  obj2 <- interp.surface.grid(obj,obj2)
  
  cat("Interp surface finished.\n")
  
  # Add a filter for whether it's a mask or not!!!
  # Temporary - could make land when ice+water touch!
  if (m) { obj2$z <- round(obj2$z,0)  }
  
  return(obj2)
}

interp.clim <- function(clima,nms=names(clima),factor=2)
{
  nvar <- length(nms)
  
  if ( "m2" %in% names(clima) )   mask <- clima$m2
  if ( "mask" %in% names(clima) ) mask <- clima$mask
  nx <- dim(mask)[1]
  ny <- dim(mask)[2]
  
  xy <- get.xy(dims=c(nx,ny))
  x <- xy$x;   y <- xy$y
  
  obj <- list(x=x,y=y)
  
  for (i in c(1:nvar)) {
  
    obj$z <- clima[[c(nms[i])]]
    
    if (nms[i] %in% c("m2","mask")) {
      obj2 <- interp2d(obj,factor,m=TRUE)
    } else {
      obj2 <- interp2d(obj,factor)
    }
    
    clima[[c(nms[i])]] <- obj2$z
  }
  
  cat("Finished interpolating clima.\n")
  
  return(clima)
}

idw <- function(var,dist,eps=1e-4,fac=8,dmax=NA,na.rm=FALSE,return.wts=FALSE)
{
  # Scale distance to normalize it
  if (is.na(dmax)) { dmax <- 1; dist <- dist/sum(dist) }
  
  # Make sure dist is small, but > 0
  dist[dist < eps] <- eps
  
  # Determine the weights and normalize
  #wts <- 1 / dist^2 
  wts <- exp(-(dist/dmax)*fac)
  
  if (na.rm) {
    wts[is.na(var)] <- 0
    var[is.na(var)] <- 0
  }
  
  # Normalize the weights
  wts <- wts / sum(wts)
  
  if ( !return.wts) {
    return( sum(var*wts) )   # Return the actual value of the variable at this point
  } else {
    return(wts)              # Just return the calculated weights
  }
}

### DETERMINE IDW interpolated field for a lower-resolution grid...
tolores <- function(old,new,nm="zs",dmax=NA,frac=20,latlon=FALSE)
{
  # old (0): hi-res; new (1): lo-res
  
  nx0 <- dim(old$xx)[1]; ny0 <- dim(old$xx)[2]
  dx0 <- diff(range(old$xx)) / (nx0-1)
  
  nx1 <- dim(new$xx)[1]; ny1 <- dim(new$xx)[2]
  dx1 <- diff(range(new$xx)) / (nx1-1)
  
  if (latlon) { dx0 <- 10; dx1 <- 100 }
  
  #if ( is.na(dmax) ) dmax <- dx1
  ratio <- dx1/dx0
  
  n <- ceiling(dx1/dx0)
  if (ratio==1) n <- 4      # If grid resolutions are the same, add some neighbors
  if (n %% 2 == 0) n <- n+1 # Make sure n is odd
         
  dist <- matrix(0,nrow=n,ncol=n)
  
  nn <- (n-1)/2
  
  for (ii in -nn:nn) {
    for (jj in -nn:nn) {
      i <- ii+nn+1; j <- jj+nn+1
      dist[i,j] <- sqrt(ii^2+jj^2)*dx0 
    }
  }
  
  # Now I have an interpolation matrix with distances, normalize
  # this using a weighting function
  ## Loop over weights to make central point have desired weight (eg 10%, 20%, 50%...)
  fac = 200*ratio
  for ( i in c(1:100) ) {
    weights <- as.matrix(idw(var=NA,dist=dist,dmax=dmax,fac=fac,return.wts=TRUE),nrow=n,ncol=n)
    diff1 <- 100*max(weights)-frac
    
    if ( abs( diff1 ) < 0.5 ) break   # If weighting is within 0.5%, stop iterating
    fac = round(fac-0.4*diff1,1)  # Increase or decrease fac to reach better fit
    fac = max(0,fac)              # Make sure fac is not less than zero
  }
  
  cat("iter",i,": fac",fac,": central point weighting",round(100*max(weights),0),"(goal:",frac,")\n")
  
  #weights <- as.matrix(idw(var=NA,dist=dist,dmax=dmax,fac=fac,return.wts=TRUE),nrow=n,ncol=n)
  wtts <<- weights
  
  # Determine number of variables to be interpolated
  nv <- length(nm)
  
  # Set up a lo-res matrix for output
  all <- list()
  for (q in 1:nv) {
  
    out <- matrix(NA,nrow=nx1,ncol=ny1) 
    in0 <- old[[nm[q]]]
    
    # Loop over the low-res matrix
    for ( i1 in 1:nx1 ) {
      for ( j1 in 1:ny1 ) {
        
        # Get the current index in the hi-res matrix
        i0 <- (i1-1)*ratio+1; j0 <- (j1-1)*ratio+1
        
        # Reset the neighbor array
        neighbors <- weights * 0
        ki <- 1; kj <- 1
        for ( ii in (i0-nn):(i0+nn)) {
          if ( ii < 1 ) ii <- 1; if ( ii >nx0 ) ii <- nx0
          kj <- 1
          for ( jj in (j0-nn):(j0+nn)) {
            if ( jj < 1 ) jj <- 1; if ( jj >ny0 ) jj <- ny0
            
              neighbors[ki,kj] <- in0[ii,jj]
            
            kj <- kj+1  
          }
          ki <- ki+1
        }
        
        out[i1,j1] <- sum(neighbors*weights)
      }
    }
  
    # Store the out matrix as the variable in an output array
    all[[nm[q]]] <- out
    cat("Aggregated",nm[q],"\n")
  }
  
  if (length(nm) == 1) {
   return(all[[1]])
  } else {
    return(all)
  }
}

tohires <- function(old,new,nm="zs")
{
  # old (0): lo-res; new (1): hi-res
  
  nx0 <- dim(old$xx)[1]; ny0 <- dim(old$xx)[2]
  dx0 <- diff(range(old$xx)) / (nx0-1)
  
  nx1 <- dim(new$xx)[1]; ny1 <- dim(new$xx)[2]
  dx1 <- diff(range(new$xx)) / (nx1-1)
  
  ratio <- dx1/dx0
  
  n <- ceiling(dx1/dx0); if (n %% 2 == 0) n <- n+1
  dist <- matrix(0,nrow=n,ncol=n)
  
  nn <- (n-1)/2
  
  # Set up a hi-res matrix for output
  out <- matrix(NA,nrow=nx1,ncol=ny1)
  
  # Loop over the hi-res matrix
  for ( i1 in 1:nx1 ) {
    for ( j1 in 1:ny1 ) {
      
      i0 <- floor((i1-1)/ratio)+1; j0 <- floor((j1-1)/ratio)+1
      i0n <- i0+1; j0n <- j0+1
      
      if (i0 > nx0)  i0  <- nx0; if (j0 > ny0)  j0  <- ny0
      if (i0n > nx0) i0n <- nx0; if (j0n > ny0) j0n <- ny0
      
      ifrac <- (((i1-1)/ratio+1) - i0) / (ratio+1)
      jfrac <- (((j1-1)/ratio+1) - j0) / (ratio+1)
      
      if (ifrac < 0 | jfrac < 0 |
          ifrac > 1 | jfrac > 1 ) {
        cat("Error, out of bounds...\n")
        cat("ifrac : ",ifrac,"\n","jfrac : ",jfrac,"\n")
        break
      } else {
        cat("ifrac : ",ifrac,"\n","jfrac : ",jfrac,"\n")
        cat("i0:i0n",i0,":",i0n,"\n")
        cat("j0:j0n",j0,":",j0n,"\n")
      }
      f0 <- old[[nm]][i0,j0] *(1-ifrac) + old[[nm]][i0n,j0n]*ifrac
      f1 <- old[[nm]][i0,j0n]*(1-ifrac) + old[[nm]][i0n,j0n]*ifrac
      
      out[i1,j1] <- f0*(1-jfrac) + f1*jfrac
    
    }
  }
  

  return(out)
}

interpw2 <- function(old,grid,nms=c("zs"),n=4,fac=6,latlon=TRUE)
{
  if ( latlon ) {
    p0 <- data.frame( lat=as.vector(old$lat), lon=as.vector(old$lon) )
    p1 <- data.frame( lat=as.vector(grid$lat),lon=as.vector(grid$lon) )
  } else {
    p0 <- data.frame( lat=as.vector(old$yy), lon=as.vector(old$xx) )
    p1 <- data.frame( lat=as.vector(grid$yy),lon=as.vector(grid$xx) )
  }
  
  n0 <- dim(p0)[1]; n1 <- dim(p1)[1]

  # Create a new list to hold new variables, populate with empty values
  new <- list(); nvar <- length(nms)
  for ( q in 1:nvar) {
    new[[nms[q]]] <- grid$lat*NA
  }
  
  checks <- c(1,seq(from=10,to=100,by=10))
  cat("Performing interpolations",n1,"...")
  
  wtts <<- numeric(0)

  # Loop over the new grid to get new values at each point
  for ( i in 1:n1 ) {
    
    per <- 100*i/(round(n1/100)*100)
    if ( per %in% checks) cat(per,"% ",sep="")
      
    # Get the weights and indices of neighbors
    p <- bilin(p1[i,],p0,n=n,latlon=latlon); ii <- p$ii
    wts <- idw(NA,dist=p$d,fac=fac,return.wts=TRUE)
    
    wtts <<- c(wtts,wts)
    
    # Loop over variable and scale each one by weights
    for ( q in 1:nvar) {
      nm <- nms[q]
      new[[nm]][i] <- sum(old[[nm]][ii]*wts)
      
#       if ( nm == "t2m" ) {
#         cat("wts:",round(wts,3),"\n")
#       }
    }

  }
  
  cat("\n")
  
  return(new)
}

interpw3 <- function(old,grid,nms=c("zs"),n=4,fac=6,latlon=TRUE,ii=NA,weights=NA,weights.out=FALSE)
{ ## 3d variables can be input (nx,ny,ntime)
  ## If 3rd dimension is present, loop over it
  
  if (is.na(ii[1])) ii <- c(1:length(old$lat))
  
  # Make sure some values exist
  if (!"xx" %in% names(grid)) grid$xx <- NA
  if (!"yy" %in% names(grid)) grid$yy <- NA
  
  if ( latlon ) {
    p0 <- data.frame( lat=as.vector(old$lat[ii]), lon=as.vector(old$lon[ii]) )
    p1 <- data.frame( lat=as.vector(grid$lat),lon=as.vector(grid$lon) )
  } else {
    p0 <- data.frame( lat=as.vector(old$yy[ii]), lon=as.vector(old$xx[ii]) )
    p1 <- data.frame( lat=as.vector(grid$yy),lon=as.vector(grid$xx) )
  }
  
  # Get lengths of vectors
  n0 <- dim(p0)[1]; n1 <- dim(p1)[1]
  
  
  if ( length(weights)==1 ) {  # weights==NA => perform the interpolations
  
    # Set up "percent complete" indicators
    checks <- c(1,seq(from=10,to=100,by=10))
    cat("Performing interpolations",n1,"...")
    
    p <- list()

    # Loop over the new grid to get new values at each point
    for ( i in 1:n1 ) {
      
      per <- 100*i/(round(n1/100)*100)
      if ( per %in% checks) cat(per,"% ",sep="")
        
      # Get the weights and indices of neighbors
      p[[i]]     <- bilin(p1[i,],p0,n=n,latlon=latlon,idw=fac)
          
    }
    
  } else {
  
    p <- weights
    cat("Weights obtained from input argument...")
  
  }
  
  if (weights.out) {
    return(p)
  } else {
    cat("\n","Now apply weights to variables...\n")
    
    # Create a new list to hold new variables, populate with empty/default values
  #   new <- list(lat=grid$lat,lon=grid$lon,xx=grid$xx,yy=grid$yy)
    nvar <- length(nms); nt <- rep(1,nvar)
    
    # Loop over variable names
    for ( q in 1:nvar) {
      nm <- nms[q]
      
      if ( length(dim(old[[nm]])) > 2 ) nt[q] <- dim(old[[nm]])[3]

      # Make sure old matrix is 3d (for consistency - it will convert back later)
      dim(old[[nm]]) <- c(dim(old[[nm]])[1:2],nt[q])
      
      cat(nm,"(",dim(old[[nm]]),")\n")
      
      # Make some empty variables
      new0 <- array(NA,dim=dim(grid$lat))
      new1 <- array(NA,dim=c(dim(new0),nt[q]))
    
      # Loop over time dimension
      for ( t in 1:nt[q] ) {
        
        # Store current 2d array (of points of interest)
        now <- old[[nm]][,,t][ii]
        
        # Loop over points in 2d array and scale each one by weights
        for ( i in 1:n1 ) {
          jj <- p[[i]]$ii; wts <- p[[i]]$wts
          new0[i] <- sum(now[jj]*wts)
        }
        
        # Store 2d array in 3d array
        new1[,,t] <- new0
        
      }
      
      # If only 2d variable, reduce to two dimensions
      if ( nt[q] == 1 ) dim(new1) <- dim(new1)[1:2]
      
      # Finally store 2d/3d array into output list
  #     new[[nm]] <- new1
      grid[[nm]] <- new1
      
    }
    
    cat("Done.\n")
  
    #   return(new)
    return(grid) 
  }

}

bilin <- function(pnow,p0,n=4,latlon=TRUE,idw=NA)
{
  
#  p01    p11
#     pnow
#  p00    p10

  latdiff <- pnow$lat - p0$lat
  londiff <- pnow$lon - p0$lon
  
  xydist  <- sqrt(latdiff^2 + londiff^2)
  
  # Get the closest n points
  ii <- order(xydist)[1:n]
  #ii <- which(rank(xydist) <= n)
  
  # Get the distances (in meters)
  dist <- xydist[ii]; k <- 1
  
  if ( latlon ) {
    for (i in ii) { 
      dist[k] <- calcdistw(pnow,list(lat=p0$lat[i],lon=p0$lon[i]))
      k <- k+1 
    }
  }
  
  d <- dist/sum(dist)
  
  wts <- NA
  if (!is.na(idw)) wts <- idw(NA,dist=d,fac=idw,return.wts=TRUE)
  
  # Return the indices of neighbors and the normalized weights
  return( list( ii=ii,d=d,wts=wts ) )

}


interpw <- function(new_in,old,nm=c("zs"),nmax=6,degmax=0.35,mask=TRUE,distance="latlon",
                    return.dist=FALSE)
{
  # Store the new grid dimensions
  nxy <- dim(new_in$lat)
  
  # Store the location information in separate vectors
  new <- new_in
  if (mask) {
    mnm = "m2"
    if ("mask" %in% names(new_in) ) { mnm <- "mask" }
    
    ii_mask <- which(new_in[[mnm]] <= 0 )
    for (q in 1:length(new_in) ) {
      new[[q]] <- new_in[[q]][ii_mask]
    }
    
  }
  
  if ( distance == "latlon" ) {
    p1 <- data.frame( lat=as.vector(new$lat),lon=as.vector(new$lon))
    p0 <- data.frame( lat=as.vector(old$lat),lon=as.vector(old$lon))
  
  } else if ( distance == "xy" ) {
  
    p1 <- data.frame( lat=as.vector(new$lat),lon=as.vector(new$lon),
                      x=as.vector(new$xx),y=as.vector(new$yy) )
    p0 <- data.frame( lat=as.vector(old$lat),lon=as.vector(old$lon),
                      x=as.vector(old$xx),y=as.vector(old$yy) )
  }
    
    
  n0 <- dim(p0)[1]
  n1 <- dim(p1)[1]
  
  # Store the variables we are interested in
  nvar <- length(nm)
  
  # Set up the output list
  out <- out.dist <- list()
  for (q in 1:nvar) { out[[q]] <- numeric(n1) }
  names(out) <- nm
  
  # Get distances of all points relative to the current point; loop over all points
  checks <- c(1,seq(from=10,to=100,by=10))
  cat("Performing interpolations",n1,"...")
  
  ns <- ns0 <- numeric(n1)*NA
  
  if ( distance=="latlon" ) {
  
    for (i1 in 1:n1) {
      
      per <- 100*i1/(round(n1/100)*100)
      if ( per %in% checks) cat(per,"% ",sep="")
      #cat(i1," ")
      
      # Store the current position
      pnow <- p1[i1,]
      
      # Determine indices that are within a certain distance (lat,lon)
      # to current point (so that not all great cicle distances for p0 have to be calculated#)
      degdist <- sqrt( (p0$lat-pnow$lat)^2 + (p0$lon-pnow$lon)^2 )
      degdist[degdist > degmax] <- NA; ii_sort <- which(!is.na(degdist))
      
      n <- length(ii_sort)
      ns0[i1] <- n
      
      if ( n == 0 ) {   # No neighbors within distance boundary, set variable values to NA
      
        for (q in 1:nvar) { out[[q]][i1] <- NA  }
        
      } else {     # Neighbors exist within distance boundary, calculate distances
        
        degdist1   <- degdist[ii_sort]
        ii_sorted1 <- sort.list(degdist1,method="quick",na.last=NA)
        n <- min(n,nmax); ii_sorted <- ii_sort[ii_sorted1][1:n]
        ns[i1] <- n  # For output of average neighbors

        dists <- numeric(n)
        for (i0 in 1:n) { inow <- ii_sorted[i0]; dists[i0] <- calcdistw(pnow,p0[inow,]) }
        
        # Store the indices and distances used for calculation
        out.dist[[i1]] <- data.frame(ii=ii_sorted,dist=dists)
        
        if ( !return.dist ) {
        
          # Call inverse distance weighting interpolation routine
          # for each variable of interest
          for (q in 1:nvar) { 
            var <- nm[q]; out[[q]][i1] <- idw(old[[var]][ii_sorted],dists[1:n],na.rm=TRUE)   
          }
        
        } else {
          
          tmp <- numeric(n)+1  # vector of 1's
          out.dist[[i1]]$wt.idw <- idw(tmp,dists[1:n],na.rm=TRUE,return.wts=TRUE)
          
        }
        
      }
      
    }
  
  } else if ( distance == "xy" ) {
    
    for (i1 in 1:n1) {
      
      per <- 100*i1/(round(n1/100)*100)
      if ( per %in% checks) cat(per,"% ",sep="")
      #cat(i1," ")
      
      # Store the current position
      pnow <- p1[i1,]
 
      # Calculate the distance based on projected x,y values
      # (only for values in the region)
      ii    <- which( abs(p0$lat-pnow$lat) <= degmax & abs(p0$lon-pnow$lon) <= degmax )
      dists <- numeric(dim(p0)[1])*NA
      dists[ii] <- calcdist(pnow,p0[ii,])
      
      # Determine indices of neighbors that are within boundary
      ii_sort <- which(!is.na(dists)); n <- ns0[i1] <- length(ii_sort)
      
      ii_sorted1 <- sort.list(dists,na.last=NA)
      n <- ns[i1] <- min(n,nmax); ii_sorted <- ii_sort[ii_sorted1][1:n]

      # Store the indices and distances used for calculation
      out.dist[[i1]] <- data.frame(ii=ii_sorted,dist=dists[ii_sorted])
      
      if ( !return.dist ) {
        
        # Call inverse distance weighting interpolation routine
        # for each variable of interest
        for (q in 1:nvar) { 
          var <- nm[q]; out[[q]][i1] <- idw(old[[var]][ii_sorted],dists[ii_sorted],na.rm=TRUE)   
        }
      
      } else {
          
        tmp <- numeric(n)+1  # vector of 1's
        out.dist[[i1]]$wt.idw <- idw(tmp,dists[ii_sorted],na.rm=TRUE,return.wts=TRUE)
        
      }
    
    }
  
  }
  
  cat("\n")
  cat("Average # neighbors     :",mean(ns0,na.rm=TRUE),"+/-",sd(ns0,na.rm=TRUE),"\n")
  cat("Average # neighbors used:",mean(ns,na.rm=TRUE), "+/-",sd(ns,na.rm=TRUE), "\n")
  
  if ( !return.dist ) {
    
    if (mask) {   # Reset out so that dimensions fit with the entire grid, fill in points that were interpolated
      tmp <- out
      for (q in 1:nvar) {
        var <- nm[q]
        out[[var]] <- new_in[[var]]*NA
        out[[var]][ii_mask] <- tmp[[var]]
      }
    }

    # Store original information about grid
    out$lat <- new_in$lat; out$lon <- new_in$lon
    out$xx  <- new_in$xx;  out$yy <- new_in$yy
    
    # Return vectors to desired array shape
    for (q in 1:length(out)) {
      dim(out[[q]]) <- nxy
    }
  
  } else {
    
    out <- out.dist
    
  }
  
  return(out)
}

### PROJECTION FUNCTIONS

proj_plane <- function(lon,lat,convert=1e3)
{
  lambda = lon*torads
  phi    = lat*torads

  e = e_earth
  mc  = cosphi0/sqrt(1.0-e*e*sinphi0*sinphi0)
  tc  = sqrt(((1.0-sinphi0)/(1.0+sinphi0))*((1.0+e*sinphi0)/(1.0-e*sinphi0))^e)

  t   = sqrt(((1.0-sin(phi))/(1.0+sin(phi)))*((1.0+e*sin(phi))/(1.0-e*sin(phi)))^e)
  rho = a_earth*mc*t/tc

  x =  rho*sin(lambda-lambda0)
  y = -rho*cos(lambda-lambda0)
  
  out <- data.frame(x=x/convert,y=y/convert,lon=lon,lat=lat)
}

#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Subroutine :  p l a n e _ e l l i p s o i d
#  Purpose    :  Computation of position (x,y) for latitude phi 
#                and longitude lambda in a polar stereographic
#                projection with reference to the WGS84 ellipsoid.
#  Author     :  Reinhard Calov
#  updated    :  Alex Robinson (17. Mar 2008)
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plane_ellipsoid <- function(coord=NA, xy=NA,convert=1e3)
{
  proj <- "plane"
  if (!is.na(coord[1])) { lambda = coord[1]*torads; phi = coord[2]*torads }
  if (!is.na(xy[1])) { x = xy[1]*convert; y = xy[2]*convert; proj <- "ellipsoid" }
  
  e = e_earth

  mc  = cosphi0/sqrt(1.0-e*e*sinphi0*sinphi0)
  tc  = sqrt(((1.0-sinphi0)/(1.0+sinphi0))*((1.0+e*sinphi0)/(1.0-e*sinphi0))^e)

  if ( proj == "plane" ) {
    # Perform stereographic projection: from ellipsoid to plane

    t   = sqrt(((1.0-sin(phi))/(1.0+sin(phi)))*((1.0+e*sin(phi))/(1.0-e*sin(phi)))^e)
    rho = a_earth*mc*t/tc

    x =  rho*sin(lambda-lambda0)
    y = -rho*cos(lambda-lambda0)
  
  } else if ( proj == "ellipsoid" ) {
    # Perform inverse projection: from plane to ellipsoid

    rho = sqrt(x^2+y^2)
    t   = rho*tc/(a_earth*mc)

    lambda=lambda0+atan(x/(-y))

    #  fix point iteration
    phi_p = 0.5*pi-2.0*atan(t)
    l     = 0
    eps   = 3600.0
    while(eps >= 1e-9) {
      l     = l+1
      phi   = 0.5*pi-2.0*atan(t*((1.0-e*sin(phi_p))/(1.0+e*sin(phi_p)))^(0.5*e))
      eps   = abs(phi-phi_p)
      phi_p = phi
    }
    
  }
  
  out <- data.frame(x=x/convert,y=y/convert,lat=phi*todegs,lon=lambda*todegs)
  return(out)
  
}
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Subroutine :  g e o _ c o o r d
#  Purpose    :  Computation of position (x,y) for latitude phi 
#                and longitude lambda in a polar stereographic
#                projection.
#  Author     :  Ralf Greve
#  updated    :  Alex Robinson (06. Jan 2010)
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
geo_coord <- function(coord=NA, xy=NA,convert=1e3)
{
  
  proj <- "plane"
  if (!is.na(coord[1])) { lambda = coord[1]*torads; phi = coord[2]*torads }
  if (!is.na(xy[1])) { x = xy[1]*convert; y = xy[2]*convert; proj <- "ellipsoid" }
  
  fac <- 0.5                     # Northern hemisphere
  if (phi0 < 0) fac <- -fac   # Southern hemisphere

  K = (cos(0.25*pi+fac*phi0))^2
  
  phi    <- fac*pi-2.0*atan(sqrt(x^2+y^2)/(2.0*R_earth*K))
  lambda <- lambda0 + atan2(y,x)+0.5*pi

  # Adjust lambda so that it has the appropriate range
#   if (lambda < 0.0) {
#     lambda <- lambda + 2.0*pi
#   } else if (lambda >= (2.0*pi)) {
#     lambda <- lambda - 2.0*pi
#   }
#   
#   if (lambda < -pi) lambda = lambda + 2*pi
  
  out <- data.frame(x=x/convert,y=y/convert,lat=phi*todegs,lon=lambda*todegs)
  return(out)
}

creategrid <- function(x=NULL,y=NULL,x0=NULL,nx=NULL,y0=NULL,ny=NULL,dx=NULL,convert=1e3,latlon=TRUE)
{        
  cat("Creating polar stereographic / cartesian grid...")
  
  
  if (is.null(x)) x <- seq(from=x0,by=dx,length.out=nx)
  if (is.null(y)) y <- seq(from=y0,by=dx,length.out=ny)
  
  nx = length(x)
  ny = length(y) 

  xx <- yy <- lats <- lons <- matrix(NA,nrow=nx,ncol=ny)
  
  coord <- list(lat=NA,lon=NA)
  
  # Perform transformation to sphere/ellipsoid from plane
  # for all points on the grid
  for (j in 1:ny) {
    for ( i in 1:nx) {
      
      if (latlon) coord <- plane_ellipsoid(xy=c(x[i],y[j]),convert=convert)   # Slower
      #coord <- geo_coord(xy=c(x[i],y[j]),convert=convert)          # Faster, broken?
      
      lats[i,j] <- coord$lat; lons[i,j] <- coord$lon
      xx[i,j] <- x[i]; yy[i,j] <- y[j]
    }
  }
  
  out <- list(x=x,y=y,lat=lats,lon=lons,xx=xx,yy=yy)
  
  cat("done\n")
  
  return(out)
}

creategrid.latlon <- function(lon0,nx,lat0,ny,dx,dy=dx,convert=1e3)
{        
  cat("Creating lat/lon grid...")
  
  
  lon <- seq(from=lon0,by=dx,length.out=nx)
  lat <- seq(from=lat0,by=dy,length.out=ny)
  
  xx <- yy <- lats <- lons <- matrix(NA,nrow=nx,ncol=ny)
  
  # Perform transformation from sphere/ellipsoid to plane
  # for all points on the grid
  for (j in 1:ny) {
    for ( i in 1:nx) {
      coord <- plane_ellipsoid(coord=c(lon[i],lat[j]))
      lats[i,j] <- coord$lat; lons[i,j] <- coord$lon
      xx[i,j] <- coord$x; yy[i,j] <- coord$y
    }
  }
  
  out <- list(xx=xx/convert,yy=yy/convert,lat=lats,lon=lons)
  
  cat("done\n")
  
  return(out)
}


