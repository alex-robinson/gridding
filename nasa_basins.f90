

program nasa_basins 
    ! This program will take nasa input polygon borders
    ! from a file (lon,lat,basin) and interpolate it to
    ! a set of high res lon, lat points (lon,lat,basin)

    use polygons 

    implicit none 

    character(len=256) :: domain
    character(len=256) :: file_in, file_out, tmp 
    
    type basin_type 
        double precision, allocatable :: lon(:), lat(:), basin(:)
    end type 

    type(basin_type) :: inb, outb 
    integer :: nh, nl, np, npo, i, j, k  

    double precision :: dlon, dlat, lon0, lat0 
    integer :: nlat, nlon 

    ! ## USER DEFINED OPTIONS ##
    domain = "Antarctica" 
    dlon   = 0.5d0
    dlat   = 0.5d0 

    ! Get input data (from polygon files)
    if (trim(domain) .eq. "Antarctica") then 

        ! Define the input filenames
        file_in  = "data/Antarctica/nasa_basins/Ant_Full_DrainageSystem_Polygons.txt"
        file_out = "data/Antarctica/nasa_basins/Ant_Full_DrainageSystem_points.txt"

        nh = 7        ! Header length
        nl = 901329   ! File length 
        np = nl - nh  ! Number of data points 

        allocate(inb%lon(np),inb%lat(np),inb%basin(np))

        ! File format: lat, lon, basin 
        open(2,file=trim(file_in),status="old")
        do i = 1, nh
            read(2,"(a)") tmp 
        end do 
        do i = 1, np 
            read(2,*) inb%lat(i), inb%lon(i), inb%basin(i) 
        end do 
        close(2)

        write(*,*) "lon: ",minval(inb%lon),maxval(inb%lon)
        write(*,*) "lat: ",minval(inb%lat),maxval(inb%lat)
        write(*,*) "var: ",minval(inb%basin),maxval(inb%basin)

    else if (trim(domain) .eq. "Greenland") then 

        nh = 7 
        nl = 272972
        np = nl - nh  ! Number of data points 

        ! File format: basin, lat, lon 
    else 

        write(*,*) "nasa_basins:: error: "
        write(*,*) "Domain not recognized: ",trim(domain)
        stop 

    end if 

    ! Determine how many lon and lat values will be output
    nlon = int( ceiling( (maxval(inb%lon)-minval(inb%lon)) / dlon ) + 1 )
    nlat = int( ceiling( (maxval(inb%lat)-minval(inb%lat)) / dlat ) + 1 )
    npo  = nlon*nlat 

    ! Allocate output vectors 
    allocate(outb%lon(npo),outb%lat(npo),outb%basin(npo))

    write(*,*) "np, nlon, nlat, npo: ", np, nlon, nlat, npo

    ! Fill latlon values 
    lon0 = minval(inb%lon)
    lat0 = minval(inb%lat)

    k = 0 
    do i = 1, nlon 
        do j = 1, nlat 
            k = k+1 
            outb%lon(k) = lon0 + (i-1)*dlon
            outb%lat(k) = lat0 + (j-1)*dlat
        end do 
    end do 

    write(*,*) "lon: ",minval(outb%lon),maxval(outb%lon)
    write(*,*) "lat: ",minval(outb%lat),maxval(outb%lat)

    ! Go through each output point and determine if it fits inside a polygon 
    do k = 1, npo 


    end do 

! contains

end program nasa_basins 
