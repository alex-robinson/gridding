

program nasa_basins 
    ! This program will take nasa input polygon borders
    ! from a file (lon,lat,basin) and interpolate it to
    ! a set of high res lon, lat points (lon,lat,basin)

    use polygons 

    implicit none 

    character(len=256) :: domain
    character(len=256) :: file_in, file_out, tmp 
    
    type basin_type 
        real(4), allocatable :: lon(:), lat(:), basin(:)
    end type 

    type(basin_type) :: inb, outb 
    integer :: nh, nl, np, npo, i, j, k, k0, k1 

    real(4) :: dlon, dlat, lon0, lat0  
    integer :: nlat, nlon 

    real(4), allocatable :: basins(:) 
    integer, allocatable :: inds(:)
    integer :: nb, q 
    logical :: in_basin 

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

    ! Determine the input basins available 
    call unique(basins,inb%basin)
    nb = size(basins)

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

        ! Initialize basin info for current point
        in_basin      = .FALSE.
        outb%basin(k) = 0.0 

        ! Loop over basins and check point in polygon
        do q = 1, nb 
            call which(inb%basin==basins(q),inds)
            in_basin = point_in_polygon(outb%lon(k), outb%lat(k), inb%lon(inds), inb%lat(inds))
            if (in_basin) exit 
        end do 

        ! If basin was found, save it
        if (in_basin) outb%basin(k) = basins(q)

        if (mod(k,1000) .eq. 0) & 
            write(*,*) minval(outb%basin(1:k)), maxval(outb%basin(1:k)), k, "/", npo
    end do 

    ! Write output points to file 
    ! File format: lon, lat, basin 
    open(2,file=trim(file_out),status="unknown")
    write(*,"(3a12)") "lon", "lat", "basin"
    do i = 1, npo 
        write(2,*) outb%lon(i), outb%lat(i), outb%basin(i) 
    end do 
    close(2)


contains

    subroutine which(x,ind)

        implicit none 

        logical :: x(:)
        integer, allocatable :: ind(:) 
        integer :: n 

        n = count(x)
        if (allocated(ind)) deallocate(ind)
        allocate(ind(n))

        n = 0
        do i = 1, size(x) 
            if (x(i)) then 
                n = n+1
                ind(n) = i 
            end if
        end do 

        return 

    end subroutine which 

!     subroutine get_polygon(xx,yy,bb,b)

!         implicit none 


!         return 

!     end subroutine get_polygon 

    subroutine unique(xu,x)
        ! Return only the unique values of a vector
        ! http://rosettacode.org/wiki/Remove_duplicate_elements#Fortran

        implicit none 

        real(4) :: x(:)          ! The input
        real(4) :: res(size(x))  ! The unique values
        real(4), allocatable :: xu(:)  ! The output 
        integer :: k                   ! The number of unique elements
        integer :: i, j
        real(4), parameter :: tol = 1d-5
        logical :: found 

        k = 1
        res(1) = x(1)
        do i=2,size(x)
            found = .FALSE.
            do j=1,k
                if (abs(res(j)-x(i)) .le. tol) then 
                   ! Found a match so start looking again
                   found = .TRUE. 
                   cycle 
                end if
            end do
            ! No match found so add it to the output
            if (.not. found) then 
                k = k + 1
                res(k) = x(i)
            end if 
        end do

        ! Store output in properly sized output vector
        if(allocated(xu)) deallocate(xu)
        allocate(xu(k))
        xu = res(1:k)

        return 

    end subroutine unique 

end program nasa_basins 
