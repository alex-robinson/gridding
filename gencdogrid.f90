program gencdogrid
    ! This program will produce a grid description file
    ! for use with `cdo remapcon` commands etc. 
    ! However, it can be generated automatically using 
    ! `cdo -griddes infile` if the projection is well
    ! defined. 
    
    use coord 
    use ncio 
    use control 
    
    implicit none 

    integer :: nx, ny 
    real(4), allocatable :: lon2D(:,:) 
    real(4), allocatable :: lat2D(:,:) 
    real(4), allocatable :: lon(:) 
    real(4), allocatable :: lat(:) 
    logical :: is_projection 

    character(len=512) :: filename_in, filename_out

    !filename_in  = "ice_data/Greenland/GRL-32KM/GRL-32KM_REGIONS.nc"
    !filename_out = "ice_data/Greenland/GRL-32KM/GRL-32KM_cdodesc.txt"

    filename_in  = "/Users/robinson/GoogleDriveUCM/wrk/mypapers/pmip-mis3/wrk/latlon-05deg_REGIONS.nc"
    filename_out = "/Users/robinson/GoogleDriveUCM/wrk/mypapers/pmip-mis3/wrk/grid_latlon-05deg.txt"

    is_projection = .TRUE. 

    if (is_projection) then 

        nx = nc_size(filename_in,"xc")
        ny = nc_size(filename_in,"yc")
    
        allocate(lon2D(nx,ny))
        allocate(lat2D(nx,ny))
    
        call nc_read(filename_in,"lon2D",lon2D)
        call nc_read(filename_in,"lat2D",lat2D)
    
    else 
        ! latlon grid 

        nx = nc_size(filename_in,"lon")
        ny = nc_size(filename_in,"lat")
    
        allocate(lon2D(nx,ny))
        allocate(lat2D(nx,ny))
        allocate(lon(nx),lat(ny))

        call nc_read(filename_in,"lon",lon)
        call nc_read(filename_in,"lat",lat)
    
        call gen_latlon2D(lon2D,lat2D,lon,lat)

    end if 

    call write_cdo_gridfile(lon2D,lat2D,filename_out)

    write(*,*) "Grid description file written: "//trim(filename_out)

contains 
    
    subroutine write_cdo_gridfile(lon2D,lat2D,filename)

        implicit none 

        real(4), intent(IN) :: lon2D(:,:) 
        real(4), intent(IN) :: lat2D(:,:) 
        character(len=*), intent(IN) :: filename 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, jm1, ip1, jp1 
        integer :: fnum
        real(4) :: bnds(4) 

        fnum = 98 

        nx = size(lon2D,1)
        ny = size(lon2D,2)

        open(fnum,file=filename,status='unknown',action='write')

        write(fnum,"(a)")     "gridtype = curvilinear"
        write(fnum,"(a,i10)") "gridsize = ", nx*ny 
        write(fnum,"(a,i10)") "xsize    = ", nx
        write(fnum,"(a,i10)") "ysize    = ", ny

        ! x values 
        write(fnum,*) ""
        write(fnum,"(a)") "# Longitudes"
        write(fnum,"(a)") "xvals = "
        do j = 1, ny 
            write(fnum,"(50000f10.3)") lon2D(:,j)
        end do 

        write(fnum,*) ""
        write(fnum,"(a)") "# Longitudes of cell corners"
        write(fnum,"(a)") "xbounds = "
        do j = 1, ny 
        do i = 1, nx 

            im1 = max(1,i-1)
            jm1 = max(1,j-1)
            ip1 = min(nx,i+1)
            jp1 = min(ny,j+1)

            ! Determine bounds (lower-right, upper-right, upper-left, lower-left)
            ! ie, get ab-nodes from aa-nodes
            bnds(1) = 0.25*(lon2D(i,j)+lon2D(ip1,j)+lon2D(i,jm1)+lon2D(ip1,jm1))
            bnds(2) = 0.25*(lon2D(i,j)+lon2D(ip1,j)+lon2D(i,jp1)+lon2D(ip1,jp1))
            bnds(3) = 0.25*(lon2D(i,j)+lon2D(im1,j)+lon2D(i,jp1)+lon2D(im1,jp1))
            bnds(4) = 0.25*(lon2D(i,j)+lon2D(im1,j)+lon2D(i,jm1)+lon2D(im1,jm1))
            
            write(fnum,"(4f10.3)") bnds 

        end do 
        end do 

        ! y values 
        write(fnum,*) ""
        write(fnum,"(a)") "# Latitudes"
        write(fnum,"(a)") "yvals = "
        do j = 1, ny 
            write(fnum,"(50000f10.3)") lat2D(:,j)
        end do 

        write(fnum,*) ""
        write(fnum,"(a)") "# Latitudes of cell corners"
        write(fnum,"(a)") "ybounds = "
        do j = 1, ny 
        do i = 1, nx 

            im1 = max(1,i-1)
            jm1 = max(1,j-1)
            ip1 = min(nx,i+1)
            jp1 = min(ny,j+1)

            ! Determine bounds (lower-right, upper-right, upper-left, lower-left)
            ! ie, get ab-nodes from aa-nodes
            bnds(1) = 0.25*(lat2D(i,j)+lat2D(ip1,j)+lat2D(i,jm1)+lat2D(ip1,jm1))
            bnds(2) = 0.25*(lat2D(i,j)+lat2D(ip1,j)+lat2D(i,jp1)+lat2D(ip1,jp1))
            bnds(3) = 0.25*(lat2D(i,j)+lat2D(im1,j)+lat2D(i,jp1)+lat2D(im1,jp1))
            bnds(4) = 0.25*(lat2D(i,j)+lat2D(im1,j)+lat2D(i,jm1)+lat2D(im1,jm1))
            
            write(fnum,"(4f10.3)") bnds 

        end do 
        end do 

        close(fnum)

        return 

    end subroutine write_cdo_gridfile

    subroutine gen_latlon2D(lon2D,lat2D,lon,lat)

        implicit none 

        real(4), intent(OUT) :: lon2D(:,:) 
        real(4), intent(OUT) :: lat2D(:,:) 
        real(4), intent(IN)  :: lon(:) 
        real(4), intent(IN)  :: lat(:) 

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(lon,1)
        ny = size(lat,1)

        do j = 1, ny 
            lon2D(:,j) = lon 
        end do 

        do i = 1, nx 
            lat2D(i,:) = lat 
        end do 

        return 

    end subroutine gen_latlon2D

end program gencdogrid