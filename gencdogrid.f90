program gencdogrid

    use coord 
    use ncio 
    use control 
    
    implicit none 

    integer :: nx, ny 
    real(4), allocatable :: lon2D(:,:) 
    real(4), allocatable :: lat2D(:,:) 
    
    character(len=512) :: filename_in, filename_out

    filename_in  = "ice_data/Greenland/GRL-32KM/GRL-32KM_REGIONS.nc"
    filename_out = "ice_data/Greenland/GRL-32KM/GRL-32KM_cdodesc.txt"

    nx = nc_size(filename_in,"xc")
    ny = nc_size(filename_in,"yc")
    
    allocate(lon2D(nx,ny))
    allocate(lat2D(nx,ny))
    
    call nc_read(filename_in,"lon2D",lon2D)
    call nc_read(filename_in,"lat2D",lat2D)
    
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

end program gencdogrid