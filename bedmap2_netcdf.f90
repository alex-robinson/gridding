
! Bedmap2: surface, etc
!ncols         6667
!nrows         6667
!xllcorner     -3333500
!yllcorner     -3333500
!cellsize      1000
!NODATA_value  -9999

! Bedmap2: arthern accum.
!ncols         7899
!nrows         8300
!xllcorner     -3949500
!yllcorner     -3949500
!cellsize      1000
!NODATA_value  -9999

! Bedmap2: Rignot velocity
!ncols         5602
!nrows         5602
!xllcorner     -2800500
!yllcorner     -2801500
!cellsize      1000
!NODATA_value  -9999

program bedmap2_netcdf

    use ncio
    use coordinates 

    implicit none 

    double precision, parameter :: mv = -9999.0 

    type(grid_class) :: grid 

    double precision, allocatable, dimension(:) :: x, y
    double precision, allocatable, dimension(:,:) :: var0, var 
    integer :: i, j, h, k

    character(len=256) :: fnm, filename_topo, filename_acc, filename_vel 

    filename_topo = "output/Antarctica/ANT-1KM_BEDMAP2_topo.nc"
!     filename_acc  = "output/Antarctica/ANT-1KM_BEDMAP2_acc.nc"
!     filename_vel  = "output/Antarctica/ANT-1KM_BEDMAP2_vel.nc"

    ! Initialize the output grid
    call grid_init(grid,name="ANT-1KM",mtype="polar stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=1.d0,nx=6667,dy=1.d0,ny=6667, &
                       lambda=0.d0,phi=-90.d0,alpha=19.0d0)

    ! Allocate grid variable
    call grid_allocate(grid,var)

    ! Write grid information to output file
    if (.TRUE.) then

        fnm = filename_topo
        call nc_create(fnm)
            
        ! Add grid axis variables to netcdf file
        call nc_write_dim(fnm,"xc",x=grid%G%x,units=grid%units)
        call nc_write_dim(fnm,"yc",x=grid%G%y,units=grid%units)

        ! Add projection information if needed
        call nc_write_map(fnm,grid%mtype,grid%proj%lambda,phi=grid%proj%phi, &
                          x_e=grid%proj%x_e,y_n=grid%proj%y_n)
     
    !     call nc_write(fnm,"x2D",grid%x,dim1=xnm,dim2=ynm,grid_mapping=grid%name)
    !     call nc_write(fnm,"y2D",grid%y,dim1=xnm,dim2=ynm,grid_mapping=grid%name)
        call nc_write(fnm,"lon2D",real(grid%lon),dim1="xc",dim2="yc",grid_mapping=grid%name)
        call nc_write(fnm,"lat2D",real(grid%lat),dim1="xc",dim2="yc",grid_mapping=grid%name)

    end if 

    ! Allocate bedmap dimensions and data array
    call bedmap2_dims(x,y,var0,x0=-3333.d0,dx=1.d0,nx=6667,y0=-3333.d0,dy=1.d0,ny=6667)
    
    ! Surface elevation
    call bedmap2_read("data/Antarctica/Bedmap2/bedmap2_ascii/bedmap2_surface.txt", &
                      "zs",var0,missing_value=mv)
!     var = var0
    call bedmap2_align(var0=var0,x0=x(1),y0=y(1),var1=var,x1=grid%G%x(1),y1=grid%G%y(1),missing_value=mv)
    call nc_write(filename_topo,"zs",real(var),dim1="xc",dim2="yc",missing_value=real(mv), &
                  units="m",long_name="Surface elevation")
    
    ! Bedrock elevation
    call bedmap2_read("data/Antarctica/Bedmap2/bedmap2_ascii/bedmap2_bed.txt", &
                      "zb",var0,missing_value=mv)
    var = var0
    call nc_write(filename_topo,"zb",real(var),dim1="xc",dim2="yc",missing_value=real(mv), &
                  units="m",long_name="Bedrock elevation")
    
    ! Ice thickness
    call bedmap2_read("data/Antarctica/Bedmap2/bedmap2_ascii/bedmap2_thickness.txt", &
                      "H",var0,missing_value=mv)
    var = var0
    call nc_write(filename_topo,"H",real(var),dim1="xc",dim2="yc",missing_value=real(mv), &
                  units="m",long_name="Ice thickness")
    
    ! Ice-shelf mask
    call bedmap2_read("data/Antarctica/Bedmap2/bedmap2_ascii/bedmap2_icemask_grounded_and_shelves.txt", &
                      "mask_ice",var0,missing_value=mv)
    var = var0
    call nc_write(filename_topo,"mask_ice",int(var),dim1="xc",dim2="yc",missing_value=int(mv), &
                  units="-",long_name="Mask (ice-shelf)")
    
!     ! Rock mask
!     call bedmap2_read("data/Antarctica/Bedmap2/bedmap2_ascii/bedmap2_rockmask.txt", &
!                       "mask_rock",var0,missing_value=mv)
!     var = var0
!     call nc_write(filename_topo,"mask_rock",int(var),dim1="xc",dim2="yc",missing_value=int(mv), &
!                   units="-",long_name="Mask (rock)")
    
!     ! Vostock mask
!     call bedmap2_read("data/Antarctica/Bedmap2/bedmap2_ascii/bedmap2_lakemask_vostok.txt", &
!                       "mask_lake",var0,missing_value=mv)
!     var = var0
!     call nc_write(filename_topo,"mask_lake",int(var),dim1="xc",dim2="yc",missing_value=int(mv), &
!                   units="-",long_name="Mask (lake)")
    
    
    ! ==== Rignot velocities

    ! Bedmap2: Rignot velocity
    !ncols         5602
    !nrows         5602
    !xllcorner     -2800500
    !yllcorner     -2801500
    !cellsize      1000
    !NODATA_value  -9999

!     ! Allocate bedmap dimensions and data array
!     call bedmap2_dims(x,y,var0,x0=-2800.d0,dx=1.d0,nx=5602,y0=-2801.d0,dy=1.d0,ny=5602)
    
!     ! Surface elevation
!     call bedmap2_read("data/Antarctica/Bedmap2/bedmap2_ascii/rignot_velocity_bedmap2_grid.txt", &
!                       "uv",var0,missing_value=mv)
! !     var = var0
!     call bedmap2_align(var0=var0,x0=x(1),y0=y(1),var1=var,x1=grid%G%x(1),y1=grid%G%y(1),missing_value=mv)
!     call nc_write(filename_topo,"uv",real(var),dim1="xc",dim2="yc",missing_value=real(mv), &
!                   units="m**s^-1",long_name="Velocity")
    

contains 

    subroutine bedmap2_dims(x,y,var0,x0,dx,nx,y0,dy,ny)

        implicit none 

        double precision, allocatable :: x(:), y(:), var0(:,:)
        double precision :: x0, dx, y0, dy 
        integer          :: nx, ny 

        if (allocated(x))    deallocate(x)
        if (allocated(y))    deallocate(y)
        if (allocated(var0)) deallocate(var0)
        
        allocate(x(nx),y(ny),var0(nx,ny))

        do i = 1, nx
            x(i) = x0 + (i-1)*dx
        end do 
        do j = 1, ny
            y(j) = y0 + (j-1)*dx
        end do 

        write(*,*) "Bedmap2 input array information:"
        write(*,*) "nx, ny:   ", nx, ny 
        write(*,*) "range(x): ", minval(x), maxval(x)
        write(*,*) "range(y): ", minval(y), maxval(y)
        
        return 

    end subroutine bedmap2_dims

    subroutine bedmap2_read(filename,name,var2D,missing_value)

        implicit none 

        character(len=*) :: filename, name 
        double precision :: var2D(:,:), missing_value 
        character(len=50) :: tmpc 
        double precision  :: tmpd, missing_val0 
        integer :: i, nrow 

        write(*,*) "Reading: "//trim(filename)//": "//trim(name)

        ! Open the data file for reading
        open(166,file=filename,status="old")

        ! Read the header, extract the missing value from last line
        do i = 1, 6
            read(166,*) tmpc, tmpd 
        end do 
        missing_val0 = tmpd 

        ! Loop over all rows in file and store data in array
        nrow = size(var2D,2)
        var2D = missing_value
        do i = nrow, 1, -1 
            read(166,*) var2D(:,i)
        end do 

        ! Overwrite missing data with the desired missing value
        where(var2D .eq. missing_val0) var2D = missing_value 

        write(*,*) "Done."

        return

    end subroutine bedmap2_read

    subroutine bedmap2_align(var0,x0,y0,var1,x1,y1,missing_value)
        
        implicit none 

        double precision :: var0(:,:), var1(:,:)
        double precision :: x0, y0, x1, y1
        integer :: nx0,ny0, nx1, ny1, i0, j0, i, j 
        double precision :: missing_value 

        var1 = missing_value 

        nx0 = size(var0,1)
        ny0 = size(var0,2)
        nx1 = size(var1,1)
        ny1 = size(var1,2)

        if (x1 >= x0 .and. y1 >= y0) then 
            ! Handles case where destination grid starts after the original grid
            ! ie, x1 > x0 and y1 > y0
            i0 = nint(x1-x0)
            j0 = nint(y1-y0)

            do i = 1, nx0 
                if (i > nx1) exit 
                do j = 1, ny0
                    if (j > ny1) cycle 
                    var1(i,j) = var0(i0+i,j0+j)
                end do 
            end do 

        else 
            ! Handles case where destination grid starts before the original grid
            ! ie, x0 > x1 and y0 > y1 
            i0 = nint(x0-x1)
            j0 = nint(y0-y1)

            do i = 1, nx1 
                if (i > nx0) exit 
                do j = 1, ny1
                    if (j > ny0) cycle 
                    var1(i0+i,j0+j) = var0(i,j)
                end do 
            end do 

        end if 
        
        return

    end subroutine bedmap2_align 

end program 