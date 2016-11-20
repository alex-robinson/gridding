module generic
    ! Use this module / routine to perform some generic interpolations 

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: generic_to_grid_nn
    public :: generic_to_grid_3D_nn
    public :: find_nearest_grid
    public :: nearest_to_grid 

contains 

    subroutine generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac)
        ! Convert the variables to the desired grid format and write to file
        ! using nearest neighbor interpolation only 

        implicit none 
 
        type(grid_class), intent(IN) :: grid0, grid1 
        character(len=*), intent(IN) :: outfldr, dataset, path_in
        character(len=*), intent(IN) :: vname(:)        ! List of variable names to interpolate
        character(len=*), intent(IN) :: vname_int(:)    ! List of variable names that are integers
        integer, intent(IN), optional :: thin_fac

        ! Local variables
        character(len=512) :: filename
        type(map_class) :: map
        integer, parameter :: max_neighbors = 4 
        double precision :: lat_lim, dist_max 
        integer :: q, k 
        double precision, allocatable :: var0(:,:), var1(:,:), tmp(:,:)
        integer,          allocatable :: mask1(:,:)
        character(len=56) :: vname_now, units, long_name   
        logical :: is_int 
        integer :: thin_by, nx00, ny00  

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid1%name)//"_"//trim(dataset)//".nc"

        ! Determine whether the input variables should be thinned
        ! to match the grid0 definition (default: no thinning)
        thin_by = 1 
        if (present(thin_fac)) thin_by = thin_fac 

        ! First get dist_max in meters
        dist_max = 2.0d0* max(sqrt((grid0%G%dx*grid0%xy_conv)**2+(grid0%G%dy*grid0%xy_conv)**2), &
                              sqrt((grid1%G%dx*grid1%xy_conv)**2+(grid1%G%dy*grid1%xy_conv)**2))
        
        ! Determine lat_lim assuming 1deg = 100km
        lat_lim  = dist_max / 100d3 

        ! Convert dist_max to units of target grid 
        dist_max = dist_max / grid1%xy_conv 


        ! Initialize mapping
        call map_init(map,grid0,grid1,max_neighbors=max_neighbors, &
                      lat_lim=lat_lim,dist_max=dist_max,fldr="maps",load=.TRUE.)


        ! Allocate the input and output grid variables
        call grid_allocate(grid0,var0)
        call grid_allocate(grid1,var1)
        call grid_allocate(grid1,mask1)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        if (thin_by .gt. 1) then 
            nx00 = (grid0%G%nx-1)*thin_by+1
            ny00 = (grid0%G%ny-1)*thin_by+1
            allocate(tmp(nx00,ny00))
        else 
            call grid_allocate(grid0,tmp)
        end if 

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid1%G%x,units=trim(grid1%units))
        call nc_write_dim(filename,"yc",   x=grid1%G%y,units=trim(grid1%units))
        call grid_write(grid1,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        write(*,*) "dim(tmp): ", size(tmp,1), size(tmp,2)

        ! Loop over the variables and perform the interpolation 
        do q = 1, size(vname)

            ! Get current variable name
            vname_now   = vname(q)

            ! Check if writing integer or float
            is_int = .FALSE.
            do k = 1, size(vname_int)
                if (trim(vname_int(k)) .eq. trim(vname_now)) is_int = .TRUE. 
            end do 

            ! Read in the original variable
            call nc_read(trim(path_in),vname_now,tmp,missing_value=mv)

            if (thin_by .gt. 1) then 
                if (is_int .or. thin_by .eq. 2) then 
                    call thin(var0,tmp,by=thin_by,missing_value=mv)
                else 
                    call thin_ave(var0,tmp,by=thin_by,missing_value=mv)
                end if 
            else 
                var0 = tmp 
            end if 

            write(*,*) "range(var):         ", minval(tmp), maxval(tmp)
            if (thin_by .gt. 1) &
                write(*,*) "range(var_thinned): ", minval(var0), maxval(var0)
            
            ! Map the variable
            var1 = mv 
            call map_field(map,vname_now,var0,var1,mask1,method="nn", &
                           radius=dist_max,fill=.FALSE.,missing_value=mv)
            
            if (is_int) then 
                call nc_write(filename,vname_now,nint(var1),dim1="xc",dim2="yc",missing_value=int(mv))
            else
                call nc_write(filename,vname_now,real(var1),dim1="xc",dim2="yc",missing_value=real(mv))
            end if 
            
            ! Transfer attributes too 
            call nc_read_attr(trim(path_in),vname_now,"units",units)
            call nc_read_attr(trim(path_in),vname_now,"long_name",long_name)

            call nc_write_attr(filename,vname_now,"units",units)
            call nc_write_attr(filename,vname_now,"long_name",long_name)

        end do 



        return 

    end subroutine generic_to_grid_nn

    subroutine generic_to_grid_3D_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac)
        ! Convert the variables to the desired grid format and write to file
        ! using nearest neighbor interpolation only 

        implicit none 
 
        type(grid_class), intent(IN) :: grid0, grid1 
        character(len=*), intent(IN) :: outfldr, dataset, path_in
        character(len=*), intent(IN) :: vname(:)        ! List of variable names to interpolate
        character(len=*), intent(IN) :: vname_int(:)    ! List of variable names that are integers
        integer, intent(IN), optional :: thin_fac

        ! Local variables
        character(len=512) :: filename
        type(map_class) :: map
        integer, parameter :: max_neighbors = 4 
        double precision :: lat_lim, dist_max 
        integer :: q, k 
        double precision, allocatable :: var0(:,:), var1(:,:), tmp(:,:)
        integer,          allocatable :: mask1(:,:)
        character(len=56) :: vname_now, units, long_name   
        logical :: is_int 
        integer :: thin_by, nx00, ny00  

        character(len=56), allocatable :: dnames(:) 
        integer, allocatable :: dims(:)
        character(len=56) :: units3
        real(4), allocatable :: x3(:) 

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid1%name)//"_"//trim(dataset)//".nc"

        ! Determine whether the input variables should be thinned
        ! to match the grid0 definition (default: no thinning)
        thin_by = 1 
        if (present(thin_fac)) thin_by = thin_fac 

        ! First get dist_max in meters
        dist_max = 2.0d0* max(sqrt((grid0%G%dx*grid0%xy_conv)**2+(grid0%G%dy*grid0%xy_conv)**2), &
                              sqrt((grid1%G%dx*grid1%xy_conv)**2+(grid1%G%dy*grid1%xy_conv)**2))
        
        ! Determine lat_lim assuming 1deg = 100km
        lat_lim  = dist_max / 100d3 

        ! Convert dist_max to units of target grid 
        dist_max = dist_max / grid1%xy_conv 


        ! Initialize mapping
        call map_init(map,grid0,grid1,max_neighbors=max_neighbors, &
                      lat_lim=lat_lim,dist_max=dist_max,fldr="maps",load=.TRUE.)


        ! Allocate the input and output grid variables
        call grid_allocate(grid0,var0)
        call grid_allocate(grid1,var1)
        call grid_allocate(grid1,mask1)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        if (thin_by .gt. 1) then 
            nx00 = (grid0%G%nx-1)*thin_by+1
            ny00 = (grid0%G%ny-1)*thin_by+1
            allocate(tmp(nx00,ny00))
        else 
            call grid_allocate(grid0,tmp)
        end if 

        ! Determine 3rd dimension name and size 
        ! Note: assumes all variables have same 3rd dimension
        call nc_dims(path_in,vname(1),dnames,dims)
        call nc_read_attr(path_in,trim(dnames(3)),"units",units3)
        allocate(x3(dims(3)))
        call nc_read(path_in,trim(dnames(3)),x3)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid1%G%x,units=trim(grid1%units))
        call nc_write_dim(filename,"yc",   x=grid1%G%y,units=trim(grid1%units))
        call nc_write_dim(filename,trim(dnames(3)),x=x3,units=trim(units3))
        call grid_write(grid1,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        write(*,*) "dim(tmp): ", size(tmp,1), size(tmp,2)

        ! Loop over the variables and perform the interpolation 
        do q = 1, size(vname)

            ! Get current variable name
            vname_now   = vname(q)

            ! Check if writing integer or float
            is_int = .FALSE.
            do k = 1, size(vname_int)
                if (trim(vname_int(k)) .eq. trim(vname_now)) is_int = .TRUE. 
            end do 

            ! Loop over 3rd dimension to map variable 
            do k = 1, dims(3)

                ! Read in the original variable
                call nc_read(trim(path_in),vname_now,tmp,missing_value=mv, &
                             start=[1,1,k],count=[dims(1),dims(2),1])

                if (thin_by .gt. 1) then 
                    if (is_int .or. thin_by .eq. 2) then 
                        call thin(var0,tmp,by=thin_by,missing_value=mv)
                    else 
                        call thin_ave(var0,tmp,by=thin_by,missing_value=mv)
                    end if 
                else 
                    var0 = tmp 
                end if 

                write(*,*) "range(var):         ", minval(tmp), maxval(tmp)
                if (thin_by .gt. 1) &
                    write(*,*) "range(var_thinned): ", minval(var0), maxval(var0)
                
                ! Map the variable
                var1 = mv 
                call map_field(map,vname_now,var0,var1,mask1,method="nn", &
                               radius=dist_max,fill=.FALSE.,missing_value=mv)
                
                if (is_int) then 
                    call nc_write(filename,vname_now,nint(var1),dim1="xc",dim2="yc",dim3=trim(dnames(3)),missing_value=int(mv), &
                                  start=[1,1,k],count=[grid1%G%nx,grid1%G%ny,1])
                else
                    call nc_write(filename,vname_now,real(var1),dim1="xc",dim2="yc",dim3=trim(dnames(3)),missing_value=real(mv), &
                                  start=[1,1,k],count=[grid1%G%nx,grid1%G%ny,1])
                end if 
            
            end do ! end 3rd dimension loop 

            ! Transfer attributes too 
            call nc_read_attr(trim(path_in),vname_now,"units",units)
            call nc_read_attr(trim(path_in),vname_now,"long_name",long_name)

            call nc_write_attr(filename,vname_now,"units",units)
            call nc_write_attr(filename,vname_now,"long_name",long_name)

        end do 



        return 

    end subroutine generic_to_grid_3D_nn

    subroutine find_nearest_grid(ii,jj,x,y,xout,yout,latlon,max_dist,lat_lim)
        ! Return the indices (i,j) of the x(nx), y(ny)
        ! grid point that is closest to the desired xout/yout values
        ! Distances measured in [m]

        implicit none 

        integer, intent(INOUT) :: ii(:,:), jj(:,:) 
        real(4), intent(IN)  :: x(:), y(:)
        real(4), intent(IN)  :: xout(:,:), yout(:,:) 
        logical, intent(IN) :: latlon 
        real(4), intent(IN), optional :: max_dist 
        real(4), intent(IN), optional :: lat_lim 
        
        ! Local variables 
        real(4) :: max_distance
        real(4) :: lat_limit 
        integer :: i0, j0, i1, j1, k   
        real(4) :: dist, dist_min 
        real(4) :: xout_now, yout_now 
        integer :: nx1, ny1 
        integer :: j00, j01 

        integer :: dists_nx
        integer, allocatable :: inds_x(:)
        real(4), allocatable :: dists_x(:) 

        ! Planet parameters ("WGS84")
        real(8), parameter :: a = 6378137.0d0             ! Equatorial ellipsoid radius,   a in Snyder, WGS84
        real(8), parameter :: f = 1.0d0/298.257223563d0   ! Flattening of the ellipsoid
                
        nx1 = size(xout,1)
        ny1 = size(xout,2)

        lat_limit = 10.0 ! 10 degrees by default 
        if (present(lat_lim)) lat_limit = lat_lim 

        ! Confirm maximum distance of interest 
        max_distance = 1e10
        if (present(max_dist)) max_distance = max_dist 

        ! Set up a filter for the x-distances, only sample every nth point first
        ! then zoom in to x-points with small distances 
        dists_nx = max(2, floor((size(x)-1)/10.0))
        allocate(inds_x(dists_nx),dists_x(dists_nx))

        k = 0 
        do i0 = 1, size(x), max(1,int(size(x)/dists_nx))
            k = k+1
            inds_x(k) = i0
            if (k .eq. dists_nx) exit 
        end do 

        if (k .lt. dists_nx) then 
            write(*,*) "inds_x not filled yet."
            write(*,*) dists_nx, k, i0, int((size(x)-1)/10.0), size(x)
            stop 
        end if 

        ! Summary of low-res calc setup 
        write(*,*) "Low res x: nx = ", dists_nx, " / ", size(x) 
        write(*,*) "   range(inds): ", minval(inds_x), maxval(inds_x)
        
        ! Initialize to missing indices everywhere
        ii = -1
        jj = -1 

        do j1 = 1, ny1

            do i1 = 1, nx1
            
                ! Loop over target grid and find all nn indices 

                ! Define current target point of interest
                xout_now = xout(i1,j1)
                yout_now = yout(i1,j1)

                ! Reset the minimum distance
                dist_min = 1e10 

                ! Find index of nearest latitude row
                j0 = minloc(abs(yout_now-y),1)

!                 ! Loop over grid and find nearest neighbor indices 
!                 do j0 = 1, size(y)
                    
                    if (abs(yout_now-y(j0)) .lt. lat_limit) then 
                        ! Only check here, if the y-point is within range 

                        ! First check filtered x-points 
                        do i0 = 1, size(inds_x)
                            k = inds_x(i0)

                            if (latlon) then
                                ! Use planetary (latlon) values
                                dists_x(i0) = planet_distance(a,f,x(k),y(j0),xout_now,yout_now)

                            else
                                ! Use cartesian values to determine distance
                                dists_x(i0) = cartesian_distance(x(k),y(j0),xout_now,yout_now)                    

                            end if 

                        end do 


                        do i0 = 1, size(x)

                            k = minloc(abs(inds_x-i0),1)
                            if (dists_x(k) .lt. 5.0*max_distance) then 
                                ! Only check x-points with low-res distances in reasonable range

                                if (latlon) then
                                    ! Use planetary (latlon) values
                                    dist = planet_distance(a,f,x(i0),y(j0),xout_now,yout_now)

                                else
                                    ! Use cartesian values to determine distance
                                    dist = cartesian_distance(x(i0),y(j0),xout_now,yout_now)                    

                                end if 

                                if (dist .lt. dist_min .and. dist .lt. max_distance) then 
                                    ii(i1,j1) = i0 
                                    jj(i1,j1) = j0 
                                    dist_min = dist
                                end if 

                            end if 

                        end do 

                    end if 

!                 end do 

            end do 

            ! Output every 5th column to check progress
            if (j1 .eq. 1 .or. mod(j1,5) .eq. 0 .or. j1 .eq. ny1) then 
                write(*,"(a,i10,a3,i12,a5,g12.4)") "  ",j1, " / ", ny1,"   : ", dist_min
            end if 
        end do 

        return 

    end subroutine find_nearest_grid

    subroutine nearest_to_grid(zout,z,ii,jj)

        implicit none 

        real(4), intent(INOUT) :: zout(:,:)
        real(4), intent(IN) :: z(:,:)
        integer, intent(IN) :: ii(:,:), jj(:,:)  
        
        ! Local variables 
        integer :: i, j, inow, jnow 
        real(4) :: xout,  yout 

        ! Loop over target grid and fill in available nearest neighbors
        zout = mv 
 
        do i = 1, size(zout,1) 
        do j = 1, size(zout,2) 

            inow = ii(i,j)
            jnow = jj(i,j) 

            ! Only update output array if valid neighbor was found 
            if (inow .gt. 0 .and. jnow .gt. 0) then 
                zout(i,j) = z(inow,jnow)
            end if 

        end do
        end do 

        return 

    end subroutine nearest_to_grid


!     subroutine var_to_grid(grid0,grid1,path_in,path_out,vdef,vtype,method,thin_fac)
!         ! Convert the variables to the desired grid format and write to file
!         ! using nearest neighbor interpolation only 

!         implicit none 
 
!         type(grid_class), intent(IN) :: grid0, grid1 
!         character(len=*), intent(IN) :: outfldr, dataset, path_in
!         type(var_defs),   intent(IN) :: vdef                       ! Input/output variable definitions
!         character(len=*), intent(IN) :: vtype                      ! int, real, dble
!         character(len=*), intent(IN), optional :: method ! Interpolation method 
!         integer, intent(IN), optional :: thin_fac

!         ! Local variables
!         type(map_class) :: map
!         character(len=56) :: interp_method 
!         integer, parameter :: max_neighbors = 100 
!         double precision :: lat_lim, dist_max, dmax0, dmax1 
!         integer :: q, k 
!         double precision, allocatable :: var0(:,:), var1(:,:), tmp(:,:)
!         integer,          allocatable :: mask1(:,:)
!         logical :: is_int 
!         integer :: thin_by, nx00, ny00  

!         ! Determine interpolation method to be used 
!         interp_method = "con"   ! Conservative by default 
!         if (present(method)) interp_method = trim(method)
!         if (trim(vtype) .eq. "int") interp_method = "nn"    ! Ensure consistency

!         ! Determine whether the input variables should be thinned
!         ! to match the grid0 definition (default: no thinning)
!         thin_by = 1 
!         if (present(thin_fac)) thin_by = thin_fac 

!         ! First get dist_max in meters 
        
!         if (grid0%is_cartesian) then 
!             dmax0 = sqrt((grid0%G%dx*grid0%xy_conv)**2+(grid0%G%dy*grid0%xy_conv)**2)
!         else 
!             ! dx is in degrees, roughly convert to meters (for lat_lim)
!             dmax0 = sqrt((grid0%G%dx)**2+(grid0%G%dy)**2)*100d3
!         end if 

!         if (grid1%is_cartesian) then 
!             dmax1 = sqrt((grid1%G%dx*grid1%xy_conv)**2+(grid1%G%dy*grid1%xy_conv)**2)
!         else 
!             ! dx is in degrees, roughly convert to meters (for lat_lim)
!             dmax1 = sqrt((grid1%G%dx)**2+(grid1%G%dy)**2)*100d3
!         end if 

!         ! dist_max is the higher dmax calculated for each grid (multiplied by a factor)
!         dist_max = 2.0d0* max(dmax0,dmax1)

!         ! Determine lat_lim from dist_max, assuming 1deg = 100km
!         lat_lim  = dist_max / 100d3 

!         ! Convert dist_max to units of target grid 
!         dist_max = dist_max / grid1%xy_conv 

!         ! Initialize mapping
!         call map_init(map,grid0,grid1,max_neighbors=max_neighbors, &
!                       lat_lim=lat_lim,dist_max=dist_max,fldr="maps",load=.TRUE.)

!         ! Allocate the input and output grid variables
!         call grid_allocate(grid0,var0)
!         call grid_allocate(grid1,var1)
!         call grid_allocate(grid1,mask1)

!         ! Allocate tmp array to hold full data 
!         ! (that will be trimmed to smaller size if thin_by > 1)
!         if (thin_by .gt. 1) then 
!             nx00 = (grid0%G%nx-1)*thin_by+1
!             ny00 = (grid0%G%ny-1)*thin_by+1
!             allocate(tmp(nx00,ny00))
!         else 
!             call grid_allocate(grid0,tmp)
!         end if 

!         ! Read in the original variable
!         call nc_read(trim(path_in),vdef%nm_in,tmp,missing_value=mv)

!         if (thin_by .gt. 1) then 
!             if (trim(interp_method) .eq. "nn" .or. thin_by .eq. 2) then 
!                 call thin(var0,tmp,by=thin_by,missing_value=mv)
!             else 
!                 call thin_ave(var0,tmp,by=thin_by,missing_value=mv)
!             end if 
!         else 
!             var0 = tmp 
!         end if 

!         ! Make sure units match output units 
!         where (var0 .ne. mv) var0 = var0*vdef%conv 

!         ! Map the variable
!         var1 = mv 

!         if (trim(interp_method) .eq. "con") then 
!             call map_field_conservative_map1(map%map,vdef%nm_out,var0,var1)
!         else 
!             call map_field(map,vdef%nm_out,var0,var1,mask1,method=interp_method, &
!                            radius=dist_max,fill=.FALSE.,missing_value=mv)
!         end if 

!         ! Write the variable to the file 
!         select case(trim(vtype))

!             case("int") 
!                 call nc_write(path_out,vdef%nm_out,nint(var1),dim1="xc",dim2="yc",missing_value=int(mv))
            
!             case("real")
!                 call nc_write(path_out,vdef%nm_out,real(var1),dim1="xc",dim2="yc",missing_value=real(mv))

!             case("dble")
!                 call nc_write(path_out,vdef%nm_out,real(var1),dim1="xc",dim2="yc",missing_value=real(mv))

!             case DEFAULT 
!                 write(*,*) "var_to_grid:: error: vtype must be one of (int,real,dble): ", trim(vtype)
!                 stop 

!         end select 

!         ! Write attributes too 
!         call nc_write_attr(filename,vdef%nm_out,"units",    vdef%units_out)
!         call nc_write_attr(filename,vdef%nm_out,"long_name",vdef%long_name)

!         return 

!     end subroutine var_to_grid


end module generic
