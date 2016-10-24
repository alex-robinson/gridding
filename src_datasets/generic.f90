module generic
    ! Use this module / routine to perform some generic interpolations 

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: generic_to_grid_nn
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
                           radius=dist_max*grid1%xy_conv,fill=.FALSE.,missing_value=mv)
            
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


end module generic
