module generic
    ! Use this module / routine to perform some generic interpolations 

    use gridding_datasets
    use coord 

    implicit none 

    private 
    public :: generic_to_grid_nn
contains 

    subroutine generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int)
        ! Convert the variables to the desired grid format and write to file
        ! using nearest neighbor interpolation only 

        implicit none 
 
        type(grid_class), intent(IN) :: grid0, grid1 
        character(len=*), intent(IN) :: outfldr, dataset, path_in
        character(len=*), intent(IN) :: vname(:)        ! List of variable names to interpolate
        character(len=*), intent(IN) :: vname_int(:)    ! List of variable names that are integers

        ! Local variables
        character(len=512) :: filename
        type(map_class) :: map
        integer, parameter :: max_neighbors = 4 
        double precision :: lat_lim, dist_max 
        integer :: q, k 
        double precision, allocatable :: var0(:,:), var1(:,:)
        integer,          allocatable :: mask1(:,:)
        character(len=56) :: vname_now, units, long_name   
        logical :: is_int 

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid1%name)//"_"//trim(dataset)//".nc"


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

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid1%G%x,units=trim(grid1%units))
        call nc_write_dim(filename,"yc",   x=grid1%G%y,units=trim(grid1%units))
!         call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call grid_write(grid1,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
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
            call nc_read(trim(path_in),vname_now,var0,missing_value=mv)

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
