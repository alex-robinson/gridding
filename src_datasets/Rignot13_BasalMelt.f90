module Rignot13_BasalMelt

    use gridding_datasets
    use coordinates 
    use planet 
    use oblimap_projection_module 
    use interp2D
    use ncio 
    
    implicit none 

    private 
    public :: Rignot13_BasalMelt_to_grid
    
contains 

    subroutine Rignot13_BasalMelt_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       BASAL MELTING DATA on bedmap2 grid
        !       Reference: 
        !           Rignot, E., Jacobs, S., Mouginot, J. and 
        !           Scheuchl, B., Ice-Shelf Melting Around Antarctica, 
        !           Science, doi: 10.1126/science.1235798, 2013
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename, infldr, prefix  
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid0
        character(len=256) :: file_in
        type(var_defs), allocatable :: vars(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:)
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp1(:,:), tmp2(:,:), tmp3(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: xax(:), yax(:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 

            ! Define the input filenames
            infldr  = "/data/sicopolis/data/Antarctica/"
            file_in = trim(infldr)//"Ant_MeltingRate.nc"

            desc    = "Ice shelf basal melting dataset"
            ref     = "Rignot, E., Jacobs, S., Mouginot, J. and Scheuchl, B.: &
                      &Ice-Shelf Melting Around Antarctica, Science, &
                      &doi: 10.1126/science.1235798, 2013."

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_BMELT-R13.nc"

            ! Define topography (BEDMAP2/rignot) grid and input variable field
            call grid_init(grid0,name="rignot-10KM",mtype="polar_stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-2800.d0,dx=10.d0,nx=561,y0=2800.d0,dy=-10.d0,ny=561, &
                   lambda=0.d0,phi=-90.d0,alpha=19.0d0)

        ! Define the variables to be mapped 
        allocate(vars(2))
        call def_var_info(vars(1),file_in,"melt_actual","bm_actual",units="m*a-1", &
                          long_name="Basal melt rate, actual present day")
        call def_var_info(vars(2),file_in,"melt_steadystate","bm_equil",units="m*a-1", &
                          long_name="Basal melt rate, shelf equilibrium")

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp1(5601,5601))  ! bedmap2-rignot array

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")

        ! Write the standard grid information (x2D,y2D,lon2D,lat2D,area,etc)
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write dataset meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ## INVARIANT FIELDS ##
        do i = 1, size(vars)
            var_now = vars(i) 
            call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=mv)
            call thin(invar,tmp1,by=10)
!             where( invar .eq. missing_value ) invar = 0.d0 
            where( invar .eq. 0.d0 ) invar = mv

            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,20.d3, &
                          fill=.TRUE.,missing_value=mv)
            call fill_mean(outvar,missing_value=mv)

            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine Rignot13_BasalMelt_to_grid



    subroutine map_field_conservative(grid0,grid,varname,var,missing_value)

        implicit none 

        type(grid_class), intent(IN)  :: grid0  ! Original grid information
        type(grid_class), intent(IN)  :: grid   ! New grid 
        character(len=*), intent(IN)  :: varname      ! Name of the variable being mapped
        double precision, intent(OUT) :: var(:,:)     ! Array of variable being mapped
        double precision, intent(IN)  :: missing_value  ! Points not included in mapping

        double precision :: x1, y1, x2, y2 
        integer :: i, j 

        ! Helpful fields available from grid_class object:
        ! grid%x    : 2D array of projected x-values [km]
        ! grid%y    : 2D array of projected y-values [km]
        ! grid%lon  : 2D array of lon-values [degrees]
        ! grid%lat  : 2D array of lat-values [degrees]
        ! grid%area : 2D array of grid-cell areas [m]
        ! grid%G%x  : vector of x-values that defines the x-axis [km]
        ! grid%G%y  : vector of y-values that defines the y-axis [km]
        ! grid%G%nx : length of x-axis 
        ! grid%G%ny : length of y-axis 

        ! Check if the grids are compatible for this mapping routine
        if (.not. same_projection(grid0%proj,grid%proj)) then 
            write(*,*) "map_field_conservative:: error:  &
                       &Currently this subroutine only interpolates between grids &
                       &on the same projection. Try again."
            stop 
        end if 

        ! From the coordinates package, there are some distance calculation routines
        ! already available:
        write(*,*) "Calculating distance..."

        x1 = grid0%G%x(1)
        y1 = grid0%G%y(1)
        x2 =  grid%G%x(1)
        y2 =  grid%G%y(1)
        write(*,*) "Distance [km]: ",x1, y1, " => ", x2, y2, " = ", cartesian_distance(x1,y1,x2,y2)


        ! A loop to get started 
        do i = 1, grid%G%nx 
            do j = 1, grid%G%ny 



            end do 
        end do 

        return 

    end subroutine map_field_conservative 


end module Rignot13_BasalMelt
