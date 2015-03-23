module ETOPO 

    use gridding_datasets
    use coordinates 
    use ncio 
    
    implicit none 

    private 
    public :: etopo1_to_grid
    
contains 

    subroutine etopo1_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       ETOPO1 DATA
        !
        ! =========================================================
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:)
            double precision, allocatable :: var(:,:)
        end type 

        type(inp_type)     :: inp
        integer :: nx, ny 
        type(grid_class)   :: grid0
        character(len=256) :: fldr_in, file_in_1, file_in_2 
        type(var_defs), allocatable :: vars(:)

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        ! Define the input filenames
        fldr_in    = "/data/sicopolis/data/ETOPO/"
        file_in_1  = trim(fldr_in)//"ETOPO1-ICE-020deg.nc"
        file_in_2  = trim(fldr_in)//"ETOPO1-BED-020deg.nc"

        desc    = "ETOPO1 present-day Earth topography data"
        ref     = "http://ngdc.noaa.gov/mgg/global/global.html"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_TOPO-ETOPO1.nc"

        ! Get the input dimensions
        nx = nc_size(file_in_1,"longitude")
        ny = nc_size(file_in_1,"latitude")
        allocate(inp%lon(nx),inp%lat(ny))
        call nc_read(file_in_1,"latitude",inp%lon)
        call nc_read(file_in_1,"latitude",inp%lat)
        inp%lat = inp%lat(ny:1)

        ! Define the input grid
        call grid_init(grid0,name="ETOPO1-020deg",mtype="latlon",units="degrees",lon180=.TRUE., &
                         x=inp%lon,y=inp%lat)
        
        ! Allocate the input array
        call grid_allocate(grid0,inp%var)

        ! Define the variables to be mapped 
        allocate(vars(2))
        call def_var_info(vars( 1),trim(file_in_1),"z","zs",units="m", &
                          long_name="Surface elevation",method="nn")
        call def_var_info(vars( 2),trim(file_in_2),"z","zb",units="m", &
                          long_name="Bedrock elevation",method="nn")

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)     

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            ! Read in current variable
            call nc_read(trim(var_now%filename),var_now%nm_in,inp%var,missing_value=missing_value)
            where(abs(inp%var) .ge. 1d8) inp%var = missing_value 

            ! Map variable to new grid
            call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
                          fill=.TRUE.,sigma=40.d0,missing_value=missing_value)

            ! Write output variable to output file
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine etopo1_to_grid


end module ETOPO 

