module badgeley2020

        use gridding_datasets
        use coord
        use ncio
        use gaussian_filter

        implicit none

        private
        public :: badgeley_to_grid

contains

        subroutine badgeley_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !      Badgeley et al., 2020 Retreat data
        !
        ! =========================================================

        implicit none

        character(len=*) :: domain, outfldr, subfldr, path_in
        type(grid_class) :: grid
        integer :: max_neighbors
        double precision :: sigma, lat_lim
        character(len=512) :: filename
        character(len=1024) :: desc, ref

        type inp_type
            double precision, allocatable :: lon(:), lat(:), var(:,:,:), var1d(:)
            double precision, allocatable :: time(:)
            double precision :: lapse = 8.0d-3
        end type

        type(inp_type)     :: inp
        type(grid_class)   :: grid0
        character(len=256) :: fldr_in, file_in_tas, file_in_pr
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, nt

        type(map_class)  :: map
        type(var_defs) :: var_now
        double precision, allocatable :: outvar(:,:), outzs(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var

        ! Define the input filenames
        fldr_in      = trim(path_in)
        file_in_tas  = trim(fldr_in)//"tas_main_Badgeley_etal_2020.nc"
        file_in_pr  = trim(fldr_in)//"pr_main_Badgeley_etal_2020.nc"

        desc    = "Badgeley et al. 2020 simulation output"
        ref     = "source: https://arcticdata.io/catalog/view/doi:10.18739/A2599Z26M"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_CLIM-RECON-B20.nc"

        ! Load the domain information 
        nx = nc_size(file_in_tas,"lon")
        ny = nc_size(file_in_tas,"lat")
        nt = nc_size(file_in_tas,"time")
        
        allocate(inp%time(nt))
        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny,nt))
        allocate(inp%var1d(nt))
        
        call nc_read(file_in_tas,"time",inp%time)
        call nc_read(file_in_tas,"lon",inp%lon)   
        call nc_read(file_in_tas,"lat",inp%lat)

        ! Define grid points and input variable field
        call grid_init(grid0,name="B20-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat)
        
        ! Define the variables to be mapped 
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in_tas),"time","time",units="ka",long_name="number of years before 1950 CE")
        call def_var_info(vars(2),trim(file_in_tas),"posterior_mean","tas",units="degrees Celsius",long_name="surface air temperature anomaly",method="nng")
        call def_var_info(vars(3),trim(file_in_pr),"posterior_mean","pr",units="",long_name="fractional precipitation",method="nng")

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)
        
        ! Initialize mappings
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"time", x=-inp%time,units="years")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)

        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)


        ! ##  Write the time as a 1d variable
        var_now = vars(1)

        call nc_read(trim(var_now%filename),var_now%nm_in,inp%var1d,missing_value=mv,start=[1],count=[nt])
        call nc_write(filename,"Time",-inp%var1d, dim1="time",start=[1], count=[nt])

        ! Write variable metadata
        call nc_write_attr(filename,"Time","units","years")
        call nc_write_attr(filename,"Time","long_name","time")


        ! ## Map climatological gridded variables ##
        
        ! Loop over variables
        do i = 2, size(vars)
            var_now = vars(i)

            ! Read in current variable
            call nc_read(trim(var_now%filename),var_now%nm_in,inp%var,missing_value=mv,start=[1,1,1],count=[nx,ny,nt])
       
            do k = 1, nt
            
                ! Map variable to new grid
                call map_field(map,var_now%nm_in,inp%var(:,:,k),outvar,outmask,"nn",fill=.TRUE.,missing_value=mv,sigma=sigma)
                
                outvar = outvar

                ! Smooth it out via sigma
                call filter_gaussian(var=outvar,sigma=sigma,dx=grid0%G%dx)
            
                ! Write output variable to output file
                call nc_write(filename,var_now%nm_out,real(outvar), &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,k],count=[grid%G%nx,grid%G%ny,1])
                
                write(*,"(a,f10.2,1x,a2)") trim(var_now%nm_out), inp%time(k), "years BP"
                   
            end do

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        end do


        end subroutine badgeley_to_grid

end module badgeley2020
