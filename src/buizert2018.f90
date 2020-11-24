module buizert2018

        use gridding_datasets
        use coord
        use ncio
        use gaussian_filter

        implicit none

        private
        public :: buizert_to_grid

contains

        subroutine buizert_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !      Buizert et al., 2018 Retreat data
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
            double precision, allocatable :: lon(:), lat(:), var(:,:,:,:), var1d(:)
            double precision, allocatable :: time(:)
            double precision, allocatable :: zs(:,:,:)
            double precision :: lapse = 8.0d-3
        end type

        type(inp_type)     :: inp
        type(grid_class)   :: grid0
        character(len=256) :: fldr_in, file_in
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, nm, nt

        type(map_class)  :: map
        type(var_defs) :: var_now
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var

        ! Define the input filenames
        fldr_in  = trim(path_in)
        file_in  = trim(fldr_in)//"GLand_22ka_recon_Buizert_20161228_anom.nc"

        desc    = "Buizert et al., 2018 simulation output"
        ref     = "source: https://www.ncdc.noaa.gov/paleo/study/23430"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_CLIM-RECON-B18.nc"

        ! Load the domain information 
        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")
        nm = nc_size(file_in,"month")
        nt = nc_size(file_in,"time")
        
        allocate(inp%time(nt))
        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny,nm,nt))
        allocate(inp%var1d(nt))
        allocate(inp%zs(nx,ny,nt))
        
        call nc_read(file_in,"time",inp%time)
        call nc_read(file_in,"lon",inp%lon)   
        call nc_read(file_in,"lat",inp%lat)

        ! Define grid points and input variable field
        call grid_init(grid0,name="B18-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat)
        
        ! Define the variables to be mapped 
        allocate(vars(4))
        call def_var_info(vars(1),trim(file_in),"time","time",units="Calendar years before present (present = 1950 C.E.)",long_name="Age")
        call def_var_info(vars(2),trim(file_in),"tas","tas",units="K",long_name="Surface air temperature at 2m height (anomaly)",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr","pr",units="m/s water equivalent",long_name="Precipitation fraction",method="nng")
        call def_var_info(vars(4),trim(file_in),"DEM","zs",units="m above sea level",long_name="Ice elevation")

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)
        
        ! Initialize mappings
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"month", x=1,nx=12,units="1")
        call nc_write_dim(filename,"time", x=inp%time,units="years")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)

        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)


        ! ## Write the time as a 1d variable
        var_now = vars(1)

        call nc_read(trim(var_now%filename),var_now%nm_in,inp%var1d,missing_value=mv,start=[1],count=[nt])
        call nc_write(filename,"Time",inp%var1d, dim1="time",start=[1], count=[nt])

        ! Write variable metadata
        call nc_write_attr(filename,"Time","units","years")
        call nc_write_attr(filename,"Time","long_name","time")


        ! ## Map climatological gridded variables ##
        
        ! Loop over variables
        do i = 2, size(vars) - 1
            var_now = vars(i)

            ! Read in current variable
            call nc_read(trim(var_now%filename),var_now%nm_in,inp%var,missing_value=mv,start=[1,1,1,1],count=[nx,ny,nm,nt])
       
            do k = 1, nt
                do m = 1, nm
                    ! Map variable to new grid
                    call map_field(map,var_now%nm_in,inp%var(:,:,m,k),outvar,outmask,"nn",fill=.TRUE.,missing_value=mv,sigma=sigma)
                
                    outvar = outvar

                    ! Smooth it out via sigma
                    call filter_gaussian(var=outvar,sigma=sigma,dx=grid0%G%dx)
            
                    ! Write output variable to output file
                    call nc_write(filename,var_now%nm_out,real(outvar), &
                            dim1="xc",dim2="yc",dim3="month",dim4="time",start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])
                
                end do

                write(*,"(a,f10.2,1x,a2)") trim(var_now%nm_out), inp%time(k), "years"   

            end do

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        end do


        ! ## Map topography 

        ! Read in topography
        var_now = vars(4)
        call nc_read(trim(var_now%filename),var_now%nm_in,inp%zs,missing_value=mv,start=[1,1,1],count=[nx,ny,nt])

        ! Write topography
        do k = 1, nt

            ! Map variable to new grid
            call map_field(map,var_now%nm_in,inp%zs(:,:,k),outvar,outmask,"nn",fill=.TRUE.,missing_value=mv,sigma=sigma)

            outvar = outvar

            ! Write output variable to output file
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="time", start=[1,1,k],count=[grid%G%nx,grid%G%ny,1])

            write(*,"(a,f10.2,1x,a2)") trim(var_now%nm_out), inp%time(k), "ka"
        end do

        ! Write variable metadata
        call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
        call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
        call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        end subroutine buizert_to_grid

end module buizert2018
