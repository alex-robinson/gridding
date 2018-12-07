module vavrus2018 

    use gridding_datasets
    use coord
    use ncio 
    use gaussian_filter
    
    implicit none 

    private 
    public :: vavrus2018_atm_to_grid
    public :: vavrus2018_ocn_to_grid

contains 

    subroutine vavrus2018_atm_to_grid(outfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       Vavrus et al (2018) MIS-19 climate DATA
        !
        ! =========================================================
        implicit none 

        character(len=*) :: domain, outfldr, path_in 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: sigma, lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:), var(:,:), var1D(:) 
            double precision :: lapse_ann    = 8.0d-3 
            double precision :: lapse_summer = 6.5d-3
        end type 

        type(inp_type)     :: inp
        type(grid_class)   :: grid0
        character(len=256) :: grid0_nm 
        character(len=256) :: file_in
        type(var_defs)     :: var_info
        integer :: nx, ny, np 

        type(map_class)  :: map
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        double precision, parameter :: sec_year = 365.0*24.0*3600.0

        desc    = "Vavrus et al. (2018) MIS-19 simulation output"
        ref     = "Vavrus et al.: Glacial Inception in Marine Isotope Stage 19: &
        &An Orbital Analog for a Natural Holocene Climate. Sci. Rep., 8:10213, &
        &doi:10.1038/s41598-018-28419-5."

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_CLIM-MIS19_V18.nc"

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)     
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",x=grid%G%y,units="kilometers")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ======================================================================
        ! Note: 
        ! Treat all variables separately, since grids for temp and precip and different
        ! topo and precip grid are the same, but consider separately for simplicity. 
        ! ======================================================================

        ! ======================================================================
        ! 1. Surface elevation 

        ! Define the variables to be mapped 
        call def_var_info(var_info,trim(file_in),"PHIS","z_srf",units="m",conv=real(1.0/9.81,kind(1d0)), &
                          long_name="Surface elevation",method="nng")

        ! Define input filename
        file_in = trim(path_in)//"/"//"PHIS.nc"

        ! Load domain information 
        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")
        grid0_nm = "vavrus2018_topo"
        np = nx*ny 

        ! Allocate input variable arrays
        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        ! Read domain variables 
        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat) 
        
        ! Define input grid based on input axes
        call grid_init(grid0,name=grid0_nm,mtype="latlon",units="degrees",lon180=.TRUE.,x=inp%lon,y=inp%lat)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        
        ! Load variable
        call nc_read(file_in,var_info%nm_in,inp%var,missing_value=mv) 
        write(*,*) "input : ", trim(var_info%nm_out), minval(inp%var), maxval(inp%var)

        ! Adjust units 
        inp%var = inp%var * var_info%conv 

        ! Map var to new grid
        call map_field(map,var_info%nm_out,inp%var,outvar,outmask,var_info%method, &
                       fill=.TRUE.,missing_value=mv,sigma=sigma)
        call nc_write(filename,var_info%nm_out,real(outvar),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,var_info%nm_out,"units",var_info%units_out)
        call nc_write_attr(filename,var_info%nm_out,"long_name",var_info%long_name)
        call nc_write_attr(filename,var_info%nm_out,"coordinates","lat2D lon2D")
        
        ! Also generate a land mask
        where(inp%var .gt. 1.d0) inp%var = 1.d0 
        where(inp%var .lt. 1.d0) inp%var = 0.d0 
        
        ! Map mask to new grid
        call map_field(map,"mask_land",inp%var,outvar,outmask,"nn", &
                          fill=.TRUE.,missing_value=mv,sigma=sigma)
        call nc_write(filename,"mask_land",nint(outvar),dim1="xc",dim2="yc")
           
        ! Write variable metadata
        call nc_write_attr(filename,"mask_land","units","1")
        call nc_write_attr(filename,"mask_land","long_name","Land mask (land=1)")
        call nc_write_attr(filename,"mask_land","coordinates","lat2D lon2D")


        ! == PRECIP (on same grid as topo) ==

        ! Define the variables to be mapped 
        call def_var_info(var_info,trim(file_in),"PRECT","pr",units="m/a",conv=sec_year, &
                          long_name="Precipitation",method="nng")

        ! Define input filename
        file_in = trim(path_in)//"/"//"PRECT.1371-1420_ave.nc"

        ! Load variable
        call nc_read(file_in,var_info%nm_in,inp%var,missing_value=mv) 
        write(*,*) "input : ", trim(var_info%nm_out), minval(inp%var), maxval(inp%var)

        ! Adjust units 
        inp%var = inp%var * var_info%conv 
        
        ! Map var to new grid
        call map_field(map,var_info%nm_out,inp%var,outvar,outmask,var_info%method, &
                       fill=.TRUE.,missing_value=mv,sigma=sigma)
        call nc_write(filename,var_info%nm_out,real(outvar),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,var_info%nm_out,"units",var_info%units_out)
        call nc_write_attr(filename,var_info%nm_out,"long_name",var_info%long_name)
        call nc_write_attr(filename,var_info%nm_out,"coordinates","lat2D lon2D")
        
        
        stop 

        ! ## Map climatological gridded variables ##
        
        var_now = var_info 

            ! Read in current variable
            call nc_read(trim(var_now%filename),var_now%nm_in,inp%var,missing_value=missing_value)
            where(abs(inp%var) .ge. 1d10) inp%var = missing_value 

!             ! Scale to sea-level temperature for interpolation
!             if (trim(var_now%nm_out) .eq. "t2m_sum") &
!                 inp%var = inp%var + inp%lapse_summer*inp%zs 
!             if (trim(var_now%nm_out) .eq. "t2m_ann") &
!                 inp%var = inp%var + inp%lapse_ann*inp%zs 

            ! Map variable to new grid
            call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
                          fill=.TRUE.,missing_value=missing_value,sigma=sigma)
            
!             ! Re-scale to near-surface temp for writing to file
!             if (trim(var_now%nm_out) .eq. "t2m_sum") &
!                 outvar = outvar - inp%lapse_summer*outzs 
!             if (trim(var_now%nm_out) .eq. "t2m_ann") &
!                 outvar = outvar - inp%lapse_ann*outzs 

            ! Write output variable to output file
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D") 

        return 

    end subroutine vavrus2018_atm_to_grid

    subroutine vavrus2018_ocn_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CLIMBER3-alpha oceanic DATA
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
            double precision, allocatable :: lon(:), lat(:), depth(:)
            double precision, allocatable :: var(:,:)
            integer,          allocatable :: mask(:,:)
        end type 

        type(inp_type)     :: inp
        type(grid_class)   :: grid0
        character(len=256) :: fldr_in, file_in 
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, nz 

        type(map_class)  :: map
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        character(len=256) :: c3_grid

        ! Define the input filenames
        fldr_in      = trim(path_in)
        file_in      = trim(fldr_in)//trim(domain)//".cdf"

        desc    = "CLIMBER-3alpha simulation output"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                              trim(grid%name)//"_"//trim(domain)//".nc"

        ! Load the domain information 
        if (nc_exists_var(file_in,"XT_I")) then 
            nx = nc_size(file_in,"XT_I")
            ny = nc_size(file_in,"YT_J")
            c3_grid = "climber3a-ocn"
        else 
            nx = nc_size(file_in,"lon")
            ny = nc_size(file_in,"lat")
            c3_grid = "climber3a-ocn-hi"
        end if 
        nz = nc_size(file_in,"ZT_K")

        allocate(inp%lon(nx),inp%lat(ny),inp%depth(nz))

        if (nc_exists_var(file_in,"XT_I")) then 
            call nc_read(file_in,"XT_I",inp%lon)
            call nc_read(file_in,"YT_J",inp%lat)
        else 
            call nc_read(file_in,"lon",inp%lon)
            call nc_read(file_in,"lat",inp%lat)
        end if
        call nc_read(file_in,"ZT_K",inp%depth)

        ! Define CLIMBER3a points and input variable field
        call grid_init(grid0,name="climber3a-ocn",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat )
        call grid_allocate(grid0,inp%var)
        call grid_allocate(grid0,inp%mask)

        ! Define the variables to be mapped 
        allocate(vars(2))
        call def_var_info(vars( 1),trim(file_in),"TEMP","to",units="degrees Celcius", &
                          long_name="Potential temperature (annual mean)",method="nng")
        call def_var_info(vars( 2),trim(file_in),"mask","mask_ocn",units="1", &
                          long_name="Land-ocean mask (0=land, 1=ocean)",method="nn")

        ! Initialize mappings
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x, units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y, units="kilometers")
        call nc_write_dim(filename,"depth",x=inp%depth,units="meters")

        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ## Map variable ##
        
        ! Map variable for each depth level
        do k = 1, nz 

            var_now = vars(1)

            ! Read in current variable (starting from last to reverse depth vector)
            call nc_read(var_now%filename,var_now%nm_in,inp%var,missing_value=mv, &
                         start=[1,1,k],count=[nx,ny,1])
            where(abs(inp%var) .ge. 1d10) inp%var = mv 

            call map_field(map, var_now%nm_in,inp%var,outvar,outmask,"nng", &
                           fill=.TRUE.,missing_value=mv,sigma=sigma)

            ! Clean up infinite values or all missing layers
            ! (eg, for deep bathymetry levels for GRL domain)
            where(outvar .ne. outvar .or. &
                  count(outvar.eq.mv) .eq. grid%npts) outvar = 1.d0

            ! Write output variable to output file
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="depth", &
                          start=[1,1,k],count=[grid%G%nx,grid%G%ny,1],missing_value=real(mv))

            ! Also map mask 
            var_now = vars(2) 

            ! Define topo mask
            inp%mask = 1
            where(inp%var == missing_value) inp%mask = 0 

            call map_field(map,var_now%nm_in,dble(inp%mask),outvar,outmask,var_now%method, &
                          fill=.TRUE.,missing_value=mv)

            ! Write output mask to output file
            call nc_write(filename,var_now%nm_out,int(outvar),dim1="xc",dim2="yc",dim3="depth", &
                          start=[1,1,k],count=[grid%G%nx,grid%G%ny,1],missing_value=int(mv))
        
        end do 

        ! Write variable metadata
        var_now = vars(1)
        call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
        call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
        call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        var_now = vars(2)
        call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
        call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
        call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        return 

    end subroutine vavrus2018_ocn_to_grid

end module vavrus2018 

