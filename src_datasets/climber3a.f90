module climber3a 

    use gridding_datasets
    use coordinates 
    use interp2D 
    use ncio 
    
    implicit none 

    private 
    public :: climber3a_atm_to_grid
    public :: climber3a_ocn_to_grid

    ! ocean directory: /home/jalvarez/outputs_clim3/ocean_marisa

contains 

    subroutine climber3a_atm_to_grid(outfldr,subfldr,grid,domain,path_in,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CLIMBER3-alpha atmospheric DATA
        !
        ! =========================================================
        implicit none 

        character(len=*) :: domain, outfldr, subfldr, path_in 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:), var(:,:)
            double precision, allocatable :: zs(:,:) 
            double precision :: lapse_ann    = 8.0d0 
            double precision :: lapse_summer = 6.5d0 
            double precision, allocatable :: var_hi(:,:)
            integer,          allocatable :: mask_hi(:,:)
        end type 

        type(inp_type)     :: inp
        type(grid_class)   :: grid0, grid0hi
        character(len=256) :: fldr_in, file_in_topo, file_in
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, np 

        type(map_class)  :: map, map0hi, maphi
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), outzs(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        ! Define the input filenames
        fldr_in      = trim(path_in)
        file_in_topo = trim(fldr_in)//"horo.cdf"
        file_in      = trim(fldr_in)//trim(domain)//".cdf"

        desc    = "CLIMBER-3alpha simulation output"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                              trim(grid%name)//"_"//trim(domain)//".nc"

        ! Load the domain information 
        nx = nc_size(file_in_topo,"XT_I")
        ny = nc_size(file_in_topo,"YT_J")
        np = nx*ny 

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))
        allocate(inp%zs(nx,ny))

        call nc_read(file_in_topo,"XT_I",inp%lon)
        call nc_read(file_in_topo,"YT_J",inp%lat)
        
        ! Define CLIMBER3a points and input variable field
        call grid_init(grid0,name="climber3a-atmos",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat )

        ! Define CLIMBER3a hi resolution intermediate interpolation grid
        call grid_init(grid0hi,name="climber3a-5deg",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x0=-180.d0,dx=5.d0,nx=73,y0=-90.d0,dy=5.d0,ny=37)
        call grid_allocate(grid0hi,inp%var_hi)
        call grid_allocate(grid0hi,inp%mask_hi)
        
        ! Define the variables to be mapped 
        allocate(vars(3))
        call def_var_info(vars( 1),trim(file_in),"TS_ANN","t2m_ann",units="Kelvin", &
                          long_name="Near-surface temperature (2-m), annual mean",method="quadrant")
        call def_var_info(vars( 2),trim(file_in),"TS_JJA","t2m_jja",units="Kelvin", &
                          long_name="Near-surface temperature (2-m), summer mean",method="quadrant")
        call def_var_info(vars( 3),trim(file_in),"PRC_ANN","pr_ann",units="mm*d**-1", &
                          long_name="Precipitation, annual mean",method="quadrant")

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        call grid_allocate(grid,outzs) 

        ! Initialize mappings
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)
        call map_init(map0hi,grid0,grid0hi,max_neighbors=10,lat_lim=8.d0,fldr="maps",load=.TRUE.)
        call map_init(maphi,grid0hi,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

!         ! Initialize mapping
!         call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Load reference topography in order to adjust temps to sea-level temps 
        call nc_read(file_in_topo,"HORO_PRESENT",inp%zs,missing_value=missing_value)
        write(*,*) "input zs : ", minval(inp%zs), maxval(inp%zs)

!         ! Map zs to new grid too
!         call map_field(map,"zs",inp%zs,outzs,outmask,"quadrant", &
!                           fill=.TRUE.,missing_value=missing_value)
        
        ! Map zs to new grid (two-step interpolation)
        call map_field(map0hi,var_now%nm_in,inp%zs,inp%var_hi,inp%mask_hi,"quadrant", &
                       missing_value=missing_value)
        call map_field(maphi, var_now%nm_in,inp%var_hi,outzs,outmask,"quadrant", &
                       missing_value=missing_value)

        ! Write output elevation to output file
        call nc_write(filename,"zs",real(outzs),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,"zs","units","m")
        call nc_write_attr(filename,"zs","long_name","Surface elevation")
        call nc_write_attr(filename,"zs","coordinates","lat2D lon2D")
            
        ! ## Map climatological gridded variables ##
        
        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            ! Read in current variable
            call nc_read(trim(var_now%filename),var_now%nm_in,inp%var,missing_value=missing_value)
            where(abs(inp%var) .ge. 1d10) inp%var = missing_value 

            ! Scale to sea-level temperature for interpolation
            if (trim(var_now%nm_out) .eq. "t2m_jja") &
                inp%var = inp%var + inp%lapse_summer*inp%zs 
            if (trim(var_now%nm_out) .eq. "t2m_ann") &
                inp%var = inp%var + inp%lapse_ann*inp%zs 

!             ! Map variable to new grid
!             call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
!                           fill=.TRUE.,missing_value=missing_value)
            
            ! Map variable to new grid (two-step interpolation)
            call map_field(map0hi,var_now%nm_in,inp%var,inp%var_hi,inp%mask_hi,"quadrant", &
                           missing_value=missing_value)
            call map_field(maphi, var_now%nm_in,inp%var_hi,outvar,outmask,"quadrant", &
                           missing_value=missing_value)

            ! Re-scale to near-surface temp for writing to file
            if (trim(var_now%nm_out) .eq. "t2m_jja") &
                outvar = outvar - inp%lapse_summer*outzs 
            if (trim(var_now%nm_out) .eq. "t2m_ann") &
                outvar = outvar - inp%lapse_ann*outzs 

            ! Write output variable to output file
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine climber3a_atm_to_grid

    subroutine climber3a_ocn_to_grid(outfldr,subfldr,grid,domain,path_in,max_neighbors,lat_lim)
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
        double precision :: lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:), z_ocn(:), depth(:)
            double precision, allocatable :: var(:,:)
            integer,          allocatable :: mask(:,:)
            double precision, allocatable :: var_hi(:,:)
            integer,          allocatable :: mask_hi(:,:)
        end type 

        type(inp_type)     :: inp
        type(grid_class)   :: grid0, grid0hi
        character(len=256) :: fldr_in, file_in 
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, nz 

        type(map_class)  :: map, map0hi, maphi
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        ! Define the input filenames
        fldr_in      = trim(path_in)
        file_in      = trim(fldr_in)//trim(domain)//".cdf"

        desc    = "CLIMBER-3alpha simulation output"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                              trim(grid%name)//"_"//trim(domain)//".nc"

        ! Load the domain information 
        nx = nc_size(file_in,"XT_I")
        ny = nc_size(file_in,"YT_J")
        nz = nc_size(file_in,"ZT_K")

        allocate(inp%lon(nx),inp%lat(ny),inp%z_ocn(nz),inp%depth(nz))

        call nc_read(file_in,"XT_I",inp%lon)
        call nc_read(file_in,"YT_J",inp%lat)
        call nc_read(file_in,"ZT_K",inp%depth)

        ! Make z negative and reverse it (rel to ocean surface)
        do k = 1, nz 
            inp%z_ocn(k) = -inp%depth(nz-k+1)
        end do  

        ! Define CLIMBER3a points and input variable field
        call grid_init(grid0,name="climber3a-ocn",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat )
        call grid_allocate(grid0,inp%var)
        call grid_allocate(grid0,inp%mask)

        ! Define CLIMBER3a hi resolution intermediate interpolation grid
        call grid_init(grid0hi,name="climber3a-2deg",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x0=-180.d0,dx=2.d0,nx=181,y0=-90.d0,dy=2.d0,ny=91)
        call grid_allocate(grid0hi,inp%var_hi)
        call grid_allocate(grid0hi,inp%mask_hi)
        
        ! Define the variables to be mapped 
        allocate(vars(2))
        call def_var_info(vars( 1),trim(file_in),"TEMP","to_ann",units="deg C", &
                          long_name="Potential temperature (annual mean)",method="quadrant")
        call def_var_info(vars( 2),trim(file_in),"mask","mask_ocn",units="1", &
                          long_name="Land-ocean mask (0=land, 1=ocean)",method="nn")

        ! Initialize mappings
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)
        call map_init(map0hi,grid0,grid0hi,max_neighbors=10,lat_lim=8.d0,fldr="maps",load=.TRUE.)
        call map_init(maphi,grid0hi,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x, units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y, units="kilometers")
        call nc_write_dim(filename,"z_ocn",x=inp%z_ocn,units="kilometers")

        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ## Map variable ##
        
        ! Map variable for each depth level
        do k = 1, nz 

            var_now = vars(1)

            ! Read in current variable (starting from last to reverse depth vector)
            call nc_read(var_now%filename,var_now%nm_in,inp%var,missing_value=missing_value, &
                         start=[1,1,nz-k+1],count=[nx,ny,1])
            where(abs(inp%var) .ge. 1d10) inp%var = missing_value 

            ! Perform two-step interpolation to higher resolution input grid,
            ! then to desired output grid 
            call map_field(map0hi,var_now%nm_in,inp%var,inp%var_hi,inp%mask_hi,"quadrant", &
                           missing_value=missing_value)
            call map_field(maphi, var_now%nm_in,inp%var_hi,outvar,outmask,"quadrant", &
                           missing_value=missing_value)

            ! Fill any additional missing values
            call fill_weighted(outvar,missing_value=missing_value)
            
            ! Clean up in case all values were missing and
            ! infinity values result from fill_weighted routine
            ! (eg, for deep bathymetry levels for GRL domain)
            where(outvar .ne. outvar) outvar = 1.d0

            ! Write output variable to output file
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="z_ocn", &
                          start=[1,1,k],count=[grid%G%nx,grid%G%ny,1])

            ! Also map mask 
            var_now = vars(2) 

            ! Define topo mask
            inp%mask = 1
            where(inp%var == missing_value) inp%mask = 0 

            call map_field(map,var_now%nm_in,dble(inp%mask),outvar,outmask,var_now%method, &
                          fill=.TRUE.,missing_value=missing_value)

            ! Write output mask to output file
            call nc_write(filename,var_now%nm_out,int(outvar),dim1="xc",dim2="yc",dim3="z_ocn", &
                          start=[1,1,k],count=[grid%G%nx,grid%G%ny,1])
        
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

    end subroutine climber3a_ocn_to_grid

end module climber3a 

