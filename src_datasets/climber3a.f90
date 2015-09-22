module climber3a 

    use gridding_datasets
    use coordinates 
    use interp2D 
    use gaussian_filter
    use ncio 
    
    implicit none 

    private 
    public :: climber3a_atm_to_grid
    public :: climber3a_ocn_to_grid
    public :: climber3a_jorge_to_grid

    ! For jorge testing 
    type g40_type  
        double precision, allocatable :: lon(:), lat(:)
        double precision, allocatable :: zs(:), tann(:), tsum(:), prec(:)
    end type 

contains 

    subroutine climber3a_atm_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
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
        double precision :: sigma, lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:), var(:,:)
            double precision, allocatable :: zs(:,:) 
            double precision :: lapse_ann    = 8.0d-3 
            double precision :: lapse_summer = 6.5d-3
        end type 

        type(inp_type)     :: inp
        type(grid_class)   :: grid0
        character(len=256) :: fldr_in, file_in_topo, file_in, nm_topo 
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, np 

        type(map_class)  :: map
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), outzs(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        ! For intermediate interpolation 
        character(len=256) :: c3_grid

        ! Define the input filenames
        fldr_in      = trim(path_in)
        file_in_topo = trim(fldr_in)//trim(domain)//"_horo.cdf"
        file_in      = trim(fldr_in)//trim(domain)//".cdf"

        desc    = "CLIMBER-3alpha simulation output"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                              trim(grid%name)//"_"//trim(domain)//".nc"

        ! Load the domain information 
        if (nc_exists_var(file_in_topo,"XT_I")) then 
            nx = nc_size(file_in_topo,"XT_I")
            ny = nc_size(file_in_topo,"YT_J")
            c3_grid = "climber3a-atmos"
        else
            nx = nc_size(file_in_topo,"lon")
            ny = nc_size(file_in_topo,"lat")
            c3_grid = "climber3a-atmos-hi"
        end if 
        np = nx*ny 

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))
        allocate(inp%zs(nx,ny))

        if (nc_exists_var(file_in_topo,"XT_I")) then 
            call nc_read(file_in_topo,"XT_I",inp%lon)
            call nc_read(file_in_topo,"YT_J",inp%lat)
        else 
            call nc_read(file_in_topo,"lon",inp%lon)
            call nc_read(file_in_topo,"lat",inp%lat)
        end if 
        
        ! Define CLIMBER3a points and input variable field
        call grid_init(grid0,name=c3_grid,mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat )

!         ! Define CLIMBER3a hi resolution intermediate interpolation grid
!         call grid_init(grid0hi,name="climber3a-5deg",mtype="latlon",units="degrees", &
!                          lon180=.TRUE.,x0=-180.d0,dx=5.d0,nx=73,y0=-90.d0,dy=5.d0,ny=37)

        ! Define the variables to be mapped 
        allocate(vars(3))
        call def_var_info(vars( 1),trim(file_in),"TS_ANN","t2m_ann",units="degrees Celcius", &
                          long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars( 2),trim(file_in),"TS_JJA","t2m_sum",units="degrees Celcius", &
                          long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars( 3),trim(file_in),"PRC_ANN","pr_ann",units="mm*d**-1", &
                          long_name="Precipitation, annual mean",method="nng")

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        call grid_allocate(grid,outzs) 
        
        ! Initialize mappings
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)
!         call map_init(map0hi,grid0,grid0hi,max_neighbors=10,lat_lim=8.d0,fldr="maps",load=.FALSE.)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Load reference topography in order to adjust temps to sea-level temps  
        call nc_read(file_in_topo,"HORO",inp%zs,missing_value=missing_value) 

        write(*,*) "input zs : ", minval(inp%zs), maxval(inp%zs)
        where(inp%zs .le. 1.d0) inp%zs = 0.d0  ! Clean-up climber ocean points 

        ! Map zs to new grid
        call map_field(map,"zs",inp%zs,outzs,outmask,"nng", &
                          fill=.TRUE.,missing_value=missing_value,sigma=sigma)
        call nc_write(filename,"zs",real(outzs),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,"zs","units","m")
        call nc_write_attr(filename,"zs","long_name","Surface elevation")
        call nc_write_attr(filename,"zs","coordinates","lat2D lon2D")
        
        ! Also generate a land mask
        inp%var = 0.d0
        where(inp%zs .gt. 1.d0) inp%var = 1.d0 

        ! Map mask to new grid
        call map_field(map,"mask_land",inp%var,outvar,outmask,"nn", &
                          fill=.TRUE.,missing_value=missing_value,sigma=sigma)
        call nc_write(filename,"mask_land",nint(outvar),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,"mask_land","units","1")
        call nc_write_attr(filename,"mask_land","long_name","Land mask (land=1)")
        call nc_write_attr(filename,"mask_land","coordinates","lat2D lon2D")

        ! ## Map climatological gridded variables ##
        
        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            ! Read in current variable
            call nc_read(trim(var_now%filename),var_now%nm_in,inp%var,missing_value=missing_value)
            where(abs(inp%var) .ge. 1d10) inp%var = missing_value 

            ! Scale to sea-level temperature for interpolation
            if (trim(var_now%nm_out) .eq. "t2m_sum") &
                inp%var = inp%var + inp%lapse_summer*inp%zs 
            if (trim(var_now%nm_out) .eq. "t2m_ann") &
                inp%var = inp%var + inp%lapse_ann*inp%zs 

            ! Map variable to new grid
            call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
                          fill=.TRUE.,missing_value=missing_value,sigma=sigma)
            
            ! Re-scale to near-surface temp for writing to file
            if (trim(var_now%nm_out) .eq. "t2m_sum") &
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

    subroutine climber3a_ocn_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
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

    end subroutine climber3a_ocn_to_grid


    subroutine climber3a_jorge_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CLIMBER3-alpha atmospheric DATA
        !       already bilinearly interpolated to the g40 grid by Jorge
        !       testing for comparison
        !
        ! =========================================================
        implicit none 

        character(len=*) :: domain, outfldr, subfldr, path_in 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: sigma, lat_lim 

        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        ! For jorge testing
        type(g40_type) :: g40 
        double precision :: lapse_ann    = 8.0d-3 
        double precision :: lapse_summer = 6.5d-3
        double precision, allocatable :: invar(:) 

        type(points_class)   :: pts0
        character(len=256) :: fldr_in 
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, np 

        type(map_class)  :: map
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), outzs(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        ! For intermediate interpolation 
        character(len=256) :: c3_grid

        ! Define the input filenames
        fldr_in      = trim(path_in)

        desc    = "CLIMBER-3alpha simulation output"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                              trim(grid%name)//"_"//trim(domain)//".nc"

        ! Load the domain information 
        np = 58081
        allocate(invar(np))
        c3_grid = "g40"
        call load_g40(g40,fldr_in,dataset=trim(domain))

        ! Define CLIMBER3a points and input variable field
        call points_init(pts0,name=c3_grid,mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=g40%lon,y=g40%lat )

        ! Define the variables to be mapped 
        allocate(vars(3))
        call def_var_info(vars( 1),"None","TS_ANN","t2m_ann",units="degrees Celcius", &
                          long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars( 2),"None","TS_JJA","t2m_sum",units="degrees Celcius", &
                          long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars( 3),"None","PRC_ANN","pr_ann",units="mm*d**-1", &
                          long_name="Precipitation, annual mean",method="nng")

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        call grid_allocate(grid,outzs) 
        
        ! Initialize mappings
        call map_init(map,pts0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        write(*,*) "input zs : ", minval(g40%zs), maxval(g40%zs)
        where(g40%zs .le. 1.d0) g40%zs = 0.d0  ! Clean-up climber ocean points 

        ! Map zs to new grid
        call map_field(map,"zs",g40%zs,outzs,outmask,"nn", &
                          fill=.TRUE.,missing_value=missing_value)
        call nc_write(filename,"zs",real(outzs),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,"zs","units","m")
        call nc_write_attr(filename,"zs","long_name","Surface elevation")
        call nc_write_attr(filename,"zs","coordinates","lat2D lon2D")
        
        ! Also generate a land mask
        invar = 0.d0
        where(g40%zs .gt. 1.d0) invar = 1.d0 

        ! Map mask to new grid
        call map_field(map,"mask_land",invar,outvar,outmask,"nn", &
                          fill=.TRUE.,missing_value=missing_value)
        call nc_write(filename,"mask_land",nint(outvar),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,"mask_land","units","1")
        call nc_write_attr(filename,"mask_land","long_name","Land mask (land=1)")
        call nc_write_attr(filename,"mask_land","coordinates","lat2D lon2D")

        ! ## Map climatological gridded variables ##
        
        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            select case(trim(var_now%nm_out))
                case("t2m_ann")
                    invar = g40%tann 
                case("t2m_sum")
                    invar = g40%tsum 
                case("pr_ann")
                    invar = g40%prec 
                case default 
                    write(*,*) "CASE NOT FOUND: "//trim(var_now%nm_out)
                    stop 

            end select 

            ! Scale to sea-level temperature for interpolation
            if (trim(var_now%nm_out) .eq. "t2m_sum") &
                invar = invar + lapse_summer*g40%zs 
            if (trim(var_now%nm_out) .eq. "t2m_ann") &
                invar = invar + lapse_ann*g40%zs 

            ! Map variable to new grid
            call map_field(map,var_now%nm_in,invar,outvar,outmask,"nn", &
                          fill=.TRUE.,missing_value=missing_value)
            
            ! Re-scale to near-surface temp for writing to file
            if (trim(var_now%nm_out) .eq. "t2m_sum") &
                outvar = outvar - lapse_summer*outzs 
            if (trim(var_now%nm_out) .eq. "t2m_ann") &
                outvar = outvar - lapse_ann*outzs 

            ! Write output variable to output file
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine climber3a_jorge_to_grid

    subroutine load_g40(g40,fldr,dataset)
        implicit none 

        type(g40_type)     :: g40 
        character(len=*)   :: fldr, dataset 
        character(len=512) :: filename 

        integer :: n = 58081   ! Number of points 

        allocate(g40%lon(n))
        allocate(g40%lat(n))
        allocate(g40%zs(n))
        allocate(g40%tann(n))
        allocate(g40%tsum(n))
        allocate(g40%prec(n))
        
        filename = trim(fldr)//"/coord-nord-40km.dat"
        g40%lon = read_vector(filename,n=n,col=5,skip=2)
        g40%lat = read_vector(filename,n=n,col=6,skip=2)

        filename = trim(fldr)//"/"//trim(dataset)//".g40"
        g40%zs    = read_vector(filename,n=n,col=4,skip=2)
        g40%tann  = read_vector(filename,n=n,col=2,skip=2)
        g40%tsum  = read_vector(filename,n=n,col=3,skip=2)
        g40%prec  = read_vector(filename,n=n,col=1,skip=2)
        
        write(*,*) "Loaded g40 dataset: "//trim(dataset)
        write(*,*) "lon: ", minval(g40%lon), maxval(g40%lon)
        write(*,*) "lat: ", minval(g40%lat), maxval(g40%lat)
        write(*,*) "zs: ",  minval(g40%zs), maxval(g40%zs)
        write(*,*) "tann: ",minval(g40%tann), maxval(g40%tann)
        write(*,*) "tsum: ",minval(g40%tsum), maxval(g40%tsum)
        write(*,*) "prec: ",minval(g40%prec), maxval(g40%prec)
        
        return 

    end subroutine load_g40 

    function read_vector(filename,n,col,skip) result(var)
        implicit none 

        character(len=*) :: filename 
        integer :: n, col, skip 
        real(8) :: var(n), tmp(50)
        character(len=10) :: tmpc
        integer :: i 

        open(16,file=trim(filename),status="old")
        do i = 1, skip
            read(16,*) tmpc 
        end do 

        do i = 1, n 
            read(16,*) tmp(1:col-1), var(i)
        end do 

        close(16)

        return
    end function read_vector

end module climber3a 

