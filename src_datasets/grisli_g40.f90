module grisli_g40

    use gridding_datasets
    use coordinates 
    use interp2D 
    use gaussian_filter
    use ncio 
    
    implicit none 

    private 
    public :: g40_topo_to_grid
    public :: g40_climber3a_to_grid

    ! For jorge testing 
    type g40_type  
        double precision, allocatable :: lon(:), lat(:)
        double precision, allocatable :: zs(:), tann(:), tsum(:), prec(:)
    end type 


contains

    subroutine g40_topo_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       TOPOGRAPHIC DATA
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
        call load_g40(g40,fldr_in,dataset=trim(domain),load_topo=.FALSE.)

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

    end subroutine g40_topo_to_grid

    subroutine g40_climber3a_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim, &
                                       load_topo)
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
        logical :: load_topo 
        
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
        call load_g40(g40,fldr_in,dataset=trim(domain),load_topo=load_topo)

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

    end subroutine g40_climber3a_to_grid


    subroutine load_g40(g40,fldr,dataset,load_topo)
        implicit none 

        type(g40_type)     :: g40 
        character(len=*)   :: fldr, dataset 
        character(len=512) :: filename 
        logical :: load_topo 

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
        g40%prec  = read_vector(filename,n=n,col=1,skip=2)
        g40%tann  = read_vector(filename,n=n,col=2,skip=2)
        g40%tsum  = read_vector(filename,n=n,col=3,skip=2)
        
        g40%zs = 0.d0 
        if (load_topo) g40%zs    = read_vector(filename,n=n,col=4,skip=2)
        
        write(*,*) "Loaded g40 dataset: "//trim(dataset)
        write(*,*) "lon: ", minval(g40%lon), maxval(g40%lon)
        write(*,*) "lat: ", minval(g40%lat), maxval(g40%lat)
        write(*,*) "zs: ",  minval(g40%zs), maxval(g40%zs)
        write(*,*) "tann: ",minval(g40%tann), maxval(g40%tann)
        write(*,*) "tsum: ",minval(g40%tsum), maxval(g40%tsum)
        write(*,*) "prec: ",minval(g40%prec), maxval(g40%prec)
        
        return 

    end subroutine load_g40 

end module grisli_g40


