module climberx
    
    use gridding_datasets
    use coord
    use ncio 
    use gaussian_filter 

    implicit none 

    private 
    public :: climberx_to_grid 
    
contains

    subroutine climberx_to_grid(outfldr,grid,expname,infldr,sigma)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !    CLIMBER-X data files (default information of interest)
        !
        ! =========================================================

        implicit none 

        character(len=*), intent(IN) :: outfldr
        type(grid_class), intent(IN) :: grid 
        character(len=*), intent(IN) :: expname
        character(len=*), intent(IN) :: infldr
        double precision, intent(IN) :: sigma

        ! Local variables 
        integer             :: max_neighbors 
        double precision    :: lat_lim 
        character(len=512)  :: filename 
        character(len=1024) :: desc, ref 

        type(grid_class)    :: grid0 
        character(len=1024) :: filename0  

        type input_type 
            integer :: nx, ny, nz 
            double precision, allocatable :: lon(:) 
            double precision, allocatable :: lat(:) 
            double precision, allocatable :: lev(:) 
            double precision, allocatable :: var2D(:,:) 
            double precision, allocatable :: var3D(:,:,:)
            double precision, allocatable :: var4D(:,:,:,:)
        end type 

        type(input_type) :: inp

        type(map_scrip_class) :: mps 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer, allocatable :: dims(:) 
        character(len=56), allocatable :: dim_names(:) 

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                          "_CLIMBERX-"//trim(expname)//".nc"

        ! ==== Input grid and dataset information ==== 

        filename0 = "maps/climberx_grid_5x5.nc"
        inp%nx = nc_size(filename0,"lon") 
        inp%ny = nc_size(filename0,"lat") 
        allocate(inp%lon(inp%nx))
        allocate(inp%lat(inp%ny))
        allocate(inp%var2D(inp%nx,inp%ny))
        allocate(inp%var3D(inp%nx,inp%ny,13))

        call nc_read(filename0,"lon",inp%lon)
        call nc_read(filename0,"lat",inp%lat)
        
        ! Also get ocean depth (lev)
        filename0 = "maps/climberx_grid-ocn_5x5.nc"
        inp%nz = nc_size(filename0,"lev") 
        allocate(inp%lev(inp%nz))
        call nc_read(filename0,"lev",inp%lev) 
        allocate(inp%var4D(inp%nx,inp%ny,inp%nz,13))

        call grid_init(grid0,name="climberx",mtype="latlon",units="degrees", &
                                            lon180=.TRUE.,x=inp%lon,y=inp%lat)

        ! Define the input filenames
        desc    = "CLIMBER-X simulation results. Experiment = "//trim(expname)
        ref     = "Willeit et al. (2021, gmd)"

        
        ! ==== MAPPING INFORMATION ====

        ! Define input grid in grid description file
        ! Note: this grid definition file was prepared by hand
        ! as climber-x output is not yet prepared to define 
        ! grid properly for cdo.
        !call grid_write_cdo_desc_short(grid0,fldr="maps") 
        
        ! Define output grid in grid description file 
        call grid_write_cdo_desc_short(grid,fldr="maps") 
        
        ! Generate SCRIP interpolation weights 
        call map_scrip_init(mps,grid0%name,grid%name,fldr="maps",src_nc=filename0)


        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"lev",  x=inp%lev,units="m")
        call nc_write_dim(filename,"month",x=1,dx=1,nx=12,units="month")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ============================================================
            
        filename0 = trim(infldr)//"/geo.nc"

        ! == Surface elevation == 

        ! Read input data 
        call nc_read(filename0,"z_sur",inp%var2D,missing_value=mv)

        ! Perform mapping and write to file 
        call map_climberx_field_2D(mps,inp%var2D,outvar,outmask,filename,"z_srf","m", &
                        "Surface elevation",dx=grid%G%dx,sigma=sigma)

        ! == Ice thickness == 

        ! Read input data 
        call nc_read(filename0,"z_ice",inp%var2D,missing_value=mv)

        ! Perform mapping and write to file 
        call map_climberx_field_2D(mps,inp%var2D,outvar,outmask,filename,"H_ice","m", &
                        "Ice thickness",dx=grid%G%dx,sigma=sigma)

        ! == Land fraction == 

        ! Read input data 
        call nc_read(filename0,"f_lnd",inp%var2D,missing_value=mv)

        ! Perform mapping and write to file 
        call map_climberx_field_2D(mps,inp%var2D,outvar,outmask,filename,"f_lnd","1", &
                        "Land fraction",dx=grid%G%dx,sigma=sigma)
        
        ! == Ice fraction == 
        
        ! Read input data 
        call nc_read(filename0,"f_ice",inp%var2D,missing_value=mv)

        ! Perform mapping and write to file 
        call map_climberx_field_2D(mps,inp%var2D,outvar,outmask,filename,"f_ice","1", &
                        "Land fraction",dx=grid%G%dx,sigma=sigma)
        
        ! ============================================================
        
        filename0 = trim(infldr)//"/atm.nc"
        
        ! == temperature == 

        ! Read input data 
        call nc_dims(filename0,"t2a",dim_names,dims)
        call nc_read(filename0,"t2a",inp%var3D,missing_value=mv, &
                start=[1,1,1,dims(4)],count=[dims(1),dims(2),13,1])

        ! Perform mapping and write to file 
        call map_climberx_field_3D(mps,inp%var3D,outvar,outmask,filename,"t2m","K", &
                        "Near-surface air temperature",dx=grid%G%dx,sigma=sigma)

        ! ann mean 
        call map_climberx_field_2D(mps,inp%var3D(:,:,13),outvar,outmask,filename,"t2m_ann","K", &
                        "Near-surface air temperature (ANN)",dx=grid%G%dx,sigma=sigma)

        ! djf mean 
        inp%var2D = sum(inp%var3D(:,:,[12,1,2]),dim=3) / 3.0d0 
        call map_climberx_field_2D(mps,inp%var2D,outvar,outmask,filename,"t2m_djf","K", &
                        "Near-surface air temperature (DJF)",dx=grid%G%dx,sigma=sigma)

        ! jja mean 
        inp%var2D = sum(inp%var3D(:,:,6:8),dim=3) / 3.0d0 
        call map_climberx_field_2D(mps,inp%var2D,outvar,outmask,filename,"t2m_jja","K", &
                        "Near-surface air temperature (JJA)",dx=grid%G%dx,sigma=sigma)

        ! == precipitation == 

        ! Read input data 
        call nc_dims(filename0,"prc",dim_names,dims)
        call nc_read(filename0,"prc",inp%var3D,missing_value=mv, &
                start=[1,1,1,dims(4)],count=[dims(1),dims(2),13,1])

        ! Perform mapping and write to file 
        call map_climberx_field_3D(mps,inp%var3D,outvar,outmask,filename,"pr","kg/m2/day", &
                        "Precipitation",dx=grid%G%dx,sigma=sigma)

        ! ann mean 
        inp%var2D = inp%var3D(:,:,13)*360d0  ! [kg/m2/day] => [kg/m2/yr]
        call map_climberx_field_2D(mps,inp%var2D,outvar,outmask,filename,"pr_ann","kg/m2/yr", &
                        "Precipitation (ANN)",dx=grid%G%dx,sigma=sigma)

        ! == snowfall == 

        ! Read input data 
        call nc_dims(filename0,"prcs",dim_names,dims)
        call nc_read(filename0,"prcs",inp%var3D,missing_value=mv, &
                start=[1,1,1,dims(4)],count=[dims(1),dims(2),13,1])

        ! Perform mapping and write to file 
        call map_climberx_field_3D(mps,inp%var3D,outvar,outmask,filename,"sf","kg/m2/day", &
                        "Snowfall",dx=grid%G%dx,sigma=sigma)

        ! ann mean 
        inp%var2D = inp%var3D(:,:,13)*360d0  ! [kg/m2/day] => [kg/m2/yr]
        call map_climberx_field_2D(mps,inp%var2D,outvar,outmask,filename,"sf_ann","kg/m2/yr", &
                        "Snowfall (ANN)",dx=grid%G%dx,sigma=sigma)

        ! == sea-level pressure == 

        ! Read input data 
        call nc_dims(filename0,"slp",dim_names,dims)
        call nc_read(filename0,"slp",inp%var3D,missing_value=mv, &
                start=[1,1,1,dims(4)],count=[dims(1),dims(2),13,1])

        ! Perform mapping and write to file 
        call map_climberx_field_3D(mps,inp%var3D,outvar,outmask,filename,"slp","Pa", &
                        "Sea-level pressure",dx=grid%G%dx,sigma=sigma)

        ! ann mean 
        inp%var2D = inp%var3D(:,:,13)
        call map_climberx_field_2D(mps,inp%var2D,outvar,outmask,filename,"slp_ann","Pa", &
                        "Sea-level pressure (ANN)",dx=grid%G%dx,sigma=sigma)

        ! ! == sea-level temperature == 

        ! ! Read input data 
        ! call nc_dims(filename0,"tsl",dim_names,dims)
        ! call nc_read(filename0,"tsl",inp%var3D,missing_value=mv, &
        !         start=[1,1,1,dims(4)],count=[dims(1),dims(2),13,1])

        ! ! Perform mapping and write to file 
        ! call map_climberx_field_3D(mps,inp%var3D,outvar,outmask,filename,"tsl","K", &
        !                 "Sea-level temperature",dx=grid%G%dx,sigma=sigma)


        return 

    end subroutine climberx_to_grid

    subroutine map_climberx_field_2D(mps,invar,outvar,outmask,filename,nm,units,long_name,dx,sigma)

        implicit none 

        type(map_scrip_class), intent(IN) :: mps
        double precision,  intent(IN)    :: invar(:,:) 
        double precision,  intent(INOUT) :: outvar(:,:) 
        integer,           intent(OUT)   :: outmask(:,:) 
        character(len=*),  intent(IN)    :: filename 
        character(len=*),  intent(IN)    :: nm 
        character(len=*),  intent(IN)    :: units 
        character(len=*),  intent(IN)    :: long_name  
        double precision,  intent(IN)    :: dx 
        double precision,  intent(IN)    :: sigma 

        ! Perform conservative interpolation 
        outvar = mv 
        call map_scrip_field(mps,"z_srf",invar,outvar,method="mean",missing_value=mv)
        
        ! Apply Gaussian smoothing to the field 
        outmask = 1
        call filter_gaussian(var=outvar,sigma=sigma,dx=dx,mask=outmask.eq.1)

        ! Write output data 
        call nc_write(filename,nm,real(outvar),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,nm,"units",units)
        call nc_write_attr(filename,nm,"long_name",long_name)
        call nc_write_attr(filename,nm,"coordinates","lat2D lon2D")
        
        return

    end subroutine map_climberx_field_2D

    subroutine map_climberx_field_3D(mps,invar,outvar,outmask,filename,nm,units,long_name,dx,sigma)

        implicit none 

        type(map_scrip_class), intent(IN) :: mps
        double precision,  intent(IN)    :: invar(:,:,:) 
        double precision,  intent(INOUT) :: outvar(:,:) 
        integer,           intent(OUT)   :: outmask(:,:) 
        character(len=*),  intent(IN)    :: filename 
        character(len=*),  intent(IN)    :: nm 
        character(len=*),  intent(IN)    :: units 
        character(len=*),  intent(IN)    :: long_name  
        double precision,  intent(IN)    :: dx 
        double precision,  intent(IN)    :: sigma 

        ! Local variables 
        integer :: k, nk 

        nk = 12 

        do k = 1, nk 

            ! Perform conservative interpolation 
            outvar = mv 
            call map_scrip_field(mps,"z_srf",invar(:,:,k),outvar,method="mean",missing_value=mv)
            
            ! Apply Gaussian smoothing to the field 
            outmask = 1
            call filter_gaussian(var=outvar,sigma=sigma,dx=dx,mask=outmask.eq.1)

            ! Write output data 
            call nc_write(filename,nm,real(outvar),dim1="xc",dim2="yc",dim3="month", &
                            start=[1,1,k],count=[size(outvar,1),size(outvar,2),1])

        end do 

        ! Write variable metadata
        call nc_write_attr(filename,nm,"units",units)
        call nc_write_attr(filename,nm,"long_name",long_name)
        call nc_write_attr(filename,nm,"coordinates","lat2D lon2D")
        
        return

    end subroutine map_climberx_field_3D
    
end module climberx