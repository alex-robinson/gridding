module pmip3 

    use gridding_datasets
    use coord
    use ncio 
    use gaussian_filter
    
    implicit none 

    type pmip_info_type

        character(len=256) :: model 
        character(len=256) :: experiment 
        character(len=256) :: pmip_case 
        character(len=256) :: file_in_suffix
        character(len=256) :: desc
        character(len=256) :: grid_name 
        logical            :: is_pd 
        
        character(len=256) :: nm_tas_ann
        character(len=256) :: nm_tas_djf 
        character(len=256) :: nm_tas_jja
        character(len=256) :: nm_pr_ann
        character(len=256) :: nm_sst_ann 

    end type 

    private
    public :: pmip3_to_grid
    public :: LGM_extensions_to_grid

contains

    subroutine pmip3_info(info,model,experiment,domain)

        implicit none 

        type(pmip_info_type), intent(INOUT) :: info 
        character(len=*),     intent(IN)    :: model 
        character(len=*),     intent(IN)    :: experiment 
        character(len=*),     intent(IN)    :: domain 

        ! Local variables  
        character(len=512) :: str 

        ! === Default values for all cases =====

        info%model       = trim(model)
        info%experiment  = trim(experiment)

        info%pmip_case      = trim(model)//"-"//trim(experiment)
        !info%file_in_suffix = trim(model)//"/"//trim(model)//"_"//trim(experiment)//".nc"
        info%desc           = "PMIP3 "//trim(model)//" - "//trim(experiment)
        info%grid_name      = trim(model)//"-grid"

        if (trim(experiment) .eq. "lgm") then 
            ! Assume LGM topography 
            info%is_pd = .FALSE. 
        else
        ! Assume present-day topography 
            info%is_pd = .TRUE. 
        end if

        info%nm_tas_ann = "tas_spatialmean_ann" 
        info%nm_tas_djf = "tas_spatialmean_djf" 
        info%nm_tas_jja = "tas_spatialmean_jja" 
        info%nm_pr_ann  = "pr_spatialmean_ann" 
        info%nm_sst_ann = "sst_spatialmean_ann" 
        
        select case(trim(info%pmip_case))
            ! 8 models with all three experiments of interest (piControl,midHolocene,lgm)

            case("CCSM4-piControl")
                str = "CCSM4_0_piControl.cvdp_data.250-1300.nc"
            case("CCSM4-midHolocene")
                str = "CCSM4_0_midHolocene.cvdp_data.1000-1300.nc"
            case("CCSM4-lgm")
                str = "CCSM4_0_lgm.cvdp_data.1800-1900.nc"
            
            case("CNRM-CM5-piControl")
                str = "CNRM-CM5_0_piControl.cvdp_data.2400-2699.nc"
            case("CNRM-CM5-midHolocene")
                str = "CNRM-CM5_0_midHolocene.cvdp_data.1950-2149.nc"
            case("CNRM-CM5-lgm")
                str = "CNRM-CM5_0_lgm.cvdp_data.1800-1999.nc"

            case("FGOALS-g2-piControl")
                str = "FGOALS-g2_0_piControl.cvdp_data.701-900.nc"
            case("FGOALS-g2-midHolocene")
                str = "FGOALS-g2_0_midHolocene.cvdp_data.801-1000.nc"
            case("FGOALS-g2-lgm")
                str = "FGOALS-g2_0_lgm.cvdp_data.550-649.nc"
            
            case("GISS-E2-R-piControl")
                str = "GISS-E2-R_0_piControl.cvdp_data.3981-4530.nc"
            case("GISS-E2-R-midHolocene")
                str = "GISS-E2-R_0_midHolocene.cvdp_data.2500-2599.nc"
            case("GISS-E2-R-lgm")
                str = "GISS-E2-R_0_lgm.cvdp_data.3000-3099.nc"

            case("IPSL-CM5A-LR-piControl")
                str = "IPSL-CM5A-LR_0_piControl.cvdp_data.1800-2799.nc"
            case("IPSL-CM5A-LR-midHolocene")
                str = "IPSL-CM5A-LR_0_midHolocene.cvdp_data.2301-2800.nc"
            case("IPSL-CM5A-LR-lgm")
                str = "IPSL-CM5A-LR_0_lgm.cvdp_data.2601-2800.nc"

            case("MIROC-ESM-piControl")
                str = "MIROC-ESM_0_piControl.cvdp_data.1800-2429.nc"
            case("MIROC-ESM-midHolocene")
                str = "MIROC-ESM_0_midHolocene.cvdp_data.2330-2429.nc"
            case("MIROC-ESM-lgm")
                str = "MIROC-ESM_0_lgm.cvdp_data.4600-4699.nc"

            case("MPI-ESM-P-piControl")
                str = "MPI-ESM-P_0_piControl.cvdp_data.2400-3000.nc"
            case("MPI-ESM-P-midHolocene")
                str = "MPI-ESM-P_0_midHolocene.cvdp_data.1850-1949.nc"
            case("MPI-ESM-P-lgm")
                str = "MPI-ESM-P_0_lgm.cvdp_data.1850-1949.nc"
            
            case("MRI-CGCM3-piControl")
                str = "MRI-CGCM3_0_piControl.cvdp_data.2001-2350.nc"
            case("MRI-CGCM3-midHolocene")
                str = "MRI-CGCM3_0_midHolocene.cvdp_data.1951-2050.nc"
            case("MRI-CGCM3-lgm")
                str = "MRI-CGCM3_0_lgm.cvdp_data.2501-2600.nc"
                      
            case DEFAULT 

                write(*,*) "pmip3_info:: Error: case not recognized: "//trim(info%pmip_case)
                stop 

        end select 

        info%file_in_suffix = trim(str) 


        ! Write summary 
        write(*,*) "PMIP3-info: ", trim(info%pmip_case), info%is_pd, trim(info%file_in_suffix)

        return 

    end subroutine pmip3_info


    subroutine pmip3_to_grid(outfldr,grid,domain,model,experiment,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CCSM4 ATM LGM data
        !
        ! =========================================================
        implicit none
        
        character(len=*), intent(IN) :: outfldr 
        type(grid_class), intent(IN) :: grid
        character(len=*), intent(IN) :: domain
        character(len=*), intent(IN) :: model
        character(len=*), intent(IN) :: experiment
        character(len=*), intent(IN) :: path_in
        double precision, intent(IN) :: sigma
        integer,          intent(IN) :: max_neighbors
        double precision, intent(IN) :: lat_lim

        ! Local variables 
        
        type(pmip_info_type) :: info

        character(len=512) :: filename
        character(len=1024) :: desc, ref
        
        type inp_type
            double precision, allocatable :: lon(:), lat(:), var(:,:)
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

        type(map_scrip_class) :: mps

        ! Get PMIP info 
        call pmip3_info(info,model,experiment,domain)

        ! Define the input filenames
        file_in          = trim(path_in)//trim(info%file_in_suffix)

        desc    = trim(info%desc)
        ref     = "source folder: "//trim(path_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_"//trim(info%pmip_case)//".nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name=trim(info%grid_name),mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        allocate(vars(5))
        call def_var_info(vars(1),trim(file_in),trim(info%nm_tas_ann),"t2m_ann",units="degC",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),trim(info%nm_tas_djf),"t2m_djf",units="degC",long_name="Near-surface temperature (2-m), DJF mean",method="nng")
        call def_var_info(vars(3),trim(file_in),trim(info%nm_tas_jja),"t2m_jja",units="degC",long_name="Near-surface temperature (2-m), JJA mean",method="nng")
        call def_var_info(vars(4),trim(file_in),trim(info%nm_pr_ann), "pr_ann", units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
        call def_var_info(vars(5),trim(file_in),trim(info%nm_sst_ann),"sst_ann",units="degC",long_name="Sea-surface temperature, annual mean",method="nng") 
        
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        ! call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Assume SCRIP map is already written.
        ! (done manually to avoid generating a huge grid object via 
        ! cdo griddes path_to_grid_file.nc > maps/grid_GRIDNAME.txt
        !call grid_write_cdo_desc_short(grid0,fldr="maps") 
        !call grid_write_cdo_desc_cdo(grid0%name,fldr="maps",file_nc=file_in)
        call grid_write_cdo_desc_explicit_latlon(real(grid0%G%x,4),real(grid0%G%y,4),grid0%name,fldr="maps")

        ! Define output grid in grid description file 
        call grid_write_cdo_desc_short(grid,fldr="maps") 
        
        ! Generate SCRIP interpolation weights 
        call map_scrip_init(mps,grid0%name,grid%name,fldr="maps",src_nc=file_in)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="km")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="km")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)
    
        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            ! Read in current variable
            call nc_read(trim(var_now%filename),var_now%nm_in,inp%var,missing_value=mv)
            
            where((abs(inp%var) .ge. 1d10)) inp%var = mv

            ! Map variable to new grid
            ! call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
            !               fill=.TRUE.,missing_value=mv,sigma=sigma)
            call map_scrip_field(mps,var_now%nm_in,inp%var,outvar,method="mean",missing_value=mv)

            outvar = outvar
            ! Write output variable to output file
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        end do

        ! Add topography 
        call pmip3_add_zs_to_grid(filename,outfldr,grid,domain,path_in,sigma, &
                                                max_neighbors,lat_lim,is_pd=info%is_pd)

        return 

    end subroutine pmip3_to_grid

    subroutine pmip3_add_zs_to_grid(filename,outfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim,is_pd)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       LGM dzs wrt PD is given for PMIP3 (Abe-Ouchi 2015)
        !       If present day, z_srf = rtopo
        !       If LGM,         z_srf = rtopo+dzs_pmip3
        !
        ! =========================================================
        implicit none
        
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: outfldr
        type(grid_class), intent(IN) :: grid
        character(len=*), intent(IN) :: domain, path_in
        integer,          intent(IN) :: max_neighbors
        double precision, intent(IN) :: sigma, lat_lim
        logical,          intent(IN) :: is_pd           ! Is this file present day? 

        ! Local variables 

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
        double precision, allocatable :: outvar(:,:), outvar1(:,:)
        integer, allocatable          :: outmask(:,:)

        double precision, allocatable :: pd_zs(:,:)

        integer :: q, k, m, i, l, n_var

        type(map_scrip_class) :: mps

        ! For intermediate interpolation 
        character(len=256) :: pmip3_grid

        ! Define the input filenames
        fldr_in          = trim(path_in)//"/extensions/"
        file_in          = trim(fldr_in)//"pmip3_21k_v0.nc"

        select case(trim(domain))

            case("Greenland")

            ! Filename of already-processed present-day topography to combine with dzs here 
            file_in_topo = trim(outfldr)//"/../"//trim(grid%name)//"_TOPO-M17.nc"

            case DEFAULT ! Antarctica, North, Laurentide, etc...

            ! Filename of already-processed present-day topography to combine with dzs here 
            file_in_topo = trim(outfldr)//"/../"//trim(grid%name)//"_TOPO-RTOPO-2.0.1.nc"

        end select 

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="pmip3_dzs",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat) 

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outvar1)
        call grid_allocate(grid,outmask)
        call grid_allocate(grid,pd_zs)

        ! Load the present-day surface elevation field (do it now, to make sure it's available)
        call nc_read(file_in_topo,"z_srf",pd_zs)
        
        ! Smooth it out to match target smoothness via sigma 
        call filter_gaussian(var=pd_zs,sigma=sigma,dx=grid%G%dx)

        if (is_pd) then 
            ! Simply writing present-day surface elevation to file 

            ! Write output variable to output file
            call nc_write(filename,"z_srf",real(pd_zs),dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,"z_srf","units","m")
            call nc_write_attr(filename,"z_srf","long_name","Surface elevation")
            call nc_write_attr(filename,"z_srf","coordinates","lat2D lon2D")

            ! ! Load PD mask
            ! call nc_read(file_in_topo,"mask",outmask)

            ! ! Write output variable to output file
            ! call nc_write(filename,"mask",outmask,dim1="xc",dim2="yc")

            ! ! Write variable metadata
            ! call nc_write_attr(filename,"mask","units","")
            ! call nc_write_attr(filename,"mask","long_name","Mask (0:ocean, 1:land, 2:grounded ice, 3:floating ice)")
            ! call nc_write_attr(filename,"mask","coordinates","lat2D lon2D")

        else 

            ! Initialize mapping
            ! call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

            ! Assume SCRIP map is already written.
            ! (done manually to avoid generating a huge grid object via 
            ! cdo griddes path_to_grid_file.nc > maps/grid_GRIDNAME.txt
            !call grid_write_cdo_desc_short(grid0,fldr="maps") 
            
            ! Define output grid in grid description file 
            call grid_write_cdo_desc_short(grid,fldr="maps") 
            
            ! Generate SCRIP interpolation weights 
            call map_scrip_init(mps,grid0%name,grid%name,fldr="maps",src_nc=file_in)


            ! dz_srf =========

            ! Read in current variable
            call nc_read(trim(file_in),"orog_diff",inp%var,missing_value=mv)
            
            where((abs(inp%var) .ge. 1d10)) inp%var = mv

            ! Map variable to new grid
            ! call map_field(map,"dz_srf",inp%var,outvar,outmask,method="nng", &
            !               fill=.TRUE.,missing_value=mv,sigma=sigma)

            call map_scrip_field(mps,"dz_srf",inp%var,outvar,method="mean",missing_value=mv)

            ! Write output variable to output file
            call nc_write(filename,"dz_srf",real(outvar),dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,"dz_srf","units","m")
            call nc_write_attr(filename,"dz_srf","long_name","Surface elevation difference")
            call nc_write_attr(filename,"dz_srf","coordinates","lat2D lon2D")

            ! Next, add present-day surface elevation to field to get z_srf 
            outvar = outvar + pd_zs 

            ! Write output variable to output file
            call nc_write(filename,"z_srf",real(outvar),dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,"z_srf","units","m")
            call nc_write_attr(filename,"z_srf","long_name","Surface elevation")
            call nc_write_attr(filename,"z_srf","coordinates","lat2D lon2D")

            ! mask =========

            ! Read in current variable - land mask
            call nc_read(trim(file_in),"mask3",inp%var,missing_value=mv)
            
            where((abs(inp%var) .ge. 1d10)) inp%var = mv

            ! Map variable to new grid
            call map_field(map,"mask3",inp%var,outvar,outmask,method="nn", &
                          fill=.TRUE.,missing_value=mv,sigma=sigma)


            ! Read in current variable - ice mask 
            call nc_read(trim(file_in),"mask1",inp%var,missing_value=mv)
            
            where((abs(inp%var) .ge. 1d10)) inp%var = mv

            ! Map variable to new grid
            call map_field(map,"mask1",inp%var,outvar1,outmask,method="nn", &
                          fill=.TRUE.,missing_value=mv,sigma=sigma)
            
            ! Define final ice-land-mask 
            outmask = 0 
            where(outvar .gt. 0.5)  outmask = 1     ! Land 
            where(outvar1 .gt. 0.5) outmask = 2     ! Land-ice 
            where(outvar1 .gt. 0.5 .and. outvar .le. 0.5) outmask = 3  ! Floating ice 

            ! Write output variable to output file
            call nc_write(filename,"mask",outmask,dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,"mask","units","")
            call nc_write_attr(filename,"mask","long_name","Mask (0:ocean, 1:land, 2:grounded ice, 3:floating ice)")
            call nc_write_attr(filename,"mask","coordinates","lat2D lon2D")

        
        end if 


        return 

    end subroutine pmip3_add_zs_to_grid

    subroutine LGM_extensions_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       LGM extensions (Abe-Ouchi 2015)
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
        character(len=256) :: pmip3_grid

        ! Define the input filenames
        fldr_in          = trim(path_in)//"/extensions/"
        file_in          = trim(fldr_in)//"pmip3_21k_extensions_ant.nc"

        desc    = "LGM extensions"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_PMIP3_21K_extensions.nc"

        !nx = nc_size(file_in,"lon")
        !ny = nc_size(file_in,"lat")

        !allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        !call nc_read(file_in,"lon",inp%lon)
        !call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="ext_LGM-grid",mtype="polar_stereographic",units="km", &
                         lon180=.TRUE.,x0=-3040.d0,dx=32.d0,nx=191,y0=3040.d0,dy=-32.d0,ny=191, &
                   lambda=0.d0,phi=-90.d0,alpha=19.0d0) 

        ! jablasco
        allocate(vars(2))
        call def_var_info(vars(1),trim(file_in),"mask_anu","mask_anu",units="%",long_name="Mask grounded ice fraction",method="nng")
        call def_var_info(vars(2),trim(file_in),"mask_ice6g","mask_ice6g",units="%",long_name="Mask grounded ice",method="nng")
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="km")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="km")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)
    
        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            ! Read in current variable
            call nc_read(trim(var_now%filename),var_now%nm_in,inp%var,missing_value=mv)
            
            where((abs(inp%var) .ge. 1d10)) inp%var = mv

            ! Map variable to new grid
            call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
                          fill=.TRUE.,missing_value=mv,sigma=sigma)

            outvar = outvar
            ! Write output variable to output file
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        end do

        return 

    end subroutine LGM_extensions_to_grid

end module pmip3
