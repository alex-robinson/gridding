module pmip3 

    use gridding_datasets
    use coord
    use ncio 
    use gaussian_filter
    
    implicit none 

    type pmip_info_type

        character(len=256) :: nm_tas_ann
        character(len=256) :: nm_tas_sum
        character(len=256) :: nm_pr_ann
        character(len=256) :: file_suffix
        character(len=256) :: desc
        character(len=256) :: grid_name 

        logical            :: is_pd 
        
    end type 

    private
    public :: PMIP3_to_grid

!     public :: CCSM4_PD_to_grid
!     public :: CCSM4_LGM_to_grid
    public :: CNRM_CM5_PD_to_grid
    public :: CNRM_CM5_LGM_to_grid 
    public :: COSMOS_ASO_PD_to_grid
    public :: COSMOS_ASO_LGM_to_grid 
    public :: FGOALS_g2_PD_to_grid
    public :: FGOALS_g2_LGM_to_grid  
    public :: GISS_E2_R_150_PD_to_grid
    public :: GISS_E2_R_150_LGM_to_grid
    public :: GISS_E2_R_151_PD_to_grid
    public :: GISS_E2_R_151_LGM_to_grid
    public :: IPSL_CM5A_LR_PD_to_grid
    public :: IPSL_CM5A_LR_LGM_to_grid 
    public :: MIROC_ESM_PD_to_grid
    public :: MIROC_ESM_LGM_to_grid 
    public :: MPI_ESM_P_p1_PD_to_grid
    public :: MPI_ESM_P_p1_LGM_to_grid 
    public :: MPI_ESM_P_p2_PD_to_grid
    public :: MPI_ESM_P_p2_LGM_to_grid
    public :: MRI_CGCM3_PD_to_grid
    public :: MRI_CGCM3_LGM_to_grid 

    public :: LGM_dzs_to_grid
    public :: LGM_extensions_to_grid

contains

    subroutine pmip3_info(info,pmip_case,domain)

        implicit none 

        type(pmip_info_type), intent(INOUT) :: info 
        character(len=*),     intent(IN)    :: pmip_case 
        character(len=*),     intent(IN)    :: domain 

        ! Local variables 
        logical :: is_south 

        is_south = .FALSE. 
        if (trim(domain) .eq. "Antarctica") is_south = .TRUE. 

        ! === Default values for all cases =====

        info%nm_tas_ann = "tas_spatialmean_ann" 
        info%nm_pr_ann  = "pr_spatialmean_ann" 

        if (is_south) then 
            info%nm_tas_sum = "tas_spatialmean_djf" 
        else 
            info%nm_tas_sum = "tas_spatialmean_jja" 
        end if 

        ! === Model-experiment specific values =====

        select case(trim(pmip_case))

            case("CCSM4-piControl") 

                info%file_suffix = "CCSM4/CCSM4_piControl.nc"
                info%desc        = "CCSM4 PMIP3 ATM - piControl"
                info%grid_name   = "CCSM4-grid"
                info%is_pd       = .TRUE. 

            case("CCSM4-LGM") 

                info%file_suffix = "CCSM4/CCSM4_lgm.nc"
                info%desc        = "CCSM4 PMIP3 ATM - LGM"
                info%grid_name   = "CCSM4-grid"
                info%is_pd       = .FALSE. 

            case DEFAULT 

                write(*,*) "pmip3_info:: Error: case not recognized: "//trim(pmip_case)
                stop 

        end select 

        return 

    end subroutine pmip3_info 


    subroutine PMIP3_to_grid(outfldr,grid,domain,pmip_case,path_in,sigma,max_neighbors,lat_lim)
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
        character(len=*), intent(IN) :: pmip_case
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

        ! Get PMIP info 
        call pmip3_info(info,pmip_case,domain)

        ! Define the input filenames
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//trim(info%file_suffix)

        desc    = trim(info%desc)
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_"//trim(pmip_case)//".nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name=trim(info%grid_name),mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),trim(info%nm_tas_ann),"t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),trim(info%nm_tas_sum),"t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),trim(info%nm_pr_ann), "pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

        ! Add topography 
        call pmip3_add_zs_to_grid(filename,outfldr,grid,domain,path_in,sigma, &
                                                max_neighbors,lat_lim,is_pd=info%is_pd)

        return 

    end subroutine PMIP3_to_grid

    subroutine CCSM4_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CCSM4 ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"CCSM4/CCSM4_piControl.nc"

        desc    = "CCSM4 PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_CCSM4-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="CCSM4-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

        ! Add topography 
        call pmip3_add_zs_to_grid(filename,outfldr,grid,domain,path_in,sigma, &
                                                max_neighbors,lat_lim,is_pd=.TRUE.)

        return 

    end subroutine CCSM4_PD_to_grid

    subroutine CNRM_CM5_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CNRM_CM5 ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"CNRM_CM5/CNRM_CM5_lgm.nc"

        desc    = "CNRM_CM5 PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_CNRM_CM5-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="CNRM_CM5-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine CNRM_CM5_LGM_to_grid

    subroutine CNRM_CM5_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CNRM_CM5 ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"CNRM_CM5/CNRM_CM5_piControl.nc"

        desc    = "CNRM_CM5 PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_CNRM_CM5-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="CNRM_CM5-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine CNRM_CM5_PD_to_grid

    subroutine COSMOS_ASO_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       COSMOS_ASO ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"COSMOS_ASO/COSMOS_ASO_lgm.nc"

        desc    = "COSMOS_ASO PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_COSMOS_ASO-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="COSMOS_ASO-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine COSMOS_ASO_LGM_to_grid

    subroutine COSMOS_ASO_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       COSMOS_ASO ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"COSMOS_ASO/COSMOS_ASO_piControl.nc"

        desc    = "COSMOS_ASO PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_COSMOS_ASO-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="COSMOS_ASO-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine COSMOS_ASO_PD_to_grid

    subroutine FGOALS_g2_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       FGOALS_g2 ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"FGOALS_g2/FGOALS_g2_lgm.nc"

        desc    = "FGOALS_g2 PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_FGOALS_g2-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="FGOALS_g2-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine FGOALS_g2_LGM_to_grid

    subroutine FGOALS_g2_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       FGOALS_g2 ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"FGOALS_g2/FGOALS_g2_piControl.nc"

        desc    = "FGOALS_g2 PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_FGOALS_g2-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="FGOALS_g2-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine FGOALS_g2_PD_to_grid

    subroutine GISS_E2_R_150_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GISS_E2_R_150 ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"GISS_E2_R_150/GISS_E2_R_150_lgm.nc"

        desc    = "GISS_E2_R_150 PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_GISS_E2_R_150-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="GISS_E2_R_150-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine GISS_E2_R_150_LGM_to_grid

    subroutine GISS_E2_R_150_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GISS_E2_R_150 ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"GISS_E2_R_150/GISS_E2_R_150_piControl.nc"

        desc    = "GISS_E2_R_150 PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_GISS_E2_R_150-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="GISS_E2_R_150-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine GISS_E2_R_150_PD_to_grid

    subroutine GISS_E2_R_151_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GISS_E2_R_151 ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"GISS_E2_R_151/GISS_E2_R_151_lgm.nc"

        desc    = "GISS_E2_R_151 PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_GISS_E2_R_151-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="GISS_E2_R_151-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine GISS_E2_R_151_LGM_to_grid

    subroutine GISS_E2_R_151_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GISS_E2_R_151 ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"GISS_E2_R_151/GISS_E2_R_151_piControl.nc"

        desc    = "GISS_E2_R_151 PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_GISS_E2_R_151-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="GISS_E2_R_151-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine GISS_E2_R_151_PD_to_grid

    subroutine IPSL_CM5A_LR_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       IPSL_CM5A_LR ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"IPSL_CM5A_LR/IPSL_CM5A_LR_lgm.nc"

        desc    = "IPSL_CM5A_LR PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_IPSL_CM5A_LR-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="IPSL_CM5A_LR-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine IPSL_CM5A_LR_LGM_to_grid

    subroutine IPSL_CM5A_LR_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       IPSL_CM5A_LR ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"IPSL_CM5A_LR/IPSL_CM5A_LR_piControl.nc"

        desc    = "IPSL_CM5A_LR PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_IPSL_CM5A_LR-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="IPSL_CM5A_LR-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine IPSL_CM5A_LR_PD_to_grid

    subroutine MIROC_ESM_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MIROC_ESM ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"MIROC_ESM/MIROC_ESM_lgm.nc"

        desc    = "MIROC_ESM PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_MIROC_ESM-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="MIROC_ESM-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine MIROC_ESM_LGM_to_grid

    subroutine MIROC_ESM_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MIROC_ESM ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"MIROC_ESM/MIROC_ESM_piControl.nc"

        desc    = "MIROC_ESM PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_MIROC_ESM-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="MIROC_ESM-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine MIROC_ESM_PD_to_grid

    subroutine MPI_ESM_P_p1_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MPI_ESM_P_p1 ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"MPI_ESM_P_p1/MPI_ESM_P_p1_lgm.nc"

        desc    = "MPI_ESM_P_p1 PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_MPI_ESM_P_p1-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="MPI_ESM_P_p1-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine MPI_ESM_P_p1_LGM_to_grid

    subroutine MPI_ESM_P_p1_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MPI_ESM_P_p1 ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"MPI_ESM_P_p1/MPI_ESM_P_p1_piControl.nc"

        desc    = "MPI_ESM_P_p1 PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_MPI_ESM_P_p1-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="MPI_ESM_P_p1-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine MPI_ESM_P_p1_PD_to_grid

    subroutine MPI_ESM_P_p2_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MPI_ESM_P_p2 ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"MPI_ESM_P_p2/MPI_ESM_P_p2_lgm.nc"

        desc    = "MPI_ESM_P_p2 PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_MPI_ESM_P_p2-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="MPI_ESM_P_p2-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine MPI_ESM_P_p2_LGM_to_grid

    subroutine MPI_ESM_P_p2_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MPI_ESM_P_p2 ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"MPI_ESM_P_p2/MPI_ESM_P_p2_piControl.nc"

        desc    = "MPI_ESM_P_p2 PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_MPI_ESM_P_p2-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="MPI_ESM_P_p2-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine MPI_ESM_P_p2_PD_to_grid

    subroutine MRI_CGCM3_LGM_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MRI_CGCM3 ATM LGM data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"MRI_CGCM3/MRI_CGCM3_lgm.nc"

        desc    = "MRI_CGCM3 PMIP3 LGM ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_MRI_CGCM3-LGM.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="MRI_CGCM3-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine MRI_CGCM3_LGM_to_grid

    subroutine MRI_CGCM3_PD_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MRI_CGCM3 ATM PD data
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
        fldr_in          = trim(path_in)
        file_in          = trim(fldr_in)//"MRI_CGCM3/MRI_CGCM3_piControl.nc"

        desc    = "MRI_CGCM3 PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_MRI_CGCM3-PD.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="MRI_CGCM3-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(3))
        call def_var_info(vars(1),trim(file_in),"tas_spatialmean_ann","t2m_ann",units="C",long_name="Near-surface temperature (2-m), annual mean",method="nng")
        call def_var_info(vars(2),trim(file_in),"tas_spatialmean_djf","t2m_sum",units="C",long_name="Near-surface temperature (2-m), summer mean",method="nng")
        call def_var_info(vars(3),trim(file_in),"pr_spatialmean_ann","pr_ann",units="mm*d**-1",long_name="Precipitation, annual mean",method="nng") 
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine MRI_CGCM3_PD_to_grid

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

        ! For intermediate interpolation 
        character(len=256) :: pmip3_grid

        ! Define the input filenames
        fldr_in          = trim(path_in)//"/extensions/"
        file_in          = trim(fldr_in)//"pmip3_21k_v0.nc"

        ! Filename of already-processed present-day topography to combine with dzs here 
        file_in_topo     = trim(outfldr)//"/../"//trim(grid%name)//"_TOPO-RTOPO-2.0.1.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="pmip3_dzs_grid",mtype="latlon",units="degrees", &
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

            ! Load PD mask
            call nc_read(file_in_topo,"mask",outmask)

            ! Write output variable to output file
            call nc_write(filename,"mask",outmask,dim1="xc",dim2="yc")

            ! Write variable metadata
            call nc_write_attr(filename,"mask","units","")
            call nc_write_attr(filename,"mask","long_name","Mask (0:ocean, 1:land, 2:grounded ice, 3:floating ice)")
            call nc_write_attr(filename,"mask","coordinates","lat2D lon2D")

        else 

            ! Initialize mapping
            call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

            ! dz_srf =========

            ! Read in current variable
            call nc_read(trim(file_in),"orog_diff",inp%var,missing_value=mv)
            
            where((abs(inp%var) .ge. 1d10)) inp%var = mv

            ! Map variable to new grid
            call map_field(map,"dz_srf",inp%var,outvar,outmask,method="nng", &
                          fill=.TRUE.,missing_value=mv,sigma=sigma)

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

    subroutine LGM_dzs_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       LGM dzs wrt PD (Abe-Ouchi 2015)
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
        file_in          = trim(fldr_in)//"pmip3_21k_v0.nc"

        desc    = "MRI_CGCM3 PMIP3 PD ATM"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_PMIP3_21K.nc"

        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)

        call grid_init(grid0,name="dzs_LGM-grid",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat ) 

        ! jablasco
        allocate(vars(2))
        call def_var_info(vars(1),trim(file_in),"dzs","dzs",units="m",long_name="Surface elevation change LGM wrt PD",method="nng")
        call def_var_info(vars(2),trim(file_in),"mask","mask_bed",units="",long_name="Mask grounded ice",method="nn")
 
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

    end subroutine LGM_dzs_to_grid

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

        call grid_init(grid0,name="ext_LGM-grid",mtype="polar_stereographic",units="kilometers", &
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
