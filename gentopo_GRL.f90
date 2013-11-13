
program gentopo

    use ncio 
    use coordinates

    implicit none

    type(grid_class) :: gice, gclim
    type(grid_class) :: gTOPO, gECMWF, gMAR
    type(map_class)  :: mTOPO_ice,  mTOPO_clim 
    type(map_class)  :: mECMWF_ice, mECMWF_clim
    type(map_class)  :: mMAR_ice,   mMAR_clim 

    type topo_vars 
        double precision, dimension(:,:), allocatable :: zs, zb, H
        integer,          dimension(:,:), allocatable :: mask
    end type 
    type(topo_vars) :: topo 

    type ecmwf_vars 
        double precision, dimension(:,:), allocatable :: zs

    end type 
    type(ecmwf_vars) :: ecm 

    type mar_vars 
        double precision, dimension(:,:), allocatable :: zs

    end type 
    type(mar_vars) :: mar 

    type ice_vars 
        double precision, dimension(:,:), allocatable :: zs, zb, H
        integer,          dimension(:,:), allocatable :: mask_interp
    end type 
    type(ice_vars) :: ice 

    character(len=256) :: file_input, file_output, file_test, var_name
    character(len=256) :: file_ice, file_topo
    integer :: t, ttot

    double precision, parameter :: missing_value = -9999.d0 

    ! =======================================================================
    !
    ! Step 1: Define global input grid and load data that will be used here
    !
    ! =======================================================================

    ! Define file names for input and output of global grids  
    file_ice       = "output/GRL-20KM_topo.nc"

    file_topo      = "data/Greenland/Greenland_bedrock_topography_V3.nc"

    ! Define MAR grid and input variable field
!     call grid_init(gMAR,name="MAR-25KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
!                    x0=-750.d0,dx=25.d0,nx=58,y0=-1200.d0,dy=25.d0,ny=108, &
!                    lambda=320.d0,phi=72.d0,alpha=7.5d0)
    

    ! Define ice grid and output variable field
    call grid_init(gice,name="GRL-20KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   dx=20.d0,nx=90,dy=20.d0,ny=150, &
                   lambda=-39.d0,phi=71.d0,alpha=7.5d0)

    ! Output new grid to file
    call nc_create(file_ice)
    call nc_write_dim(file_ice,"xc",  x=gice%G%x,units="kilometers")
    call nc_write_dim(file_ice,"yc",  x=gice%G%y,units="kilometers")
    !call nc_write_dim(file_ice,"time",x=1.d0,dx=1.d0,nx=365,units="years",calendar="365_day")
    
    call grid_write(gice,file_ice,xnm="xc",ynm="yc",create=.FALSE.)

    ! Allocate new ice variables 
    call grid_allocate(gice,ice%zs)
    call grid_allocate(gice,ice%zb)
    call grid_allocate(gice,ice%H)
    call grid_allocate(gice,ice%mask_interp)
    
    ! Define TOPO grid and input variable field
    call grid_init(gTOPO,name="TOPO-1KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-1300.d0,dx=1.d0,nx=2501,y0=-3500.d0,dy=1.d0,ny=3001, &
                   lambda=-39.d0,phi=90.d0,alpha=7.5d0)

!     call grid_write(gTOPO,"output/test_TOPO.nc",xnm="xc",ynm="yc",create=.TRUE.)
    
    ! Load data 
    call grid_allocate(gTOPO,TOPO%zs)
    call grid_allocate(gTOPO,TOPO%zb)
    call grid_allocate(gTOPO,TOPO%H)

    call nc_read(file_topo,TOPO%zs,"SurfaceElevation",missing_value=missing_value)
    call nc_read(file_topo,TOPO%zb,"BedrockElevation",missing_value=missing_value)
    call nc_read(file_topo,TOPO%H, "IceThickness",    missing_value=missing_value)
    
    ! =======================================================================
    !
    ! Step 2: Map the fields 
    !
    ! =======================================================================
    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    ! max_neighbors is the maximum neighbors to be stored for each point
    ! lat_lim is the range of latitudes relative to a given point to check neighbor distances (to speed things up)
!     call map_init(mTOPO_ice,gTOPO,gice,max_neighbors=50,lat_lim=0.4d0,fldr="maps",load=.TRUE.)

!     ! Map each field back to the SICO domain using the radius method
!     call map_field(mTOPO_ice,"zs",TOPO%zs,ice%zs,ice%mask_interp,"shepard",20.d3,fill=.TRUE.,missing_value=missing_value)
!     call map_field(mTOPO_ice,"zb",TOPO%zb,ice%zb,ice%mask_interp,"shepard",20.d3,fill=.TRUE.,missing_value=missing_value)
!     call map_field(mTOPO_ice, "H",TOPO%H, ice%H, ice%mask_interp,"shepard",20.d3,fill=.TRUE.,missing_value=missing_value)

!     ! Write new gridded ice data to grid file 
!     call nc_write(file_ice,ice%zs,  "zs", dim1="xc",dim2="yc")
!     call nc_write(file_ice,ice%zb,  "zb", dim1="xc",dim2="yc")
!     call nc_write(file_ice,ice%H,    "H", dim1="xc",dim2="yc")

end program gentopo