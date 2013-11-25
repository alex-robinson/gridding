
program gentopo

    use ncio 
    use coordinates

    implicit none

    type var_defs
        character(len=256) :: filename
        character(len=256) :: nm_in, nm_out  
        character(len=256) :: units_in, units_out 
        character(len=256) :: method
        logical :: mask 
    end type 
    type(var_defs) :: var_now

    type(grid_class) :: gice, gclim
    character(len=256) :: file_ice, file_clim
    
    ! ECMWF-related variables
    type(grid_class)   :: gECMWF 
    type(map_class)    :: mECMWF_ice, mECMWF_clim
    character(len=256) :: file_invariant, file_surface, file_pres 
    type(var_defs), allocatable :: ecmwf_invariant(:), ecmwf_surf(:), ecmwf_pres(:) 

    ! MAR-related variables
    type(grid_class)   :: gMAR
    type(map_class)    :: mMAR_ice, mMAR_clim 
    character(len=256) :: file_mar 
    type(var_defs), allocatable :: mar_invariant(:), mar_surf(:)

    ! TOPO-related variables
!     type(grid_class) :: gTOPO
!     type(map_class)  :: mTOPO_ice, mTOPO_clim 
    

    ! == Additional variables ============================= 

    double precision, parameter :: missing_value = -9999.d0 

    double precision, dimension(:,:), allocatable :: invar, icevar, climvar
    integer,          dimension(:,:), allocatable :: icemask, climmask

    integer :: i, q, k, m, nyr, nm, c
    integer :: year 

    ! =======================================================================
    !
    ! Step 1: Define grids for output data: ice and clim
    !
    ! =======================================================================

    ! ## Define ice grid and output variable field ##
    call grid_init(gice,name="GRL-20KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   dx=20.d0,nx=90,dy=20.d0,ny=150, &
                   lambda=-39.d0,phi=71.d0,alpha=7.5d0)

    ! ## Define clim grid and output variable field ##
    call grid_init(gclim,name="GRL-50KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   dx=50.d0,nx=37,dy=50.d0,ny=61, &
                   lambda=-39.d0,phi=71.d0,alpha=7.5d0)

    ! =========================================================
    !
    !       ECMWF DATA
    !
    ! =========================================================

! ###########################   
    if (.FALSE.) then 

    ! Define file names for input and output of global grids  
    file_ice       = "output/GRL-20KM_ERA-INTERIM-750Mb.nc"
    file_clim      = "output/GRL-50KM_ERA-INTERIM-750Mb.nc"
    
    ! Write ice grid to file
    call nc_create(file_ice)
    call nc_write_dim(file_ice,"xc",  x=gice%G%x,units="kilometers")
    call nc_write_dim(file_ice,"yc",  x=gice%G%y,units="kilometers")
    call nc_write_dim(file_ice,"month",x=(/1,2,3,4,5,6,7,8,9,10,11,12/),units="month")
    call nc_write_dim(file_ice,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gice,file_ice,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gice,icemask)
    call grid_allocate(gice,icevar)

    ! Write clim grid to file
    call nc_create(file_clim)
    call nc_write_dim(file_clim,"xc",   x=gclim%G%x,units="kilometers")
    call nc_write_dim(file_clim,"yc",   x=gclim%G%y,units="kilometers")
    call nc_write_dim(file_clim,"month",x=(/1,2,3,4,5,6,7,8,9,10,11,12/),units="month")
    call nc_write_dim(file_clim,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gclim,file_clim,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gclim,climmask)
    call grid_allocate(gclim,climvar)

    ! Define ECMWF input grid
    call grid_init(gECMWF,name="ECMWF-GRL075",mtype="latlon",units="kilometers",lon180=.TRUE., &
                   x0=-100.d0,dx=0.75d0,nx=161,y0=49.5d0,dy=0.75d0,ny=55)
    
    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    call map_init(mECMWF_ice, gECMWF,gice, max_neighbors=20,lat_lim=5.d0,fldr="maps",load=.TRUE.)
    call map_init(mECMWF_clim,gECMWF,gclim,max_neighbors=20,lat_lim=5.d0,fldr="maps",load=.TRUE.)


    ! Define the variables to be mapped 
    file_invariant = "data/ECMWF/NEW/ERA-INTERIM-GRL-invariant_historical_mon_197901-201212.nc"
    file_surface   = "data/ECMWF/NEW/ERA-INTERIM-GRL-surface_historical_mon_197901-201212.nc"
    file_pres      = "data/ECMWF/NEW/ERA-INTERIM-GRL-750Mb_historical_mon_197901-201212.nc"

    allocate(ecmwf_invariant(1))
    call def_var_info(ecmwf_invariant(1),trim(file_invariant),"z","zs",units="m")

    allocate(ecmwf_surf(8))
    call def_var_info(ecmwf_surf(1),trim(file_surface),"sp", "sp", units="Pa")
    call def_var_info(ecmwf_surf(2),trim(file_surface),"tcw","tcw",units="kg m**-2")
    call def_var_info(ecmwf_surf(3),trim(file_surface),"tcc","tcc",units="(0 - 1)")
    call def_var_info(ecmwf_surf(4),trim(file_surface),"u10","u10",units="m s**-1")
    call def_var_info(ecmwf_surf(5),trim(file_surface),"v10","v10",units="m s**-1")
    call def_var_info(ecmwf_surf(6),trim(file_surface),"t2m","t2m",units="K")
    call def_var_info(ecmwf_surf(7),trim(file_surface),"al", "al", units="(0 - 1)")
    call def_var_info(ecmwf_surf(8),trim(file_surface),"sst","sst",units="K")
    
    allocate(ecmwf_pres(7))
    call def_var_info(ecmwf_pres(1),trim(file_pres),"z", "p_z",units="m**2 s**-2")
    call def_var_info(ecmwf_pres(2),trim(file_pres),"t", "p_t",units="K")
    call def_var_info(ecmwf_pres(3),trim(file_pres),"q", "p_q",units="kg kg**-1")
    call def_var_info(ecmwf_pres(4),trim(file_pres),"w", "p_w",units="Pa s**-1")
    call def_var_info(ecmwf_pres(5),trim(file_pres),"r", "p_r",units="%")
    call def_var_info(ecmwf_pres(6),trim(file_pres),"u", "p_u",units="m s**-1")
    call def_var_info(ecmwf_pres(7),trim(file_pres),"v", "p_v",units="m s**-1")

    ! Allocate the input grid variable
    call grid_allocate(gECMWF,invar)

    ! ## INVARIANT FIELDS ##
    var_now = ecmwf_invariant(1) 
    call nc_read(var_now%filename,invar,var_now%nm_in)
    call map_field(mECMWF_clim,var_now%nm_in,invar,climvar,climmask,"shepard",400.d3,missing_value=missing_value)
    call nc_write(file_clim,climvar,var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)
    call map_field(mECMWF_ice, var_now%nm_in,invar,icevar, icemask, "shepard",400.d3,missing_value=missing_value)
    call nc_write(file_ice, icevar, var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)

    nyr = 2012-1979+1
    nm  = 12 

    q = 0 
    do k = 1, nyr 

        year = 1978 + k 
        write(*,*) 
        write(*,*) "=== ",year," ==="
        write(*,*)

        do m = 1, nm 
            q = q+1 

            write(*,*)
            write(*,*) "= Month ",m, " ="
            write(*,*) 

            ! ## SURFACE FIELDS ##
            do i = 1, size(ecmwf_surf)
                var_now = ecmwf_surf(i) 
                call nc_read(var_now%filename,invar,var_now%nm_in,start=(/1,1,q/),count=(/gECMWF%G%nx,gECMWF%G%ny,1/))
                call map_field(mECMWF_clim,var_now%nm_in,invar,climvar,climmask,"shepard",400.d3,missing_value=missing_value)
                call nc_write(file_clim,climvar,var_now%nm_out,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=(/1,1,m,k/),count=(/gclim%G%nx,gclim%G%ny,1,1/))
                call map_field(mECMWF_ice, var_now%nm_in,invar,icevar, icemask, "shepard",400.d3,missing_value=missing_value)
                call nc_write(file_ice,icevar,var_now%nm_out,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=(/1,1,m,k/),count=(/gice%G%nx,gice%G%ny,1,1/))
            end do 

            ! ## PRESSURE FIELDS ##
            do i = 1, size(ecmwf_pres)
                var_now = ecmwf_pres(i) 
                call nc_read(var_now%filename,invar,var_now%nm_in,start=(/1,1,q/),count=(/gECMWF%G%nx,gECMWF%G%ny,1/))
                call map_field(mECMWF_clim,var_now%nm_in,invar,climvar,climmask,"shepard",400.d3,missing_value=missing_value)
                call nc_write(file_clim,climvar,var_now%nm_out,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=(/1,1,m,k/),count=(/gclim%G%nx,gclim%G%ny,1,1/))
                call map_field(mECMWF_ice, var_now%nm_in,invar,icevar, icemask, "shepard",400.d3,missing_value=missing_value)
                call nc_write(file_ice,icevar,var_now%nm_out,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=(/1,1,m,k/),count=(/gice%G%nx,gice%G%ny,1,1/))
            end do 

        end do 
    end do 

    end if 
! ########################### 

    ! =========================================================
    !
    !       MAR (RCM) DATA
    !
    ! =========================================================

    ! Define file names for input and output of global grids  
    file_ice       = "output/GRL-20KM_MARv3.2-ERA-INTERIM_197901-201212.nc"
    file_clim      = "output/GRL-50KM_MARv3.2-ERA-INTERIM_197901-201212.nc"

    ! Write ice grid to file
    call nc_create(file_ice)
    call nc_write_dim(file_ice,"xc",  x=gice%G%x,units="kilometers")
    call nc_write_dim(file_ice,"yc",  x=gice%G%y,units="kilometers")
    call nc_write_dim(file_ice,"month",x=(/1,2,3,4,5,6,7,8,9,10,11,12/),units="month")
    call nc_write_dim(file_ice,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gice,file_ice,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gice,icemask)
    call grid_allocate(gice,icevar)

    ! Write clim grid to file
    call nc_create(file_clim)
    call nc_write_dim(file_clim,"xc",   x=gclim%G%x,units="kilometers")
    call nc_write_dim(file_clim,"yc",   x=gclim%G%y,units="kilometers")
    call nc_write_dim(file_clim,"month",x=(/1,2,3,4,5,6,7,8,9,10,11,12/),units="month")
    call nc_write_dim(file_clim,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gclim,file_clim,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gclim,climmask)
    call grid_allocate(gclim,climvar)


    ! Define MAR grid and input variable field
    call grid_init(gMAR,name="MAR-5KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-800.d0,dx=5.d0,nx=301,y0=-3400.d0,dy=5.d0,ny=561, &
                   lambda=-39.d0,phi=90.d0,alpha=7.5d0)
    
    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    call map_init(mMAR_ice, gMAR,gice, max_neighbors=20,lat_lim=1.d0,fldr="maps",load=.TRUE.)
    call map_init(mMAR_clim,gMAR,gclim,max_neighbors=20,lat_lim=1.d0,fldr="maps",load=.TRUE.)


    ! Define the variables to be mapped 
    file_mar = "data/MAR/MARv3.2/ERA-Interim_1979-2012/MARv3.2-monthly-ERA-Interim-2012.nc"

    allocate(mar_invariant(6))
    call def_var_info(mar_invariant(1),trim(file_mar),"MSK_bam01","mask_b01",units="(0 - 2)",method="nn")
    call def_var_info(mar_invariant(2),trim(file_mar),"SRF_bam01","zs_b01",  units="m")
    call def_var_info(mar_invariant(3),trim(file_mar),"MSK_bam13","mask_b13",units="(0 - 2)",method="nn")
    call def_var_info(mar_invariant(4),trim(file_mar),"SRF_bam13","zs_b13",  units="m")
    call def_var_info(mar_invariant(5),trim(file_mar),"MSK_MAR","mask_mar",units="(0 - 2)",method="nn")
    call def_var_info(mar_invariant(6),trim(file_mar),"SRF_MAR","zs_mar",  units="m")

    file_mar = "data/MAR/MARv3.2/ERA-Interim_1979-2012/MARv3.2-monthly-ERA-Interim-"

    allocate(mar_surf(23))
    call def_var_info(mar_surf( 1),trim(file_mar),"SMB", "smb", units="mm month**-1")
    call def_var_info(mar_surf( 2),trim(file_mar),"RU",  "ru",  units="mm month**-1")
    call def_var_info(mar_surf( 3),trim(file_mar),"ME",  "me",  units="mm month**-1")
    call def_var_info(mar_surf( 4),trim(file_mar),"SMB2","smb2",units="mm month**-1")
    call def_var_info(mar_surf( 5),trim(file_mar),"SF",  "sf",  units="mm month**-1")
    call def_var_info(mar_surf( 6),trim(file_mar),"RF",  "rf",  units="mm month**-1")
    call def_var_info(mar_surf( 7),trim(file_mar),"SU",  "su",  units="mm month**-1")
    call def_var_info(mar_surf( 8),trim(file_mar),"SF",  "sf",  units="mm month**-1")
    call def_var_info(mar_surf( 9),trim(file_mar),"AL",  "al",  units="mm month**-1")
    call def_var_info(mar_surf(10),trim(file_mar),"AL2", "al2", units="mm month**-1")
    call def_var_info(mar_surf(11),trim(file_mar),"SF",  "sf",  units="mm month**-1")
    call def_var_info(mar_surf(12),trim(file_mar),"ST",  "Ts",  units="mm month**-1")
    call def_var_info(mar_surf(13),trim(file_mar),"ST2", "Ts2", units="mm month**-1")
    call def_var_info(mar_surf(14),trim(file_mar),"SF",  "sf",  units="mm month**-1")
    call def_var_info(mar_surf(15),trim(file_mar),"TT",  "T3m", units="degrees Celcius")
    call def_var_info(mar_surf(16),trim(file_mar),"SWD", "swd", units="W m**-2")
    call def_var_info(mar_surf(17),trim(file_mar),"LWD", "lwd", units="W m**-2")
    call def_var_info(mar_surf(18),trim(file_mar),"SHF", "shf", units="W m**-2")
    call def_var_info(mar_surf(19),trim(file_mar),"LHF", "lhf", units="W m**-2")
    call def_var_info(mar_surf(20),trim(file_mar),"SP",  "sp",  units="hPa")
    call def_var_info(mar_surf(21),trim(file_mar),"SMBc","smbc",units="mm month**-1")
    call def_var_info(mar_surf(22),trim(file_mar),"RUc", "ruc", units="mm month**-1")
    call def_var_info(mar_surf(23),trim(file_mar),"MEc", "mec", units="mm month**-1")
    
    ! (Re)Allocate the input grid variable
    call grid_allocate(gMAR,invar)

    ! ## INVARIANT FIELDS ##
    do i = 1, size(mar_invariant)
        var_now = mar_invariant(i) 
        call nc_read(var_now%filename,invar,var_now%nm_in,missing_value=missing_value)
        call map_field(mMAR_clim,var_now%nm_in,invar,climvar,climmask,var_now%method,100.d3, &
                      fill=.TRUE.,missing_value=missing_value)
        where(invar .eq. missing_value) invar = 0.d0 
        call nc_write(file_clim,climvar,var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)
        call map_field(mMAR_ice, var_now%nm_in,invar,icevar, icemask, var_now%method,50.d3, &
                       fill=.TRUE.,missing_value=missing_value)
        where(invar .eq. missing_value) invar = 0.d0 
        call nc_write(file_ice, icevar, var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)
    end do 

! ########################### 
    if (.TRUE.) then 
       
    nyr = 2012-1979+1
    nm  = 12 

    q = 0 
    do k = 1, nyr 

        year = 1978 + k 
        write(*,*) 
        write(*,*) "=== ",year," ==="
        write(*,*)

        var_now%filename = ""
        write(var_now%filename,"(a,i4,a3)") trim(file_mar),year,".nc"
        write(*,*) "|"//trim(var_now%filename)//"|"

        do m = 1, nm 
            q = q+1 

            write(*,*)
            write(*,*) "= Month ",m, " ="
            write(*,*) 

            ! ## SURFACE FIELDS ##
            do i = 1, size(mar_surf)
                var_now = mar_surf(i) 
                call nc_read(var_now%filename,invar,var_now%nm_in,start=(/1,1,q/),count=(/gMAR%G%nx,gMAR%G%ny,1/))
                call map_field(mMAR_clim,var_now%nm_in,invar,climvar,climmask,"shepard",100.d3, &
                               fill=.TRUE.,missing_value=missing_value)
                call nc_write(file_clim,climvar,var_now%nm_out,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=(/1,1,m,k/),count=(/gclim%G%nx,gclim%G%ny,1,1/))
                call map_field(mMAR_ice, var_now%nm_in,invar,icevar, icemask, "shepard",50.d3, &
                               fill=.TRUE.,missing_value=missing_value)
                call nc_write(file_ice,icevar,var_now%nm_out,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=(/1,1,m,k/),count=(/gice%G%nx,gice%G%ny,1,1/))
            end do 

        end do 
    end do 

    end if 
! ########################### 


    ! =========================================================
    !
    !       TOPO DATA
    !
    ! =========================================================

    !     file_topo      = "data/Greenland/Greenland_bedrock_topography_V3.nc"

    ! Define TOPO grid and input variable field
!     call grid_init(gTOPO,name="TOPO-1KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
!                    x0=-1300.d0,dx=1.d0,nx=2501,y0=-3500.d0,dy=1.d0,ny=3001, &
!                    lambda=-39.d0,phi=90.d0,alpha=7.5d0)

!     call grid_write(gTOPO,"output/test_TOPO.nc",xnm="xc",ynm="yc",create=.TRUE.)
    
    ! Load data 
!     call grid_allocate(gTOPO,TOPO%zs)
!     call grid_allocate(gTOPO,TOPO%zb)
!     call grid_allocate(gTOPO,TOPO%H)

!     call nc_read(file_topo,TOPO%zs,"SurfaceElevation",missing_value=missing_value)
!     call nc_read(file_topo,TOPO%zb,"BedrockElevation",missing_value=missing_value)
!     call nc_read(file_topo,TOPO%H, "IceThickness",    missing_value=missing_value)
    
!     ! Initialize 'to' and 'fro' mappings
!     ! max_neighbors is the maximum neighbors to be stored for each point
!     ! lat_lim is the range of latitudes relative to a given point to check neighbor distances (to speed things up)
!     call map_init(mTOPO_ice,gTOPO,gice,max_neighbors=50,lat_lim=0.4d0,fldr="maps",load=.TRUE.)

!     ! Map each field back to the SICO domain using the radius method
!     call map_field(mTOPO_ice,"zs",TOPO%zs,ice%zs,ice%mask_interp,"shepard",20.d3,fill=.TRUE.,missing_value=missing_value)
!     call map_field(mTOPO_ice,"zb",TOPO%zb,ice%zb,ice%mask_interp,"shepard",20.d3,fill=.TRUE.,missing_value=missing_value)
!     call map_field(mTOPO_ice, "H",TOPO%H, ice%H, ice%mask_interp,"shepard",20.d3,fill=.TRUE.,missing_value=missing_value)

!     ! Write new gridded ice data to grid file 
!     call nc_write(file_ice,ice%zs,  "zs", dim1="xc",dim2="yc")
!     call nc_write(file_ice,ice%zb,  "zb", dim1="xc",dim2="yc")
!     call nc_write(file_ice,ice%H,    "H", dim1="xc",dim2="yc")

contains 

    subroutine def_var_info(var,filename,nm_in,nm_out,units,method,mask)
        implicit none 

        type(var_defs) :: var 
        character(len=*) :: filename,nm_in,nm_out,units
        character(len=*), optional :: method 
        logical, optional :: mask 

        var%filename  = trim(filename)
        var%nm_in     = trim(nm_in)
        var%nm_out    = trim(nm_out)
        var%units_in  = trim(units)
        var%units_out = trim(units)

        var%method = "shepard"
        if (present(method)) var%method = trim(method)

        var%mask = .FALSE. 
        if (present(mask)) var%mask = mask 

        return 

    end subroutine def_var_info

end program gentopo