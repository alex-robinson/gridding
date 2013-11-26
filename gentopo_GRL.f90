
program gentopo

    use ncio 
    use coordinates

    implicit none

    type var_defs
        character(len=256) :: filename
        character(len=256) :: nm_in, nm_out  
        character(len=256) :: units_in, units_out 
        character(len=256) :: method
        logical :: mask, dimextra
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
    type(grid_class) :: gTOPO
    type(map_class)  :: mTOPO_ice, mTOPO_clim 
    character(len=256) :: file_topo
    type(var_defs), allocatable :: topo_invariant(:)

    ! == Additional variables ============================= 

    double precision, parameter :: missing_value = -9999.d0 

    double precision, dimension(:,:), allocatable :: invar, icevar, climvar, tmp 
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
                   lambda=-40.d0,phi=72.d0,alpha=7.5d0)

    ! ## Define clim grid and output variable field ##
    call grid_init(gclim,name="GRL-50KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   dx=50.d0,nx=37,dy=50.d0,ny=61, &
                   lambda=-40.d0,phi=72.d0,alpha=7.5d0)

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
    !       MAR (RCM) DATA - MARv3.2 interpolated data on ftp site
    !
    ! =========================================================

! ########################### 
    if (.FALSE.) then 

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
!     call grid_init(gMAR,name="MAR-5KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
!                    x0=-800.d0,dx=5.d0,nx=301,y0=-3400.d0,dy=5.d0,ny=561, &
!                    lambda=-39.d0,phi=90.d0,alpha=7.5d0)
    call grid_init(gMAR,name="MAR-25KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-750.d0,dx=25.d0,nx=58,y0=-1200.d0,dy=25.d0,ny=108, &
                   lambda=320.d0,phi=72.d0,alpha=7.5d0)
    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    call map_init(mMAR_ice, gMAR,gice, max_neighbors=30,lat_lim=3.d0,fldr="maps",load=.TRUE.)
    call map_init(mMAR_clim,gMAR,gclim,max_neighbors=30,lat_lim=3.d0,fldr="maps",load=.TRUE.)


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

    allocate(mar_surf(22))
    call def_var_info(mar_surf( 1),trim(file_mar),"SMB", "smb", units="mm month**-1")
    call def_var_info(mar_surf( 2),trim(file_mar),"RU",  "ru",  units="mm month**-1")
    call def_var_info(mar_surf( 3),trim(file_mar),"ME",  "me",  units="mm month**-1")
    call def_var_info(mar_surf( 4),trim(file_mar),"SMB2","smb2",units="mm month**-1")
    call def_var_info(mar_surf( 5),trim(file_mar),"SF",  "sf",  units="mm month**-1")
    call def_var_info(mar_surf( 6),trim(file_mar),"RF",  "rf",  units="mm month**-1")
    call def_var_info(mar_surf( 7),trim(file_mar),"SU",  "su",  units="mm month**-1")
    call def_var_info(mar_surf( 8),trim(file_mar),"SF",  "sf",  units="mm month**-1")
    call def_var_info(mar_surf( 9),trim(file_mar),"AL",  "al",  units="(0 - 1)")
    call def_var_info(mar_surf(10),trim(file_mar),"AL2", "al2", units="(0 - 1)")
    call def_var_info(mar_surf(11),trim(file_mar),"SF",  "sf",  units="mm month**-1")
    call def_var_info(mar_surf(12),trim(file_mar),"ST",  "Ts",  units="degrees Celcius")
    call def_var_info(mar_surf(13),trim(file_mar),"ST2", "Ts2", units="degrees Celcius")
    call def_var_info(mar_surf(14),trim(file_mar),"TT",  "T3m", units="degrees Celcius")
    call def_var_info(mar_surf(15),trim(file_mar),"SWD", "swd", units="W m**-2")
    call def_var_info(mar_surf(16),trim(file_mar),"LWD", "lwd", units="W m**-2")
    call def_var_info(mar_surf(17),trim(file_mar),"SHF", "shf", units="W m**-2")
    call def_var_info(mar_surf(18),trim(file_mar),"LHF", "lhf", units="W m**-2")
    call def_var_info(mar_surf(19),trim(file_mar),"SP",  "sp",  units="hPa")
    call def_var_info(mar_surf(20),trim(file_mar),"SMBcorr","smbc",units="mm month**-1")
    call def_var_info(mar_surf(21),trim(file_mar),"RUcorr", "ruc", units="mm month**-1")
    call def_var_info(mar_surf(22),trim(file_mar),"MEcorr", "mec", units="mm month**-1")
    
    ! (Re)Allocate the input grid variable
    call grid_allocate(gMAR,invar)

    ! ## INVARIANT FIELDS ##
    do i = 1, size(mar_invariant)
        var_now = mar_invariant(i) 
        call nc_read(var_now%filename,invar,var_now%nm_in,missing_value=missing_value)
        call map_field(mMAR_clim,var_now%nm_in,invar,climvar,climmask,var_now%method,100.d3, &
                      fill=.TRUE.,missing_value=missing_value)
        where(climvar .eq. missing_value) climvar = 0.d0 
        call nc_write(file_clim,climvar,var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)
        call map_field(mMAR_ice, var_now%nm_in,invar,icevar, icemask, var_now%method,50.d3, &
                       fill=.TRUE.,missing_value=missing_value)
        where(icevar .eq. missing_value) icevar = 0.d0  
        call nc_write(file_ice, icevar, var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)
    end do 

    nyr = 2012-1979+1
    nm  = 12 
       
    do k = 1, nyr 

        year = 1978 + k 
        write(*,*) 
        write(*,*) "=== ",year," ==="
        write(*,*)

        q = 0 
        do m = 1, nm 
            q = q+1 

            write(*,*)
            write(*,*) "= Month ",m, " ="
            write(*,*) 

            ! ## SURFACE FIELDS ##
            do i = 1, size(mar_surf)
                var_now = mar_surf(i) 
                write(var_now%filename,"(a,i4,a3)") trim(file_mar),year,".nc"

                call nc_read(var_now%filename,invar,var_now%nm_in,missing_value=missing_value, &
                             start=(/1,1,q/),count=(/gMAR%G%nx,gMAR%G%ny,1/))
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
    !       MAR (RCM) DATA - MARv3.2 original data passed by Xavier
    !
    ! =========================================================

! ########################### 
    if (.TRUE.) then 

    ! Define file names for input and output of global grids  
    file_ice       = "output/GRL-20KM_MARv3.2-ERA-INTERIM_197901-201112.nc"
    file_clim      = "output/GRL-50KM_MARv3.2-ERA-INTERIM_197901-201112.nc"

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
    call nc_write_dim(file_clim,"time", x=1979,dx=1,nx=33,units="years",calendar="360_day")
    
    call grid_write(gclim,file_clim,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gclim,climmask)
    call grid_allocate(gclim,climvar)


    ! Define MAR grid and input variable field
    call grid_init(gMAR,name="MAR-25KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-750.d0,dx=25.d0,nx=58,y0=-1200.d0,dy=25.d0,ny=108, &
                   lambda=-39.d0,phi=71.d0,alpha=7.5d0)
    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    call map_init(mMAR_ice, gMAR,gice, max_neighbors=20,lat_lim=5.d0,fldr="maps",load=.TRUE.)
    call map_init(mMAR_clim,gMAR,gclim,max_neighbors=20,lat_lim=5.d0,fldr="maps",load=.TRUE.)


    ! Define the variables to be mapped 
    file_mar = "data/MAR/MAR_ERA-INTERIM/MARv3.2_historical_mon_197901-197912.nc"

    allocate(mar_invariant(3))
    call def_var_info(mar_invariant(1),trim(file_mar),"SRF","mask_srf",units="(1 - 4)",method="nn")
    call def_var_info(mar_invariant(2),trim(file_mar),"SOL","mask_sol",units="(0 - 12)",method="nn")
    call def_var_info(mar_invariant(3),trim(file_mar),"SH","zs",  units="m")

    file_mar = "data/MAR/MAR_ERA-INTERIM/MARv3.2_historical_mon_"

    allocate(mar_surf(23))
    call def_var_info(mar_surf( 1),trim(file_mar),"SMB", "smb", units="mm d**-1",dimextra=.TRUE.)
    call def_var_info(mar_surf( 2),trim(file_mar),"RU",  "ru",  units="mm d**-1")
    call def_var_info(mar_surf( 3),trim(file_mar),"ME",  "me",  units="mm d**-1",dimextra=.TRUE.)
    call def_var_info(mar_surf( 4),trim(file_mar),"RZ",  "rz",  units="mm d**-1",dimextra=.TRUE.)
    call def_var_info(mar_surf( 5),trim(file_mar),"SF",  "sf",  units="mm d**-1")
    call def_var_info(mar_surf( 6),trim(file_mar),"RF",  "rf",  units="mm d**-1")
    call def_var_info(mar_surf( 7),trim(file_mar),"SU",  "su",  units="mm d**-1",dimextra=.TRUE.)
    call def_var_info(mar_surf( 8),trim(file_mar),"SF",  "sf",  units="mm d**-1")
    call def_var_info(mar_surf( 9),trim(file_mar),"TT",  "t3m", units="degrees Celcius",dimextra=.TRUE.)
    call def_var_info(mar_surf(10),trim(file_mar),"QQ",  "Q",   units="g kg**-1",dimextra=.TRUE.)
    call def_var_info(mar_surf(11),trim(file_mar),"UU",  "u",   units="m s**-1",dimextra=.TRUE.)
    call def_var_info(mar_surf(12),trim(file_mar),"VV",  "v",   units="m s**-1",dimextra=.TRUE.)
    call def_var_info(mar_surf(13),trim(file_mar),"SP",  "sp",  units="hPa")
    call def_var_info(mar_surf(14),trim(file_mar),"SWD", "swd", units="W m**-2")
    call def_var_info(mar_surf(15),trim(file_mar),"LWD", "lwd", units="W m**-2")
    call def_var_info(mar_surf(16),trim(file_mar),"LWU", "lwu", units="W m**-2")
    call def_var_info(mar_surf(17),trim(file_mar),"SHF", "shf", units="W m**-2")
    call def_var_info(mar_surf(18),trim(file_mar),"LHF", "lhf", units="W m**-2")
    call def_var_info(mar_surf(19),trim(file_mar),"AL1", "al1", units="(0 - 1)")
    call def_var_info(mar_surf(20),trim(file_mar),"AL2", "al2", units="(0 - 1)")
    call def_var_info(mar_surf(21),trim(file_mar),"CC",  "cc",  units="(0 - 1)")
    call def_var_info(mar_surf(22),trim(file_mar),"STT", "ts",  units="degrees Celcius",dimextra=.TRUE.)
    call def_var_info(mar_surf(23),trim(file_mar),"SHSN2","Hs", units="m",dimextra=.TRUE.)
    

    ! (Re)Allocate the input grid variable
    call grid_allocate(gMAR,invar)

    ! ## INVARIANT FIELDS ##
    do i = 1, size(mar_invariant)
        var_now = mar_invariant(i) 
        call nc_read(var_now%filename,invar,var_now%nm_in,missing_value=missing_value)
        climvar = missing_value 
        call map_field(mMAR_clim,var_now%nm_in,invar,climvar,climmask,var_now%method,100.d3, &
                      fill=.FALSE.,missing_value=missing_value)
        call nc_write(file_clim,climvar,var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)
        icevar = missing_value 
        call map_field(mMAR_ice, var_now%nm_in,invar,icevar, icemask, var_now%method,100.d3, &
                       fill=.FALSE.,missing_value=missing_value)  
        call nc_write(file_ice, icevar, var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)
    end do 

    nyr = 2011-1979+1
    nm  = 12 
       
    do k = 1, nyr 

        year = 1978 + k 
        write(*,*) 
        write(*,*) "=== ",year," ==="
        write(*,*)
 
        do m = 1, nm 
            q = m 

            write(*,*)
            write(*,*) "= Month ",m, " ="
            write(*,*) 

            ! ## SURFACE FIELDS ##
            do i = 1, size(mar_surf)
                var_now = mar_surf(i) 
                write(var_now%filename,"(a,i4,a3,i4,a5)") trim(file_mar),year,"01-",year,"12.nc"
                if (var_now%dimextra) then 
                    call nc_read(var_now%filename,invar,var_now%nm_in,missing_value=missing_value, &
                                  start=(/1,1,1,q/),count=(/gMAR%G%nx,gMAR%G%ny,1,1/))
                else 
                    call nc_read(var_now%filename,invar,var_now%nm_in,missing_value=missing_value, &
                             start=(/1,1,q/),count=(/gMAR%G%nx,gMAR%G%ny,1/))
                end if 
                climvar = missing_value 
                call map_field(mMAR_clim,var_now%nm_in,invar,climvar,climmask,"shepard",100.d3, &
                               fill=.FALSE.,missing_value=missing_value)
                call nc_write(file_clim,climvar,var_now%nm_out,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=(/1,1,m,k/),count=(/gclim%G%nx,gclim%G%ny,1,1/))
                icevar = missing_value 
                call map_field(mMAR_ice, var_now%nm_in,invar,icevar, icemask, "shepard",100.d3, &
                               fill=.FALSE.,missing_value=missing_value)
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

! ########################### 
    if (.FALSE.) then 

    ! Define file names for input and output of global grids  
    file_ice       = "output/GRL-20KM_TOPO.nc"
    file_clim      = "output/GRL-50KM_TOPO.nc"

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


    ! Define TOPO grid and input variable field
    call grid_init(gTOPO,name="TOPO-10KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-1300.d0,dx=10.d0,nx=251,y0=-3500.d0,dy=10.d0,ny=301, &
                   lambda=-39.d0,phi=90.d0,alpha=7.5d0)

    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    call map_init(mTOPO_ice, gTOPO,gice, max_neighbors=30,lat_lim=1d0,fldr="maps",load=.TRUE.)
    call map_init(mTOPO_clim,gTOPO,gclim,max_neighbors=30,lat_lim=1d0,fldr="maps",load=.TRUE.)

    ! Define the variables to be mapped 
    file_topo = "data/Greenland/Greenland_bedrock_topography_V3.nc"

    allocate(topo_invariant(4))
    call def_var_info(topo_invariant(1),trim(file_topo),"BedrockElevation","zb",units="m")
    call def_var_info(topo_invariant(2),trim(file_topo),"SurfaceElevation","zs",units="m")
    call def_var_info(topo_invariant(3),trim(file_topo),"IceThickness",    "H", units="m")
    call def_var_info(topo_invariant(4),trim(file_topo),"LandMask",      "mask",units="(0 - 4",method="nn")

    ! (Re)Allocate the input grid variable
    call grid_allocate(gTOPO,invar)

    ! Allocate tmp array to hold full data (trim to smaller size)
    allocate(tmp(2501,3001))

    ! ## INVARIANT FIELDS ##
    do i = 1, size(topo_invariant)
        var_now = topo_invariant(i) 
        call nc_read(var_now%filename,tmp,var_now%nm_in,missing_value=missing_value)
        call thin(invar,tmp,by=10)
        if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
            where( invar .eq. missing_value ) invar = 0.d0 
        end if
        if (trim(var_now%nm_out) .eq. "zb") then 
            call fill(invar,missing_value=missing_value,fill_value=-1000.d0)
        end if 
        call map_field(mTOPO_clim,var_now%nm_in,invar,climvar,climmask,var_now%method,20.d3, &
                      fill=.TRUE.,missing_value=missing_value)
        call nc_write(file_clim,climvar,var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)
        call map_field(mTOPO_ice, var_now%nm_in,invar,icevar, icemask, var_now%method,50.d3, &
                       fill=.TRUE.,missing_value=missing_value)
        call nc_write(file_ice, icevar, var_now%nm_out,  dim1="xc",dim2="yc",units=var_now%units_out)
    end do 

    end if 
! ########################### 
    

contains 

    subroutine def_var_info(var,filename,nm_in,nm_out,units,method,mask,dimextra)
        implicit none 

        type(var_defs) :: var 
        character(len=*) :: filename,nm_in,nm_out,units
        character(len=*), optional :: method 
        logical, optional :: mask, dimextra

        var%filename  = trim(filename)
        var%nm_in     = trim(nm_in)
        var%nm_out    = trim(nm_out)
        var%units_in  = trim(units)
        var%units_out = trim(units)

        var%method = "shepard"
        if (present(method)) var%method = trim(method)

        var%mask = .FALSE. 
        if (present(mask)) var%mask = mask 

        var%dimextra = .FALSE.
        if (present(dimextra)) var%dimextra = dimextra 

        return 

    end subroutine def_var_info

    subroutine thin(var1,var,by)
        implicit none

        double precision, dimension(:,:) :: var, var1 
        integer :: by 
        integer :: i,j, nx, ny 
        integer :: i1, j1

        nx = size(var,1)
        ny = size(var,2) 

        var1 = missing_value 

        i1 = 0
        do i = 1, nx, by 
            i1 = i1+1 
            j1 = 0 
            do j = 1, ny, by  
                j1 = j1 + 1 
                var1(i1,j1) = var(i,j)
            end do 
        end do 

        return
    end subroutine thin 

    subroutine fill(var,missing_value,fill_value)
        implicit none 
        double precision, dimension(:,:) :: var 
        double precision :: missing_value 
        double precision, optional :: fill_value

        integer :: q, nx, ny, i, j 
        integer, parameter :: qmax = 50 ! Iterations 

        double precision, dimension (3,3) :: neighb, weight
        double precision :: wtot, mval 
        double precision, dimension(:,:), allocatable :: filled
        nx = size(var,1)
        ny = size(var,2) 

        allocate(filled(nx,ny))

        if (present(fill_value)) then
            where(var .eq. missing_value) var = fill_value 
        end if 

        do q = 1, qmax 

            filled = missing_value 

            do i = 2, nx-1 
                do j = 2, ny-1 
                    neighb = var(i-1:i+1,j-1:j+1)

                    weight = 0.d0 
                    where (neighb .ne. missing_value) weight = 1.d0
                    wtot = sum(weight)

                    if (wtot .gt. 0.d0) then 
                        mval = sum(neighb*weight)/wtot
                        where (neighb .eq. missing_value) neighb = mval 
                    end if 

                    filled(i-1:i+1,j-1:j+1) = neighb 

                end do 
            end do 

            where(filled .ne. missing_value) var = filled 

            write(*,*) q," : Missing values: ", count(var .eq. missing_value)
            if ( count(var .eq. missing_value) .eq. 0 ) exit 
        end do 

        return
    end subroutine fill 

end program gentopo