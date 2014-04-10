
program gentopo

    use ncio 
    use coordinates
    use vargrid 

    implicit none

    type(grid_class) :: gice, gclim
    character(len=256) :: file_ice, file_clim
    
    ! ECMWF-related variables
    type(grid_class)   :: gECMWF 
    type(map_class)    :: mECMWF_ice, mECMWF_clim
    character(len=256) :: file_invariant, file_surface
    character(len=256) :: files_pres(9)
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

    ! OLR-related variables 
    type(grid_class) :: gOLR 
    type(map_class)  :: mOLR_ice, mOLR_clim 
    character(len=256) :: file_olr 
    type(var_defs) :: olr_toa 

    ! CERES-related variables 
    type(grid_class) :: gCERES 
    type(map_class)  :: mCERES_ice, mCERES_clim 
    character(len=256) :: file_ceres 
    type(var_defs), allocatable :: ceres(:)

    ! == Additional variables =============================  

    type(var_defs) :: var_now

    double precision, dimension(:,:), allocatable :: invar, icevar, climvar, tmp 
    integer,          dimension(:,:), allocatable :: icemask, climmask

    integer :: i, q, k, m, nyr, nm, c, l 
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
    file_ice       = "output/GRL-20KM_ERA-INTERIM_mon_197901-201212.nc"
    file_clim      = "output/GRL-50KM_ERA-INTERIM_mon_197901-201212.nc"
    
    ! Write ice grid to file
    call nc_create(file_ice)
    call nc_write_dim(file_ice,"xc",  x=gice%G%x,units="kilometers")
    call nc_write_dim(file_ice,"yc",  x=gice%G%y,units="kilometers")
    call nc_write_dim(file_ice,"plev", x=[1000.d0,950.d0,850.d0,750.d0,700.d0,650.d0,600.d0,550.d0,500.d0],units="hPa")
    call nc_write_dim(file_ice,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
    call nc_write_dim(file_ice,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gice,file_ice,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gice,icemask)
    call grid_allocate(gice,icevar)

    ! Write clim grid to file
    call nc_create(file_clim)
    call nc_write_dim(file_clim,"xc",   x=gclim%G%x,units="kilometers")
    call nc_write_dim(file_clim,"yc",   x=gclim%G%y,units="kilometers")
    call nc_write_dim(file_clim,"plev", x=[1000.d0,950.d0,850.d0,750.d0,700.d0,650.d0,600.d0,550.d0,500.d0],units="hPa")
    call nc_write_dim(file_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
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
    files_pres(1)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-1000Mb_historical_mon_197901-201212.nc"
    files_pres(2)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-950Mb_historical_mon_197901-201212.nc"
    files_pres(3)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-850Mb_historical_mon_197901-201212.nc"
    files_pres(4)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-750Mb_historical_mon_197901-201212.nc"
    files_pres(5)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-700Mb_historical_mon_197901-201212.nc"
    files_pres(6)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-650Mb_historical_mon_197901-201212.nc"
    files_pres(7)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-600Mb_historical_mon_197901-201212.nc"
    files_pres(8)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-550Mb_historical_mon_197901-201212.nc"
    files_pres(9)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-500Mb_historical_mon_197901-201212.nc"

    allocate(ecmwf_invariant(1))
    call def_var_info(ecmwf_invariant(1),trim(file_invariant),"z","zs",units="m")

    allocate(ecmwf_surf(12))
    call def_var_info(ecmwf_surf( 1),trim(file_surface),"sp", "sp", units="Pa")
    call def_var_info(ecmwf_surf( 2),trim(file_surface),"tcw","tcw",units="kg m**-2")
    call def_var_info(ecmwf_surf( 3),trim(file_surface),"tclw","tclw",units="kg m**-2")
    call def_var_info(ecmwf_surf( 4),trim(file_surface),"tciw","tciw",units="kg m**-2")
    call def_var_info(ecmwf_surf( 5),trim(file_surface),"p56.162","clw",units="kg m**-2")
    call def_var_info(ecmwf_surf( 6),trim(file_surface),"p57.162","ciw",units="kg m**-2")
    call def_var_info(ecmwf_surf( 7),trim(file_surface),"tcc","tcc",units="(0 - 1)")
    call def_var_info(ecmwf_surf( 8),trim(file_surface),"u10","u10",units="m s**-1")
    call def_var_info(ecmwf_surf( 9),trim(file_surface),"v10","v10",units="m s**-1")
    call def_var_info(ecmwf_surf(10),trim(file_surface),"t2m","t2m",units="K")
    call def_var_info(ecmwf_surf(11),trim(file_surface),"al", "al", units="(0 - 1)")
    call def_var_info(ecmwf_surf(12),trim(file_surface),"sst","sst",units="K")
    
    allocate(ecmwf_pres(7))
    call def_var_info(ecmwf_pres( 1),"None","t", "p_t",units="K",         plev="plev",filenames=files_pres)
    call def_var_info(ecmwf_pres( 2),"None","q", "p_q",units="kg kg**-1", plev="plev",filenames=files_pres)
    call def_var_info(ecmwf_pres( 3),"None","r", "p_r",units="%",         plev="plev",filenames=files_pres)
    call def_var_info(ecmwf_pres( 4),"None","z", "p_z",units="m**2 s**-2",plev="plev",filenames=files_pres)
    call def_var_info(ecmwf_pres( 5),"None","w", "p_w",units="Pa s**-1",  plev="plev",filenames=files_pres)
    call def_var_info(ecmwf_pres( 6),"None","u", "p_u",units="m s**-1",   plev="plev",filenames=files_pres)
    call def_var_info(ecmwf_pres( 7),"None","v", "p_v",units="m s**-1",   plev="plev",filenames=files_pres)

    ! Allocate the input grid variable
    call grid_allocate(gECMWF,invar)

    ! ## INVARIANT FIELDS ##
    var_now = ecmwf_invariant(1) 
    call nc_read(var_now%filename,var_now%nm_in,invar)
    call map_field(mECMWF_clim,var_now%nm_in,invar,climvar,climmask,"shepard",400.d3,missing_value=missing_value)
    call nc_write(file_clim,var_now%nm_out,climvar,  dim1="xc",dim2="yc",units=var_now%units_out)
    call map_field(mECMWF_ice, var_now%nm_in,invar,icevar, icemask, "shepard",400.d3,missing_value=missing_value)
    call nc_write(file_ice, var_now%nm_out,icevar,   dim1="xc",dim2="yc",units=var_now%units_out)

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
                call nc_read(var_now%filename,var_now%nm_in,invar,start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1])
                call map_field(mECMWF_clim,var_now%nm_in,invar,climvar,climmask,"shepard",400.d3,missing_value=missing_value)
                call nc_write(file_clim,var_now%nm_out,climvar,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=[1,1,m,k],count=[gclim%G%nx,gclim%G%ny,1,1])
                call map_field(mECMWF_ice, var_now%nm_in,invar,icevar, icemask, "shepard",400.d3,missing_value=missing_value)
                call nc_write(file_ice,var_now%nm_out,icevar,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=[1,1,m,k],count=[gice%G%nx,gice%G%ny,1,1])
            end do 

            ! ## PRESSURE FIELDS ##
            do i = 1, size(ecmwf_pres)
                var_now = ecmwf_pres(i) 

                do l = 1, size(files_pres)   ! Loop over pressure layers
                    call nc_read(var_now%filenames(l),var_now%nm_in,invar,start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1])
                    call map_field(mECMWF_clim,var_now%nm_in,invar,climvar,climmask,"shepard",400.d3,missing_value=missing_value)
                    call nc_write(file_clim,var_now%nm_out,climvar,  dim1="xc",dim2="yc",dim3="plev",dim4="month",dim5="time", &
                                  units=var_now%units_out,start=[1,1,l,m,k],count=[gclim%G%nx,gclim%G%ny,1,1,1])
                    call map_field(mECMWF_ice, var_now%nm_in,invar,icevar, icemask, "shepard",400.d3,missing_value=missing_value)
                    call nc_write(file_ice,var_now%nm_out,icevar,  dim1="xc",dim2="yc",dim3="plev",dim4="month",dim5="time", &
                                  units=var_now%units_out,start=[1,1,l,m,k],count=[gice%G%nx,gice%G%ny,1,1,1])
                end do 

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
    call nc_write_dim(file_ice,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
    call nc_write_dim(file_ice,"time", x=1979,dx=1,nx=33,units="years",calendar="360_day")
    
    call grid_write(gice,file_ice,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gice,icemask)
    call grid_allocate(gice,icevar)

    ! Write clim grid to file
    call nc_create(file_clim)
    call nc_write_dim(file_clim,"xc",   x=gclim%G%x,units="kilometers")
    call nc_write_dim(file_clim,"yc",   x=gclim%G%y,units="kilometers")
    call nc_write_dim(file_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
    call nc_write_dim(file_clim,"time", x=1979,dx=1,nx=33,units="years",calendar="360_day")
    
    call grid_write(gclim,file_clim,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gclim,climmask)
    call grid_allocate(gclim,climvar)


    ! Define MAR grid and input variable field
    call grid_init(gMAR,name="MAR-25KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-750.d0,dx=25.d0,nx=58,y0=-1200.d0,dy=25.d0,ny=108, &
                   lambda=-40.d0,phi=72.d0,alpha=7.5d0)
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
        call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value)
        climvar = missing_value 
        call map_field(mMAR_clim,var_now%nm_in,invar,climvar,climmask,var_now%method,100.d3, &
                      fill=.FALSE.,missing_value=missing_value)
        call nc_write(file_clim,var_now%nm_out,climvar,  dim1="xc",dim2="yc",units=var_now%units_out)
        icevar = missing_value 
        call map_field(mMAR_ice, var_now%nm_in,invar,icevar, icemask, var_now%method,100.d3, &
                       fill=.FALSE.,missing_value=missing_value)  
        call nc_write(file_ice, var_now%nm_out,icevar,   dim1="xc",dim2="yc",units=var_now%units_out)
    end do 

    stop 

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
                    call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value, &
                                  start=[1,1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1,1])
                else 
                    call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value, &
                             start=[1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1])
                end if 
                climvar = missing_value 
                call map_field(mMAR_clim,var_now%nm_in,invar,climvar,climmask,"shepard",100.d3, &
                               fill=.FALSE.,missing_value=missing_value)
                call nc_write(file_clim,var_now%nm_out,climvar,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=[1,1,m,k],count=[gclim%G%nx,gclim%G%ny,1,1])
                icevar = missing_value 
                call map_field(mMAR_ice, var_now%nm_in,invar,icevar, icemask, "shepard",100.d3, &
                               fill=.FALSE.,missing_value=missing_value)
                call nc_write(file_ice,var_now%nm_out,icevar,  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=[1,1,m,k],count=[gice%G%nx,gice%G%ny,1,1])
            end do 

        end do 
    end do 

    end if 

! ########################### 

    ! =========================================================
    !
    !       MAR (RCM) DATA - MARv3.3 downloaded from the ftp site:
    !       ftp://ftp.climato.be/fettweis/MARv3.3/Greenland
    !       Data is available on the Bamber et al. (2001) 5km grid
    !       ERA-40 + ERA-Interim combined datasets
    !
    ! =========================================================

! ########################### 
    if (.FALSE.) then 

    ! Define file names for input and output of global grids  
    file_ice       = "output/GRL-20KM_MARv3.3-15km-monthly-ERA-Interim_197901-201312.nc"
    file_clim      = "output/GRL-50KM_MARv3.3-15km-monthly-ERA-Interim_197901-201312.nc"

    ! Write ice grid to file
    call nc_create(file_ice)
    call nc_write_dim(file_ice,"xc",  x=gice%G%x,units="kilometers")
    call nc_write_dim(file_ice,"yc",  x=gice%G%y,units="kilometers")
    call nc_write_dim(file_ice,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
    call nc_write_dim(file_ice,"time", x=1958,dx=1,nx=56,units="years",calendar="360_day")
    
    call grid_write(gice,file_ice,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gice,icemask)
    call grid_allocate(gice,icevar)

    ! Write clim grid to file
    call nc_create(file_clim)
    call nc_write_dim(file_clim,"xc",   x=gclim%G%x,units="kilometers")
    call nc_write_dim(file_clim,"yc",   x=gclim%G%y,units="kilometers")
    call nc_write_dim(file_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
    call nc_write_dim(file_clim,"time", x=1958,dx=1,nx=56,units="years",calendar="360_day")
    
    call grid_write(gclim,file_clim,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gclim,climmask)
    call grid_allocate(gclim,climvar)

    ! Define MAR (Bamber et al. 2001) grid and input variable field
    call grid_init(gMAR,name="Bamber01-5KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-800.d0,dx=5.d0,nx=301,y0=-3400.d0,dy=5.d0,ny=561, &
                   lambda=-39.d0,phi=90.d0,alpha=7.5d0)
    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    call map_init(mMAR_ice, gMAR,gice, max_neighbors=30,lat_lim=3.d0,fldr="maps",load=.TRUE.)
    call map_init(mMAR_clim,gMAR,gclim,max_neighbors=30,lat_lim=3.d0,fldr="maps",load=.TRUE.)


    ! Define the variables to be mapped 
    file_mar = "/data/sicopolis/data/MARv3.3/Greenland/ERA_1958-2013_15km/"// &
               "MARv3.3-15km-monthly-ERA-Interim-2013.nc"

    allocate(mar_invariant(2))
    call def_var_info(mar_invariant(1),trim(file_mar),"MSK_MAR","mask",units="(0 - 2)",method="nn")
    call def_var_info(mar_invariant(2),trim(file_mar),"SRF_MAR","zs",units="m")

    file_mar = "/data/sicopolis/data/MARv3.3/Greenland/ERA_1958-2013_15km/"// &
               "MARv3.3-15km-monthly-"

    allocate(mar_surf(16))
    call def_var_info(mar_surf( 1),trim(file_mar),"SMB", "smb", units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
    call def_var_info(mar_surf( 2),trim(file_mar),"RU",  "ru",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
    call def_var_info(mar_surf( 3),trim(file_mar),"ME",  "me",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
    call def_var_info(mar_surf( 4),trim(file_mar),"ST",  "ts",  units="degrees Celcius")
    call def_var_info(mar_surf( 5),trim(file_mar),"TT",  "t3m", units="degrees Celcius")
    call def_var_info(mar_surf( 6),trim(file_mar),"SF",  "sf",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
    call def_var_info(mar_surf( 7),trim(file_mar),"RF",  "rf",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
    call def_var_info(mar_surf( 8),trim(file_mar),"SU",  "su",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
    call def_var_info(mar_surf( 9),trim(file_mar),"AL",  "al",  units="(0 - 1)")
    call def_var_info(mar_surf(10),trim(file_mar),"SWD", "swd", units="W m**-2")
    call def_var_info(mar_surf(11),trim(file_mar),"LWD", "lwd", units="W m**-2")
    call def_var_info(mar_surf(12),trim(file_mar),"SHF", "shf", units="W m**-2")
    call def_var_info(mar_surf(13),trim(file_mar),"LHF", "lhf", units="W m**-2")
    call def_var_info(mar_surf(14),trim(file_mar),"SP",  "sp",  units="hPa")
    call def_var_info(mar_surf(15),trim(file_mar),"UU",  "u",   units="m s**-1")
    call def_var_info(mar_surf(16),trim(file_mar),"VV",  "v",   units="m s**-1")

    ! (Re)Allocate the input grid variable
    call grid_allocate(gMAR,invar)

    ! ## INVARIANT FIELDS ##
    do i = 1, size(mar_invariant)
        var_now = mar_invariant(i) 
        call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value)
        climvar = missing_value 
        call map_field(mMAR_clim,var_now%nm_in,invar,climvar,climmask,var_now%method,100.d3, &
                      fill=.FALSE.,missing_value=missing_value)
        if (var_now%method .eq. "nn") then 
            call nc_write(file_clim,var_now%nm_out,nint(climvar),  dim1="xc",dim2="yc",units=var_now%units_out)
        else
            call nc_write(file_clim,var_now%nm_out,real(climvar),  dim1="xc",dim2="yc",units=var_now%units_out)
        end if 
        icevar = missing_value 
        call map_field(mMAR_ice, var_now%nm_in,invar,icevar, icemask, var_now%method,100.d3, &
                       fill=.FALSE.,missing_value=missing_value)  
        if (var_now%method .eq. "nn") then 
            call nc_write(file_ice, var_now%nm_out,nint(icevar),   dim1="xc",dim2="yc",units=var_now%units_out)
        else 
            call nc_write(file_ice, var_now%nm_out,real(icevar),   dim1="xc",dim2="yc",units=var_now%units_out)
        end if 
    end do 

    nyr = 2013-1958+1
    nm  = 12 
       
    do k = 1, nyr 

        year = 1957 + k 
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
                if (year .le. 1978) then   
                    write(var_now%filename,"(a,a,i4,a3,i4,a)") trim(file_mar),"ERA-40-",year,".nc"
                else
                    write(var_now%filename,"(a,a,i4,a3,i4,a)") trim(file_mar),"ERA-Interim-",year,".nc"
                end if
                call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value, &
                         start=[1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1])
                where (invar .ne. missing_value) invar = invar*var_now%conv 
                climvar = missing_value 
                call map_field(mMAR_clim,var_now%nm_in,invar,climvar,climmask,"shepard",100.d3, &
                               fill=.FALSE.,missing_value=missing_value)
                call nc_write(file_clim,var_now%nm_out,real(climvar),  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=[1,1,m,k],count=[gclim%G%nx,gclim%G%ny,1,1])
                icevar = missing_value 
                call map_field(mMAR_ice, var_now%nm_in,invar,icevar, icemask, "shepard",100.d3, &
                               fill=.FALSE.,missing_value=missing_value)
                call nc_write(file_ice,var_now%nm_out,real(icevar),  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              units=var_now%units_out,start=[1,1,m,k],count=[gice%G%nx,gice%G%ny,1,1])
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
    call nc_write_dim(file_ice,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
    call nc_write_dim(file_ice,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gice,file_ice,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gice,icemask)
    call grid_allocate(gice,icevar)

    ! Write clim grid to file
    call nc_create(file_clim)
    call nc_write_dim(file_clim,"xc",   x=gclim%G%x,units="kilometers")
    call nc_write_dim(file_clim,"yc",   x=gclim%G%y,units="kilometers")
    call nc_write_dim(file_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
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
        call nc_read(var_now%filename,var_now%nm_in,tmp,missing_value=missing_value)
        call thin(invar,tmp,by=10)
        if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
            where( invar .eq. missing_value ) invar = 0.d0 
        end if
        if (trim(var_now%nm_out) .eq. "zb") then 
            call fill(invar,missing_value=missing_value,fill_value=-1000.d0)
        end if 
        call map_field(mTOPO_clim,var_now%nm_in,invar,climvar,climmask,var_now%method,20.d3, &
                      fill=.TRUE.,missing_value=missing_value)
        call nc_write(file_clim,var_now%nm_out,climvar,  dim1="xc",dim2="yc",units=var_now%units_out)
        call map_field(mTOPO_ice, var_now%nm_in,invar,icevar, icemask, var_now%method,50.d3, &
                       fill=.TRUE.,missing_value=missing_value)
        call nc_write(file_ice, var_now%nm_out,icevar,   dim1="xc",dim2="yc",units=var_now%units_out)
    end do 

    end if 
! ########################### 

    ! =========================================================
    !
    !       OLR DATA
    !
    ! =========================================================

! ########################### 
    if (.FALSE.) then 

    ! Define file names for input and output of global grids  
    file_ice       = "output/GRL-20KM_OLR_clim1981-2010.nc"
    file_clim      = "output/GRL-50KM_OLR_clim1981-2010.nc"

    ! Write ice grid to file
    call nc_create(file_ice)
    call nc_write_dim(file_ice,"xc",  x=gice%G%x,units="kilometers")
    call nc_write_dim(file_ice,"yc",  x=gice%G%y,units="kilometers")
    call nc_write_dim(file_ice,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
!     call nc_write_dim(file_ice,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gice,file_ice,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gice,icemask)
    call grid_allocate(gice,icevar)

    ! Write clim grid to file
    call nc_create(file_clim)
    call nc_write_dim(file_clim,"xc",   x=gclim%G%x,units="kilometers")
    call nc_write_dim(file_clim,"yc",   x=gclim%G%y,units="kilometers")
    call nc_write_dim(file_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
!     call nc_write_dim(file_clim,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gclim,file_clim,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gclim,climmask)
    call grid_allocate(gclim,climvar)

    ! Define OLR grid and input variable field
    call grid_init(gOLR,name="OLR-2.5deg",mtype="latlon",units="degrees",lon180=.FALSE., &
                   x0=0.0d0,dx=2.5d0,nx=144,y0=-90.d0,dy=2.5d0,ny=73 )

    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    call map_init(mOLR_ice, gOLR,gice, max_neighbors=10,lat_lim=4d0,fldr="maps",load=.TRUE.)
    call map_init(mOLR_clim,gOLR,gclim,max_neighbors=10,lat_lim=4d0,fldr="maps",load=.TRUE.)

    ! Define the variables to be mapped 
    file_olr = "data/NOAA/OLR/olr.mon.ltm.nc"
    call def_var_info(olr_toa,trim(file_olr),"olr","olr",units="W m**-2",method="quadrant")

    ! (Re)Allocate the input grid variable
    call grid_allocate(gOLR,invar)

    ! Loop over months and map gridded variables
    nm = 12 
    do m = 1, nm 

        write(*,*)
        write(*,*) "= Month ",m, " ="
        write(*,*) 

        var_now = olr_toa
        call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value, &
                     start=[1,1,m],count=[gOLR%G%nx,gOLR%G%ny,1])
        call map_field(mOLR_clim,var_now%nm_in,invar,climvar,climmask,var_now%method, &
                      fill=.TRUE.,missing_value=missing_value)
        call nc_write(file_clim,var_now%nm_out,climvar,  dim1="xc",dim2="yc",dim3="month", &
                      units=var_now%units_out,start=[1,1,m],count=[gclim%G%nx,gclim%G%ny,1])
        call map_field(mOLR_ice, var_now%nm_in,invar,icevar, icemask, var_now%method, &
                       fill=.TRUE.,missing_value=missing_value)
        call nc_write(file_ice, var_now%nm_out,icevar,   dim1="xc",dim2="yc",dim3="month", &
                      units=var_now%units_out,start=[1,1,m],count=[gice%G%nx,gice%G%ny,1])
    
    end do 

    end if 
! ########################### 

    ! =========================================================
    !
    !       CERES DATA
    !
    ! =========================================================

! ########################### 
    if (.FALSE.) then 

    ! Define file names for input and output of global grids  
    file_ice       = "output/GRL-20KM_CERES_clim2001-2013.nc"
    file_clim      = "output/GRL-50KM_CERES_clim2001-2013.nc"

    ! Write ice grid to file
    call nc_create(file_ice)
    call nc_write_dim(file_ice,"xc",  x=gice%G%x,units="kilometers")
    call nc_write_dim(file_ice,"yc",  x=gice%G%y,units="kilometers")
    call nc_write_dim(file_ice,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
!     call nc_write_dim(file_ice,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gice,file_ice,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gice,icemask)
    call grid_allocate(gice,icevar)

    ! Write clim grid to file
    call nc_create(file_clim)
    call nc_write_dim(file_clim,"xc",   x=gclim%G%x,units="kilometers")
    call nc_write_dim(file_clim,"yc",   x=gclim%G%y,units="kilometers")
    call nc_write_dim(file_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
!     call nc_write_dim(file_clim,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
    
    call grid_write(gclim,file_clim,xnm="xc",ynm="yc",create=.FALSE.)
    call grid_allocate(gclim,climmask)
    call grid_allocate(gclim,climvar)

    ! Define CERES grid and input variable field
    call grid_init(gCERES,name="CERES-1deg",mtype="latlon",units="degrees",lon180=.FALSE., &
                   x0=0.5d0,dx=1d0,nx=360,y0=-90.d0,dy=1d0,ny=180 )

    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    call map_init(mCERES_ice, gCERES,gice, max_neighbors=10,lat_lim=4d0,fldr="maps",load=.TRUE.)
    call map_init(mCERES_clim,gCERES,gclim,max_neighbors=10,lat_lim=4d0,fldr="maps",load=.TRUE.)

    ! Define the variables to be mapped 
    allocate(ceres(10))
    file_CERES = "data/CERES/CERES_EBAF-TOA_Ed2.8_Subset_CLIM01-CLIM12.nc"
    call def_var_info(ceres( 1),trim(file_ceres),"toa_sw_all_clim","toa_sw_all",units="W m**-2",method="radius")
    call def_var_info(ceres( 2),trim(file_ceres),"toa_sw_clr_clim","toa_sw_clr",units="W m**-2",method="radius")
    call def_var_info(ceres( 3),trim(file_ceres),"toa_lw_all_clim","toa_lw_all",units="W m**-2",method="radius")
    call def_var_info(ceres( 4),trim(file_ceres),"toa_lw_clr_clim","toa_lw_clr",units="W m**-2",method="radius")
    call def_var_info(ceres( 5),trim(file_ceres),"toa_net_all_clim","toa_net_all",units="W m**-2",method="radius")
    call def_var_info(ceres( 6),trim(file_ceres),"toa_net_clr_clim","toa_net_clr",units="W m**-2",method="radius")
    call def_var_info(ceres( 7),trim(file_ceres),"toa_cre_sw_clim","toa_cre_sw",units="W m**-2",method="radius")
    call def_var_info(ceres( 8),trim(file_ceres),"toa_cre_lw_clim","toa_cre_lw",units="W m**-2",method="radius")
    call def_var_info(ceres( 9),trim(file_ceres),"toa_cre_net_clim","toa_cre_net",units="W m**-2",method="radius")
    call def_var_info(ceres(10),trim(file_ceres),"solar_clim","solar",units="W m**-2",method="radius")

    ! (Re)Allocate the input grid variable
    call grid_allocate(gCERES,invar)

    ! Loop over months and map gridded variables
    nm = 12 
    do m = 1, nm 

        write(*,*)
        write(*,*) "= Month ",m, " ="
        write(*,*) 

        do i = 1, size(ceres)
            var_now = ceres(i)
            call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value, &
                         start=[1,1,m],count=[gCERES%G%nx,gCERES%G%ny,1])
            call map_field(mCERES_clim,var_now%nm_in,invar,climvar,climmask,var_now%method, &
                          fill=.TRUE.,missing_value=missing_value)
            call nc_write(file_clim,var_now%nm_out,climvar,  dim1="xc",dim2="yc",dim3="month", &
                          units=var_now%units_out,start=[1,1,m],count=[gclim%G%nx,gclim%G%ny,1])
            call map_field(mCERES_ice, var_now%nm_in,invar,icevar, icemask, var_now%method, &
                           fill=.TRUE.,missing_value=missing_value)
            call nc_write(file_ice, var_now%nm_out,icevar,   dim1="xc",dim2="yc",dim3="month", &
                          units=var_now%units_out,start=[1,1,m],count=[gice%G%nx,gice%G%ny,1])
        end do 
    end do 

    end if 
! ########################### 

end program gentopo

