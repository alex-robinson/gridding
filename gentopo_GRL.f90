
program gentopo

    use ncio 
    use coordinates
    use vargrid 
    use gridding_datasets

    implicit none

    type(grid_class) :: g50KM, g25KM, g20KM, g20KMb, g10KM
    character(len=256) :: file_50KM, file_25KM, file_20KM, file_20KMb, file_10KM
    character(len=256) :: outfldr 

    type(grid_class)   :: gice, gclim
    character(len=256) :: file_ice, file_clim
    
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

    ! ## Define clim grid and output variable field ##
    call grid_init(g50KM,name="GRL-50KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   dx=50.d0,nx=37,dy=50.d0,ny=61, &
                   lambda=-40.d0,phi=72.d0,alpha=7.5d0)

    ! ## Define ice grid and output variable field ##
    call grid_init(g20KM,name="GRL-20KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   dx=20.d0,nx=90,dy=20.d0,ny=150, &
                   lambda=-40.d0,phi=72.d0,alpha=7.5d0)

    ! Define Bamber et al. 2001 20KM grid and input variable field
    call grid_init(g20KMb,name="Bamber01-20KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-800.d0,dx=20.d0,nx=76,y0=-3400.d0,dy=20.d0,ny=141, &
                   lambda=-39.d0,phi=90.d0,alpha=7.5d0)

    ! ## Define ice grid and output variable field ##
    call grid_init(g10KM,name="GRL-10KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   dx=10.d0,nx=180,dy=10.d0,ny=300, &
                   lambda=-40.d0,phi=72.d0,alpha=7.5d0)

    ! For compilability
    gice  = g20KM 
    gclim = g50KM 

    write(*,*) 

    ! =========================================================
    !
    !       TOPO DATA
    !
    ! =========================================================
    if (.TRUE.) then 

        outfldr = "output/Greenland"

        ! Map to the grids of interest from Bamber et al (2013) grid
        call Bamber13_to_grid(outfldr, g50KM, "Greenland",max_neighbors=20,lat_lim=2.d0)
        call Bamber13_to_grid(outfldr, g20KM, "Greenland",max_neighbors=15,lat_lim=2.d0)
        call Bamber13_to_grid(outfldr,g20KMb, "Greenland",max_neighbors=10,lat_lim=2.d0)
        call Bamber13_to_grid(outfldr, g10KM, "Greenland",max_neighbors=10,lat_lim=2.d0)

    end if 



    ! =========================================================
    !
    !       ECMWF DATA
    !
    ! =========================================================  
    if (.FALSE.) then 

        outfldr = "output/Greenland"

        ! Map to the grids of interest from 0.75 degree ECMWF dataset
        call ecmwf_to_grid(outfldr, g50KM, "GRL075",max_neighbors=8,lat_lim=2.d0)
        call ecmwf_to_grid(outfldr, g20KM, "GRL075",max_neighbors=8,lat_lim=2.d0)
        call ecmwf_to_grid(outfldr,g20KMb, "GRL075",max_neighbors=8,lat_lim=2.d0)
        call ecmwf_to_grid(outfldr, g10KM, "GRL075",max_neighbors=8,lat_lim=2.d0)

    end if 



    ! =========================================================
    !
    !       MAR (RCM) DATA - MARv3.3 downloaded from the ftp site:
    !       ftp://ftp.climato.be/fettweis/MARv3.3/Greenland
    !       Data is available on the Bamber et al. (2001) 5km grid
    !       ERA-40 + ERA-Interim combined datasets
    !
    ! =========================================================
    if (.FALSE.) then 

        outfldr = "output/Greenland"

        ! Map to the grids of interest from 5 km Greenland grid 
        call MARv33_to_grid(outfldr, g50KM, "Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
        call MARv33_to_grid(outfldr, g20KM, "Greenland-ERA",max_neighbors=15,lat_lim=2.d0)
        call MARv33_to_grid(outfldr,g20KMb, "Greenland-ERA",max_neighbors=15,lat_lim=2.d0)
        call MARv33_to_grid(outfldr, g10KM, "Greenland-ERA",max_neighbors=10,lat_lim=2.d0)

    end if 

    ! =========================================================
    !
    !       MAR (RCM) DATA - MARv3.3 downloaded from the ftp site:
    !       ftp://ftp.climato.be/fettweis/MARv3.3/Greenland
    !       Data is available on the Bamber et al. (2001) 5km grid
    !       MIROC5 histo+rcp85 combined datasets
    !
    ! =========================================================
    if (.FALSE.) then 

        outfldr = "output/Greenland"

        ! Map to the grids of interest from 5 km Greenland grid 
        call MARv33_to_grid(outfldr, g50KM, "Greenland-MIROC5-RCP85",max_neighbors=20,lat_lim=2.d0)
        call MARv33_to_grid(outfldr, g20KM, "Greenland-MIROC5-RCP85",max_neighbors=15,lat_lim=2.d0)
        call MARv33_to_grid(outfldr,g20KMb, "Greenland-MIROC5-RCP85",max_neighbors=15,lat_lim=2.d0)
        call MARv33_to_grid(outfldr, g10KM, "Greenland-MIROC5-RCP85",max_neighbors=10,lat_lim=2.d0)

    end if 

    ! =========================================================
    !
    !       MAR (RCM) DATA - MARv3.2 original data passed by Xavier
    !
    ! =========================================================
    if (.FALSE.) then 

        outfldr = "output/Greenland"

        ! Map to the grids of interest from 5 km Greenland grid 
        call MARv32_to_grid(outfldr, g50KM, "Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
        call MARv32_to_grid(outfldr, g20KM, "Greenland-ERA",max_neighbors=15,lat_lim=2.d0)
        call MARv32_to_grid(outfldr,g20KMb, "Greenland-ERA",max_neighbors=15,lat_lim=2.d0)
        call MARv32_to_grid(outfldr, g10KM, "Greenland-ERA",max_neighbors=10,lat_lim=2.d0)

    end if 



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
    

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gentopo

