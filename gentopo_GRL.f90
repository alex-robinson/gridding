
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
    if (.FALSE.) then 

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
    !       CERES DATA
    !
    ! =========================================================
    if (.TRUE.) then 

        outfldr = "output/Greenland"

        ! Map to the grids of interest from global grid 
        call CERES_to_grid(outfldr, g50KM,"Global",max_neighbors=8,lat_lim=2.d0)
        call CERES_to_grid(outfldr, g20KM,"Global",max_neighbors=8,lat_lim=2.d0)
        call CERES_to_grid(outfldr,g20KMb,"Global",max_neighbors=8,lat_lim=2.d0)
        call CERES_to_grid(outfldr, g10KM,"Global",max_neighbors=8,lat_lim=2.d0)

    end if 

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gentopo

