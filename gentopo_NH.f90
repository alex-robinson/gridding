
program gentopo

    use ncio 
    use coordinates

    use CERES 
    use climber3a
    use ECMWF 
    use ETOPO 
    use GeothermalHeatFlux
    use sediments  
    use topo_reconstructions 
    
    implicit none

    type(grid_class)   :: grid
    character(len=256) :: gridname, outfldr 
    character(len=256) :: path 

    write(*,*) 

    ! =========================================================
    !
    ! USER DEFINITIONS
    !
    ! =========================================================
 
    gridname = "NH-40KM"
    outfldr  = "output/NH/"//trim(gridname)

    ! =========================================================
    !
    ! OUTPUT GRID DEFINITION
    !
    ! =========================================================

    select case(trim(gridname))
        case("NH-40KM")
            call grid_init(grid,name="NH-40KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=40.d0,nx=224,dy=40.d0,ny=208, &
                           lambda=-53.d0,phi=78.d0,alpha=32.7d0)

        case DEFAULT
            write(*,*) "gentopo:: error: grid name not recognized: "//trim(gridname)
            stop 

    end select

    call grid_write(grid,"output/"//trim(gridname)//".nc",xnm="xc",ynm="yc",create=.TRUE.) 

    ! =========================================================
    !
    ! DATASET TO GRID CALCULATIONS
    !
    ! =========================================================

    call CERES_to_grid(outfldr,grid,"Global",max_neighbors=4,lat_lim=2.d0)

    call ecmwf_to_grid(outfldr,grid,sigma=30.d0,max_neighbors=1,lat_lim=2.d0)
    call ecmwf_to_grid(outfldr,grid,clim_range=[1981,2010])

    call etopo1_to_grid(outfldr,grid,"NH",max_neighbors=1,lat_lim=1.d0)
  
    call ghfDavies_to_grid(outfldr,grid,"NH",max_neighbors=4,lat_lim=2.d0)
    call ghfShapiro_to_grid(outfldr,grid,"NH",max_neighbors=4,lat_lim=2.d0)

!     Paleo topography 
    call ICE6GC_to_grid(outfldr,grid,"NH",max_neighbors=4,lat_lim=2.d0)
    call ICE5G_to_grid(outfldr,grid,"NH",max_neighbors=4,lat_lim=2.d0)

    call sedLaske_to_grid(outfldr,grid,"NH",max_neighbors=4,lat_lim=2.d0)

!     ### MODEL DATA ### 

!     CLIMBER-3alpha
    path = "/data/sicopolis/data/CLIMBER3a/Montoya2008/"
    call climber3a_atm_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7strong", &
                               path_in=path,sigma=250.d0,max_neighbors=10,lat_lim=5.d0)
    call climber3a_atm_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7weak", &
                               path_in=path,sigma=250.d0,max_neighbors=10,lat_lim=5.d0)
    call climber3a_atm_to_grid(outfldr,"Montoya2008",grid,domain="present", &
                               path_in=path,sigma=250.d0,max_neighbors=10,lat_lim=5.d0)
    
    path = "/data/sicopolis/data/CLIMBER3a/Montoya2008/"
    call climber3a_ocn_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7strong_ocean", &
                               path_in=path,sigma=100.d0,max_neighbors=10,lat_lim=5.d0)
    call climber3a_ocn_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7weak_ocean", &
                               path_in=path,sigma=100.d0,max_neighbors=10,lat_lim=5.d0)
    call climber3a_ocn_to_grid(outfldr,"Montoya2008",grid,domain="present_ocean", &
                               path_in=path,sigma=100.d0,max_neighbors=10,lat_lim=5.d0)
    

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gentopo

