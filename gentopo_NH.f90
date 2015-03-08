
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

!     call etopo1_to_grid(outfldr,grid,"NH",max_neighbors=4,lat_lim=2.d0)

!     call CERES_to_grid(outfldr,grid,"Global",max_neighbors=8,lat_lim=2.d0)

!     call ecmwf_to_grid(outfldr, grid,"NH",max_neighbors=8,lat_lim=2.d0)
!     call ecmwf_to_grid( outfldr,grid,"NH",clim_range=[1981,2010])

!     call sedLaske_to_grid(outfldr,grid,"NH",max_neighbors=10,lat_lim=2.d0)
    
    call ghfDavies_to_grid(outfldr,grid,"NH",max_neighbors=10,lat_lim=2.d0)
    
    ! CLIMBER-3alpha
!     call climber3a_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7strong",max_neighbors=10,lat_lim=5.d0)
!     call climber3a_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7weak",max_neighbors=10,lat_lim=5.d0)
!     call climber3a_to_grid(outfldr,"Montoya2008",grid,domain="present",max_neighbors=10,lat_lim=5.d0)
    
    ! Paleo topography 
!     call ICE6GC_to_grid(outfldr,grid,"NH",max_neighbors=4,lat_lim=2.d0)
!     call ICE5G_to_grid(outfldr,grid,"NH",max_neighbors=4,lat_lim=2.d0)

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gentopo

