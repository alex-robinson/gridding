
program gentopo

    use ncio 
    use coordinates
    
    use AN1CRUST
    use bedmap2
    use CERES
    use climber3a 
    use ECMWF 
    use ETOPO 
    use GeothermalHeatFlux 
    use NasaBasins 
    use RACMO2 
    use Rignot13_BasalMelt  
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
 
    gridname = "ANT-40KM"
    outfldr  = "output/Antarctica/"//trim(gridname)

    ! =========================================================
    !
    ! OUTPUT GRID DEFINITION
    !
    ! =========================================================

    select case(trim(gridname))
        case("ANT-40KM")
            call grid_init(grid,name="ANT-40KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=40.d0,nx=141,dy=40.d0,ny=141, &
                           lambda=0.d0,phi=-90.d0,alpha=19.0d0)

        case("ANT-20KM")
            call grid_init(grid,name="ANT-20KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=20.d0,nx=281,dy=20.d0,ny=281, &
                           lambda=0.d0,phi=-90.d0,alpha=19.0d0)

        case DEFAULT
            write(*,*) "gentopo:: error: grid name not recognized: "//trim(gridname)
            stop 

    end select

    ! =========================================================
    !
    ! DATASET TO GRID CALCULATIONS
    !
    ! =========================================================

!     call An15litho_to_grid(outfldr,grid,"Antarctica",max_neighbors=5,lat_lim=1.0d0)

!     call bedmap2_to_grid(outfldr,   grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0)
!     call bedmap2vel_to_grid(outfldr,grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0)
!     call bedmap2acc_to_grid(outfldr,grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0)

!     call CERES_to_grid(outfldr,grid,"Global",max_neighbors=4, lat_lim=2.d0)
    
!     call ecmwf_to_grid(outfldr,grid,"NH",sigma=30.d0,max_neighbors=1,lat_lim=2.d0)
!     call ecmwf_to_grid(outfldr,grid,"NH",clim_range=[1981,2010])

!     call etopo1_to_grid(outfldr,grid,"Antarctica",max_neighbors=1,lat_lim=1.d0)
  
!     ! Paleo topography 
!     call ICE6GC_to_grid(outfldr,grid,"Antarctica",max_neighbors=4,lat_lim=2.d0)
!     call ICE5G_to_grid(outfldr,grid,"Antarctica",max_neighbors=4,lat_lim=2.d0)

!     call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",max_neighbors=20,lat_lim=0.5d0)
!     call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2000,2010])
!     call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2001,2030])
!     call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2071,2100])
    
!     call RACMO23_to_grid( outfldr, grid, "ANT27",max_neighbors=20,lat_lim=0.5d0)
!     call RACMO23_to_grid( outfldr, grid, "ANT27",clim_range=[1981,2010])

    ! Note: Antartica-c20 doesn't work because some files only contain 239 months of
    !      data while they should all have 240 months (1980-1999)
    !      This dataset is not used for now...
!     call RACMO2rot_to_grid( outfldr, grid, "Antarctica-c20",max_neighbors=20,lat_lim=0.5d0)
!     call RACMO2rot_to_grid( outfldr, grid, "Antarctica-c20",clim_range=[1980,1999])
    
    ! Rignot basal melting data
!     call Rignot13_BasalMelt_to_grid(outfldr,grid,"Antarctica",max_neighbors=10, lat_lim=1.d0)
        
    ! NASA drainage basins 
!     call nasaBasins_to_grid(outfldr,grid,"Antarctica")
    
!     call sedLaske_to_grid(outfldr,grid,"Antarctica",max_neighbors=4,lat_lim=2.d0)
!     call ghfMaule_to_grid(outfldr,grid,"Antarctica",max_neighbors=4,lat_lim=2.d0)
!     call ghfDavies_to_grid(outfldr,grid,"Antarctica",max_neighbors=4,lat_lim=2.d0)
!     call ghfShapiro_to_grid(outfldr,grid,"Antarctica",max_neighbors=4,lat_lim=2.d0)
    
    ! CLIMBER-3alpha
!     path = "/data/sicopolis/data/CLIMBER3a/Montoya2008/"
!     call climber3a_atm_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7strong", &
!                                path_in=path,sigma=250.d0,max_neighbors=10,lat_lim=5.d0)
!     call climber3a_atm_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7weak", &
!                                path_in=path,sigma=250.d0,max_neighbors=10,lat_lim=5.d0)
!     call climber3a_atm_to_grid(outfldr,"Montoya2008",grid,domain="present", &
!                                path_in=path,sigma=250.d0,max_neighbors=10,lat_lim=5.d0)
    
!     path = "/data/sicopolis/data/CLIMBER3a/Montoya2008/"
!     call climber3a_ocn_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7strong_ocean", &
!                                path_in=path,sigma=100.d0,max_neighbors=10,lat_lim=5.d0)
!     call climber3a_ocn_to_grid(outfldr,"Montoya2008",grid,domain="lgm_1p7weak_ocean", &
!                                path_in=path,sigma=100.d0,max_neighbors=10,lat_lim=5.d0)
!     call climber3a_ocn_to_grid(outfldr,"Montoya2008",grid,domain="present_ocean", &
!                                path_in=path,sigma=100.d0,max_neighbors=10,lat_lim=5.d0)
    
    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gentopo

