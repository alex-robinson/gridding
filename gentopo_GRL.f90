
program gentopo

    use ncio 
    use coordinates

    use CERES 
    use climber2
    use climber3a
    use ECMWF
    use ETOPO 
    use GeothermalHeatFlux
    use MAR 
    use NasaBasins
    use sediments 
    use stratigraphy
    use topo_reconstructions 
    use topographies_grl 

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
 
    gridname = "GRL-40KM"
    outfldr  = "output/Greenland/"//trim(gridname)

    ! =========================================================
    !
    ! OUTPUT GRID DEFINITION
    !
    ! =========================================================

    select case(trim(gridname))
        case("GRL-40KM")
            call grid_init(grid,name="GRL-40KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=40.d0,nx=45,dy=40.d0,ny=75, &
                           lambda=-40.d0,phi=72.d0,alpha=8.4d0)

        case("GRL-20KM")
            call grid_init(grid,name="GRL-20KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=20.d0,nx=90,dy=20.d0,ny=150, &
                           lambda=-40.d0,phi=72.d0,alpha=8.4d0)

        case("GRL-10KM")
            call grid_init(grid,name="GRL-10KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=10.d0,nx=180,dy=10.d0,ny=300, &
                           lambda=-40.d0,phi=72.d0,alpha=8.4d0)

        case("Bamber01-20KM")
            call grid_init(grid,name="Bamber01-20KM",mtype="polar_stereographic",units="kilometers", &
                           lon180=.TRUE.,x0=-800.d0,dx=20.d0,nx=76,y0=-3400.d0,dy=20.d0,ny=141, &
                           lambda=-39.d0,phi=90.d0,alpha=7.5d0)

        case("Bamber01-10KM")
            call grid_init(grid,name="Bamber01-10KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,x0=-800.d0,dx=10.d0,nx=151,y0=-3400.d0,dy=10.d0,ny=281, &
                           lambda=-39.d0,phi=90.d0,alpha=7.5d0)

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

    call Morlighem14_to_grid(outfldr,grid,"Greenland",max_neighbors=20,lat_lim=1.d0)


    if (.FALSE.) then 
    call Bamber13_to_grid(outfldr,grid,"Greenland",max_neighbors=10,lat_lim=2.d0)

    call CERES_to_grid(outfldr,grid,"Global",max_neighbors=4,lat_lim=2.d0)

    call ecmwf_to_grid(outfldr, grid,"Global",sigma=30.d0,max_neighbors=4,lat_lim=2.d0)
    call ecmwf_to_grid( outfldr,grid,"Global",clim_range=[1981,2010])

    call etopo1_to_grid(outfldr,grid,"Greenland",max_neighbors=1,lat_lim=1.d0)
  
    call MARv35_to_grid(outfldr,grid,"Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
    call MARv35_to_grid(outfldr,grid,"Greenland-ERA",clim_range=[1981,2010])

    call nasaBasins_to_grid(outfldr,grid,"Greenland")

    call sedLaske_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=2.d0)
    call ghfMaule_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=2.d0)
    call ghfDavies_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=2.d0)
    call ghfShapiro_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=2.d0)
    
    ! Paleo topography 
    call ICE6GC_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=2.d0)
    call ICE5G_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=2.d0)

    call LGMsimpson_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=1.d0)
    
    ! CLIMBER-3alpha
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

    ! CLIMBER2
    path = "data/climber_data/NCO2_nc/"
    call climber2_atm_to_grid(outfldr,"Ganopolski2011",grid,sim="860ka", &
                              path_in=path,max_neighbors=4,lat_lim=20.d0)
    
    end if 
    

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gentopo



! ## NOT SUPPORTED ANYMORE:

!     call MARv33_to_grid(outfldr,grid,"Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
!     call MARv33_to_grid(outfldr,grid,"Greenland-ERA",clim_range=[1981,2010])

!     call MARv33_to_grid(outfldr,grid,"Greenland-MIROC5-RCP85",max_neighbors=20,lat_lim=2.d0)
!     call MARv33_to_grid(outfldr,grid,"Greenland-MIROC5-RCP85",clim_range=[1981,2010])
!     call MARv33_to_grid(outfldr,grid,"Greenland-MIROC5-RCP85",clim_range=[2071,2100])

!     call MARv32_to_grid(  outfldr, grid, "Greenland-ERA",max_neighbors=20,lat_lim=2.d0)

