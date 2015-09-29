
program gen_gridded

    use ncio 
    use coordinates
    use control 

    ! Datasets that can be gridded by this program (alphabetical)
    use AN1CRUST
    use bedmap2
    use CERES
    use climber2
    use climber3a 
    use davini2015
    use ECMWF 
    use ETOPO 
    use GeothermalHeatFlux 
    use MAR 
    use NasaBasins 
    use RACMO2 
    use Rignot13_BasalMelt  
    use sediments 
    use stratigraphy
    use topo_reconstructions 
    use topographies_grl 

    implicit none

    type(grid_class)   :: grid
    character(len=256) :: domain, grid_name, outfldr 
    character(len=256) :: path, subfldr  
    real(8) :: sigma1, sigma2 

    write(*,*) 

    ! =========================================================
    !
    ! USER DEFINITIONS
    !
    ! =========================================================
    
    domain    = "North"
    grid_name = "NH-40KM"
    outfldr   = "output/"//trim(domain)//"/"//trim(grid_name)

    ! =========================================================
    !
    ! OUTPUT GRID DEFINITION
    !
    ! =========================================================

    call domain_definition(grid,grid_name)
    
    ! =========================================================
    !
    ! DATASET TO GRID CALCULATIONS
    !
    ! =========================================================

    ! == Global datasets - applicable to all domains ==

!     call CERES_to_grid(outfldr,     grid,"Global",max_neighbors=4,lat_lim=2.d0)
!     call etopo1_to_grid(outfldr,    grid,"Global",max_neighbors=1,lat_lim=1.d0)
!     call ICE6GC_to_grid(outfldr,    grid,"Global",max_neighbors=4,lat_lim=2.d0)
!     call ICE5G_to_grid(outfldr,     grid,"Global",max_neighbors=4,lat_lim=2.d0)
!     call sedLaske_to_grid(outfldr,  grid,"Global",max_neighbors=4,lat_lim=2.d0)
!     call ghfDavies_to_grid(outfldr, grid,"Global",max_neighbors=4,lat_lim=2.d0)
!     call ghfShapiro_to_grid(outfldr,grid,"Global",max_neighbors=4,lat_lim=2.d0)
    
!     call ecmwf_to_grid(outfldr,grid,sigma=30.d0,max_neighbors=1,lat_lim=2.d0)
!     call ecmwf_to_grid(outfldr,grid,clim_range=[1981,2010])

!     call ecmwf_ocn_to_grid(outfldr,grid,sigma=30.d0,max_neighbors=4,lat_lim=2.d0)
!     call ecmwf_ocn_to_grid(outfldr,grid,clim_range=[1981,2010])

!     call ecmwf40_to_grid(outfldr,grid,sigma=100.d0,max_neighbors=1,lat_lim=2.d0)
!     call ecmwf40_to_grid(outfldr,grid,clim_range=[1958,2001])

    path = "/data/sicopolis/data/CLIMBER3a/Montoya2008/"
    sigma1 = 250.d0 
    sigma2 = 100.d0 
    write(subfldr,"(a15,i3,a2)") "Montoya2008_sig", int(sigma1), "km"
    call system("mkdir -p "//trim(outfldr)//"/"//trim(subfldr))

    call climber3a_atm_to_grid(outfldr,subfldr,grid,domain="lgm_1p7strong", &
                               path_in=path,sigma=sigma1,max_neighbors=10,lat_lim=5.d0)
    call climber3a_atm_to_grid(outfldr,subfldr,grid,domain="lgm_1p7weak", &
                               path_in=path,sigma=sigma1,max_neighbors=10,lat_lim=5.d0)
    call climber3a_atm_to_grid(outfldr,subfldr,grid,domain="present", &
                               path_in=path,sigma=sigma1,max_neighbors=10,lat_lim=5.d0)
    
!     write(subfldr,"(a15,i3,a2)") "Montoya2008_sig", int(sigma2), "km"
!     call system("mkdir -p "//trim(outfldr)//"/"//trim(subfldr))
!     call climber3a_ocn_to_grid(outfldr,subfldr,grid,domain="lgm_1p7strong_ocean", &
!                                path_in=path,sigma=sigma2,max_neighbors=10,lat_lim=5.d0)
!     call climber3a_ocn_to_grid(outfldr,subfldr,grid,domain="lgm_1p7weak_ocean", &
!                                path_in=path,sigma=sigma2,max_neighbors=10,lat_lim=5.d0)
!     call climber3a_ocn_to_grid(outfldr,subfldr,grid,domain="present_ocean", &
!                                path_in=path,sigma=sigma2,max_neighbors=10,lat_lim=5.d0)

!     path = "/data/sicopolis/data/CLIMBER3a/Montoya2008hi/"
!     sigma1 = 10.d0 
!     call climber3a_atm_to_grid(outfldr,"Montoya2008hi",grid,domain="lgm_1p7strong", &
!                                path_in=path,sigma=sigma1,max_neighbors=10,lat_lim=5.d0)
!     call climber3a_atm_to_grid(outfldr,"Montoya2008hi",grid,domain="lgm_1p7weak", &
!                                path_in=path,sigma=sigma1,max_neighbors=10,lat_lim=5.d0)
!     call climber3a_atm_to_grid(outfldr,"Montoya2008hi",grid,domain="present", &
!                                path_in=path,sigma=sigma1,max_neighbors=10,lat_lim=5.d0)
    
!     path = "/data/sicopolis/data/CLIMBER3a/Montoya2008_jorgebilin/"
!     call climber3a_jorge_to_grid(outfldr,"Montoya2008_jorgebilin",grid,domain="lgm_1p7weak", &
!                                path_in=path,sigma=sigma1,max_neighbors=10,lat_lim=5.d0,load_topo=.TRUE.)
!     call climber3a_jorge_to_grid(outfldr,"Montoya2008_jorgebilin",grid,domain="present", &
!                                path_in=path,sigma=sigma1,max_neighbors=10,lat_lim=5.d0,load_topo=.TRUE.)
!     call climber3a_jorge_to_grid(outfldr,"Montoya2008_jorgebilin",grid,domain="climatology", &
!                                path_in=path,sigma=sigma1,max_neighbors=10,lat_lim=5.d0,load_topo=.FALSE.)
    
!     path = "data/climber_data/NCO2_nc/"
!     call climber2_atm_to_grid(outfldr,"Ganopolski2011",grid,sim="860ka", &
!                               path_in=path,max_neighbors=4,lat_lim=20.d0)

    if (trim(domain) .eq. "Antarctica") then 
        ! == Antarctica only datasets ==

!         call An15litho_to_grid(outfldr, grid,"Antarctica", max_neighbors=5, lat_lim=1.0d0)
!         call bedmap2_to_grid(outfldr,   grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0)
!         call bedmap2vel_to_grid(outfldr,grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0)
!         call bedmap2acc_to_grid(outfldr,grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0)
!         call ghfMaule_to_grid(outfldr,  grid,"Antarctica",max_neighbors=4,lat_lim=2.d0)
!         call nasaBasins_to_grid(outfldr,grid,"Antarctica")
!         call Rignot13_BasalMelt_to_grid(outfldr,grid,"Antarctica",max_neighbors=10,lat_lim=1.d0)
        
!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",max_neighbors=20,lat_lim=0.5d0)
!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2000,2010])
!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2001,2030])
!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2071,2100])

!         call RACMO23_to_grid( outfldr, grid, "ANT27",max_neighbors=20,lat_lim=0.5d0)
!         call RACMO23_to_grid( outfldr, grid, "ANT27",clim_range=[1981,2010])

        ! Note: Antartica-c20 doesn't work because some files only contain 239 months of
        !      data while they should all have 240 months (1980-1999)
        !      This dataset is not used for now...
        !call RACMO2rot_to_grid( outfldr, grid, "Antarctica-c20",max_neighbors=20,lat_lim=0.5d0)
        !call RACMO2rot_to_grid( outfldr, grid, "Antarctica-c20",clim_range=[1980,1999])
    
    end if 

    if (trim(domain) .eq. "Greenland") then 
        ! == Greenland only datasets ==

!         call Bamber13_to_grid(outfldr,grid,"Greenland",max_neighbors=10,lat_lim=2.d0)   
!         call ghfMaule_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=2.d0)
!         call Morlighem14_to_grid(outfldr,grid,"Greenland",max_neighbors=20,lat_lim=1.d0)
!         call MARv35_to_grid(outfldr,grid,"Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
!         call MARv35_to_grid(outfldr,grid,"Greenland-ERA",clim_range=[1981,2010])
!         call nasaBasins_to_grid(outfldr,grid,"Antarctica")
!         call LGMsimpson_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=1.d0)
        
!         path = "data/Davini_GreenlandAMOC/"
!         call davini2015_to_grid(outfldr,"Davini2015",grid,domain="control", &
!                                 path_in=path,sigma=40.d0,max_neighbors=4,lat_lim=5.d0)
!         call davini2015_to_grid(outfldr,"Davini2015",grid,domain="bedrock", &
!                                 path_in=path,sigma=40.d0,max_neighbors=4,lat_lim=5.d0)
    
    end if 


    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gen_gridded

