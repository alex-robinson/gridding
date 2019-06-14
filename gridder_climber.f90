
program gridder

    use coord 
    use control 
    
    ! Datasets that can be gridded by this program (alphabetical)
    use climber2
    use climber3a 

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
    grid_name = "NH-32KM"
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
    
    write(subfldr,"(a15,i3,a2)") "Montoya2008_sig", int(sigma2), "km"
    call system("mkdir -p "//trim(outfldr)//"/"//trim(subfldr))
    call climber3a_ocn_to_grid(outfldr,subfldr,grid,domain="lgm_1p7strong_ocean", &
                               path_in=path,sigma=sigma2,max_neighbors=10,lat_lim=5.d0)
    call climber3a_ocn_to_grid(outfldr,subfldr,grid,domain="lgm_1p7weak_ocean", &
                               path_in=path,sigma=sigma2,max_neighbors=10,lat_lim=5.d0)
    call climber3a_ocn_to_grid(outfldr,subfldr,grid,domain="present_ocean", &
                               path_in=path,sigma=sigma2,max_neighbors=10,lat_lim=5.d0)

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

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gridder

