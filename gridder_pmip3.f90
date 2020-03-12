program gridder

    use coord 
    use control 
    
    ! Datasets that can be gridded by this program (alphabetical)
    use pmip3 

    implicit none

    type(grid_class)   :: grid
    character(len=256) :: domain, grid_name, outfldr 
    character(len=256) :: path, subfldr, sigma_str  
    real(8) :: sigma

    write(*,*) 

    ! =========================================================
    !
    ! USER DEFINITIONS
    !
    ! =========================================================
    
    domain    = "Greenland"
    grid_name = "GRL-32KM"
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

    path = "data/PMIP3/"
    sigma = 50.0d0 
    
    write(sigma_str,*) int(sigma)
    sigma_str = trim(adjustl(sigma_str))

    write(subfldr,"(a15,a,a2)") "PMIP3_sig", trim(sigma_str), "km"
    subfldr = trim(adjustl(subfldr))
    outfldr = trim(outfldr)//"/"//trim(subfldr)
    
    call system("mkdir -p "//trim(outfldr))

!     call LGM_dzs_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!     call LGM_extensions_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                                path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
    
    call PMIP3_to_grid(outfldr,grid,domain,pmip_case="CCSM4-PiControl", &
                              path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
    call PMIP3_to_grid(outfldr,grid,domain,pmip_case="CCSM4-LGM", &
                              path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!     call CCSM4_PD_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!     call CCSM4_LGM_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call CNRM_CM5_PD_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call CNRM_CM5_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call COSMOS_ASO_PD_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call COSMOS_ASO_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call FGOALS_g2_PD_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call FGOALS_g2_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call GISS_E2_R_150_PD_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call GISS_E2_R_150_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call GISS_E2_R_151_PD_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call GISS_E2_R_151_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call IPSL_CM5A_LR_PD_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call IPSL_CM5A_LR_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call MIROC_ESM_PD_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call MIROC_ESM_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call MPI_ESM_P_p1_PD_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call MPI_ESM_P_p1_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call MPI_ESM_P_p2_PD_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call MPI_ESM_P_p2_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call MRI_CGCM3_PD_to_grid(outfldr,subfldr,grid,domain="present", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
!    call MRI_CGCM3_LGM_to_grid(outfldr,subfldr,grid,domain="lgm", &
!                               path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)

    

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gridder

