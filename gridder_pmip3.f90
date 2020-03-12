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

    character(len=256), allocatable :: models(:) 
    character(len=256) :: model 
    integer :: n, nmodels 

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
    
    nmodels = 11 

    allocate(models(nmodels))
    models(1)  = "CCSM4" 
    models(2)  = "CNRM_CM5"
    models(3)  = "COSMOS_ASO"
    models(4)  = "FGOALS_g2"
    models(5)  = "GISS_E2_R_150"
    models(6)  = "GISS_E2_R_151"
    models(7)  = "IPSL_CM5A_LR"
    models(8)  = "MIROC_ESM"
    models(9)  = "MPI_ESM_P_p1"
    models(10) = "MPI_ESM_P_p2"
    models(11) = "MRI_CGCM3"
    
    do n = 1, nmodels 
        
        model = models(n) 

        call PMIP3_to_grid(outfldr,grid,domain,model=trim(model),experiment="piControl", &
                                  path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
        call PMIP3_to_grid(outfldr,grid,domain,model=trim(model),experiment="lgm", &
                                  path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)

    end do 

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gridder

