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
    
    domain    = "Laurentide"
    grid_name = "LIS-20KM"
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

    path = "data/PMIP3/all_PMIP_cvdp_data/data/"
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
    
    nmodels = 8 

    allocate(models(nmodels))
    models(1)  = "CCSM4" 
    models(2)  = "CNRM-CM5"
    models(3)  = "FGOALS-g2"
    models(4)  = "GISS-E2-R"
    models(5)  = "IPSL-CM5A-LR"
    models(6)  = "MIROC-ESM"
    models(7)  = "MPI-ESM-P"
    models(8)  = "MRI-CGCM3"
    

    do n = 1, nmodels 
        
        model = models(n) 

        call pmip3_to_grid(outfldr,grid,domain,model=trim(model),experiment="piControl", &
                                  path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
        call pmip3_to_grid(outfldr,grid,domain,model=trim(model),experiment="midHolocene", &
                                  path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)
        call pmip3_to_grid(outfldr,grid,domain,model=trim(model),experiment="lgm", &
                                  path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5.d0)

    end do 

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gridder


! From Javi's data processing:
!     nmodels = 11 

!     allocate(models(nmodels))
!     models(1)  = "CCSM4" 
!     models(2)  = "CNRM-CM5"
!     models(3)  = "COSMOS_ASO"
!     models(4)  = "FGOALS_g2"
!     models(5)  = "GISS_E2_R_150"
!     models(6)  = "GISS_E2_R_151"
!     models(7)  = "IPSL_CM5A_LR"
!     models(8)  = "MIROC_ESM"
!     models(9)  = "MPI_ESM_P"
!     models(11) = "MRI_CGCM3"

