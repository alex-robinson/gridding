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
    grid_name = "GRL-16KM"
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

    nmodels = 7 

    allocate(models(nmodels))
    models(1)  = "CCSM4" 
    models(2)  = "CNRM-CM5"
    models(3)  = "FGOALS-g2"
    models(4)  = "IPSL-CM5A-LR"
    models(5)  = "MIROC-ESM"
    models(6)  = "MPI-ESM-P"
    models(7)  = "MRI-CGCM3"
    
    ! Excluded models
    !models(4)  = "GISS-E2-R"        ! No sst data available
    
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

