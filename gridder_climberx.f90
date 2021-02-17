
program gridder

    use coord 
    use control 
    
    ! Datasets that can be gridded by this program (alphabetical)
    use climberx

    implicit none

    type(grid_class)   :: grid
    character(len=256) :: domain, grid_name, outfldr_ref, outfldr 
    character(len=256) :: path, subfldr  
    real(8) :: sigma1

    write(*,*) 

    ! =========================================================
    !
    ! USER DEFINITIONS
    !
    ! =========================================================
    
    domain      = "North"
    grid_name   = "NH-32KM"
    outfldr_ref = "output/"//trim(domain)//"/"//trim(grid_name)

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

    ! Define smoothing radius (km)
    sigma1 = 100.d0

    ! Define output folder
    write(subfldr,"(a12,i3,a2)") "climberx_sig", int(sigma1), "km"
    outfldr = trim(outfldr_ref)//"/"//trim(subfldr) 
    call system("mkdir -p "//trim(outfldr))
    write(*,*) "outfldr = ", trim(outfldr) 

    ! Preindustrial simulation 
    path = "/Users/robinson/models/climber-x/output/bench_v098/preind"
    call climberx_to_grid(outfldr,grid,expname="preind",infldr=path,sigma=sigma1)

    ! Default LGM simulation (shallower, weaker amoc)
    path = "/Users/robinson/models/climber-x/output/bench_v098/lgm"
    call climberx_to_grid(outfldr,grid,expname="lgm-default",infldr=path,sigma=sigma1)

    ! LGM simulation (strong amoc)
    path = "/Users/robinson/models/climber-x/output/bench_v098/lgm_amoc03_5070"
    call climberx_to_grid(outfldr,grid,expname="lgm-strong",infldr=path,sigma=sigma1)

    ! LGM simulation (amoc killed)
    path = "/Users/robinson/models/climber-x/output/bench_v098/lgm_amoc02fw"
    call climberx_to_grid(outfldr,grid,expname="lgm-off",infldr=path,sigma=sigma1)

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gridder

