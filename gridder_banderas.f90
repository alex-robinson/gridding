
program gridder

    use coord 
    use control 
    
    ! Datasets that can be gridded by this program (alphabetical)
    use banderas2018

    implicit none

    type(grid_class)   :: grid
    character(len=256) :: domain, grid_name, outfldr_ref, outfldr 
    character(len=256) :: path, subfldr  

    write(*,*) 

    ! =========================================================
    !
    ! USER DEFINITIONS
    !
    ! =========================================================
    
    domain      = "Global"
    grid_name   = "latlon-05deg"
    outfldr     = "output/"//trim(domain)//"/"//trim(grid_name)

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

    ! Banderas et al (2018)
    call banderas2018_to_grid(outfldr,grid,domain)

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gridder

