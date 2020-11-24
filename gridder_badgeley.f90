program gridder
  
    use coord
    use control

    ! Datasets that can be gridded by this program (alphabetical)
    use badgeley2020

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

    path = "data/badgeley2020/"
    sigma = 50.d0 

    write(sigma_str,*) int(sigma)
    sigma_str = trim(adjustl(sigma_str))

    write(subfldr,"(a15,a,a2)") "B20_sig", trim(sigma_str), "km"
    subfldr = trim(adjustl(subfldr))
    outfldr = trim(outfldr)//"/"//trim(subfldr)

    call system("mkdir -p "//trim(outfldr))

    call badgeley_to_grid(outfldr,subfldr,grid,domain="holocene-retreat", &
                              path_in=path,sigma=sigma,max_neighbors=10,lat_lim=5d0)

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gridder
