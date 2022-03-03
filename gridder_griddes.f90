program gridder_griddes
    ! Use this program to generate cdo-compliant grid description files. 

    use coord 
    use control 
    
    use grid_to_cdo 

    implicit none

    character(len=256) :: domain, grid_name  
    type(grid_class)   :: grid
    
    write(*,*) 

    ! =========================================================
    !
    ! USER DEFINITIONS
    !
    ! =========================================================
        
    ! Default case:
    domain    = "Greenland"
    grid_name = "GRL-16KM"
    
    ! Replace with command line arguments:
    call load_command_line_args(domain,grid_name)

    ! =========================================================
    !
    ! OUTPUT GRID DEFINITION
    !
    ! =========================================================

    call domain_definition(grid,grid_name)  

    ! Write a regional mask 
    !call write_regions(outfldr,grid,domain) 
    
    call grid_cdo_write_desc_short(grid,fldr="maps")

    write(*,*) "Grid description file written." 

contains 

    subroutine load_command_line_args(domain,grid_name)
        ! Load the parameter filename from the command line 
        ! call eg: ./yelmo_test.x yelmo_Greenland.nml 

        implicit none 

        character(len=*), intent(OUT) :: domain
        character(len=*), intent(OUT) :: grid_name

        ! Local variables 
        integer :: narg 

        narg = command_argument_count()

        if (narg .ne. 2) then 
            write(*,*) "load_command_line_args:: Error: The following &
            &argument must be provided: domain, grid_name"
            stop 
        end if 

        call get_command_argument(1,domain)
        call get_command_argument(2,grid_name)

        return 

    end subroutine load_command_line_args

end program gridder_griddes