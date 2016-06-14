
module regions
    ! This module defines regions globally (lat-lon)
    ! and returns them on the coordinates of interest (eg, polar_stereographic)

    use coordinates 
    use polygons 

    implicit none 


    type region_type 
        character(len=128) :: name 
        double precision, allocatable :: lon(:), lat(:)

    end type



contains 

    

    function get_region_map(grid) result(mask)
        ! For a given grid (input), output a mask of the regions
        ! that overlap with it

        implicit 

        type(grid_class), intent(IN) :: grid 
        integer :: mask(grid%G%nx,grid%G%ny)


        ! First convert lon/lat points to grid of interest, then check point_in_polygon
!         mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
!         where (mask_reg) outmask = 2 


        return 

    end function get_region_map 

    subroutine create_regions(regs)

        implicit none 

        type(region_type), allocatable :: regs(:)

        ! Allocate the region_type to hold all regions of interest
        allocate(regs(10))

        ! === Define each region ===

        ! Ellesmere Island
        call region_init(regs(1),"Ellesmere Island",n=7)
        xp = [-59.3,-60.0,-70.8,-78.0,-100.3,-98.8,-88.7 ]
        yp = [ 85.0, 82.4, 79.7, 76.0,  80.9, 81.2, 85.0 ]
        
        ! Iceland 
        if (allocated(xp)) deallocate(xp)
        if (allocated(yp)) deallocate(yp)
        allocate(xp(5),yp(5))
        xp = [-16.6,-22.3,-25.7,-28.8,-25.1]
        yp = [ 69.5, 69.2, 67.3, 65.7, 57.7]
        mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
        where (mask_reg) outmask = 3 

        ! Svalbard
        if (allocated(xp)) deallocate(xp)
        if (allocated(yp)) deallocate(yp)
        allocate(xp(5),yp(5))
        xp = [ 40.0, 40.0, 20.0,  0.0, -10.0 ]
        yp = [ 85.0, 80.0, 73.0, 75.0,  85.0 ]
        mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
        where (mask_reg) outmask = 4 


        return 

    end subroutine create_regions


    subroutine region_init(reg,name,n)

        implicit none 

        type(region_type), intent(INOUT) :: reg 
        character(len=*),  intent(IN)    :: name 
        integer,           intent(IN)    :: n 

        reg%name = trim(name)
        if (allocated(reg%lon)) deallocate(reg%lon)
        allocate(reg%lon(n))
        if (allocated(reg%lat)) deallocate(reg%lat)
        allocate(reg%lat(n))

        return 

    end subroutine region_init 

end module regions

