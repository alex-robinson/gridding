
module regions
    ! This module defines regions globally (lat-lon)
    ! and returns them on the coordinates of interest (eg, polar_stereographic)

    use coordinates 
    use polygons 

    implicit none 

    private
    public :: get_region_map_north 

contains 

    

    function get_region_map(grid) result(mask)
        ! For a given grid (input), output a mask of the regions
        ! that overlap with it

        implicit none

        type(grid_class), intent(IN) :: grid 
        integer :: mask(grid%G%nx,grid%G%ny)


        type(points_class), allocatable :: regs(:)
        
        ! First, generate all potential regions 
!         call create_regions(regs)

        ! First map latlon points to grid of interest, then check point_in_polygon
!         mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
!         where (mask_reg) outmask = 2 


        return 

    end function get_region_map 

    function get_region_map_north(grid) result(mask)
        ! For a given grid (input), output a mask of the regions
        ! that overlap with it

        implicit none

        type(grid_class), intent(IN) :: grid 
        integer :: mask(grid%G%nx,grid%G%ny)

        type(points_class), allocatable :: regs(:)
        type(points_class) :: pts 
        logical :: in_reg(grid%G%nx,grid%G%ny)
        integer :: q 

        ! Allocate the region_type to hold all regions of interest
        allocate(regs(5))

        ! === Define each region ===

!         call points_init(regs(1),grid0=grid,name="Ellesmere Island", &
!                          x = [-59.3d0,-60.0d0,-70.8d0,-78.0d0,-100.3d0,-98.8d0,-88.7d0 ], &
!                          y = [ 85.0d0, 82.4d0, 79.7d0, 76.0d0,  80.9d0, 81.2d0, 85.0d0 ], &
!                          latlon=.TRUE.)

!         call points_init(regs(2),grid0=grid,name="Iceland", &
!                          x = [-16.6d0,-22.3d0,-25.7d0,-28.8d0,-25.1d0 ], &
!                          y = [ 69.5d0, 69.2d0, 67.3d0, 65.7d0, 57.7d0 ], &
!                          latlon=.TRUE.)
        
!         call points_init(regs(3),grid0=grid,name="Svalbard", &
!                          x = [ 40.0d0, 40.0d0, 20.0d0,  0.0d0, -10.0d0 ], &
!                          y = [ 85.0d0, 80.0d0, 73.0d0, 75.0d0,  85.0d0 ], &
!                          latlon=.TRUE.)
        
        call points_init(regs(1),grid0=grid,name="grl",filename="regions/polygon_grl.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(2),grid0=grid,name="grl",filename="regions/polygon_grl-inner.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(3),grid0=grid,name="grl",filename="regions/polygon_ellesmere.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(4),grid0=grid,name="grl",filename="regions/polygon_iceland.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(5),grid0=grid,name="grl",filename="regions/polygon_svalbard.txt",latlon=.TRUE.,skip=1)

        write(*,*)
        write(*,*)
        write(*,*)
        write(*,*) "----------------------------------------"
        write(*,*) "grl: ",minval(regs(4)%lon), maxval(regs(4)%lon)
        write(*,*) "grl: ",minval(regs(4)%lat), maxval(regs(4)%lat)
        write(*,*) "----------------------------------------"
        write(*,*)
        write(*,*)
        write(*,*)

        mask = 0 
        do q = 1, size(regs)
            in_reg = point_in_polygon(real(grid%x),real(grid%y),real(regs(q)%x),real(regs(q)%y))
            where (in_reg) mask = q  
        end do 

        return 

    end function get_region_map_north

end module regions

