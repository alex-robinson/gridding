
module regions
    ! This module defines regions globally (lat-lon)
    ! and returns them on the coordinates of interest (eg, polar_stereographic)

    use coordinates 
    use polygons 

    implicit none 

    private
    public :: get_region_map_greenland 

contains 

    function get_region_map_greenland(grid) result(mask)
        ! For a given grid (input), output a mask of the regions
        ! that overlap with it

        implicit none

        type(grid_class), intent(IN) :: grid 
        real(4) :: mask(grid%G%nx,grid%G%ny)

        type(points_class), allocatable :: regs(:)
        type(points_class) :: pts 
        logical :: in_reg(grid%G%nx,grid%G%ny)
        integer :: q 

        real(4), parameter :: mask_index = 1.0 

        ! Allocate the region_type to hold all regions of interest
        allocate(regs(5))

        ! === Define each region ===
        call points_init(regs(1),grid0=grid,name="grl",filename="regions/polygon_grl.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(2),grid0=grid,name="grl",filename="regions/polygon_grl-inner.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(3),grid0=grid,name="grl",filename="regions/polygon_ellesmere.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(4),grid0=grid,name="grl",filename="regions/polygon_svalbard.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(5),grid0=grid,name="grl",filename="regions/polygon_iceland.txt",latlon=.TRUE.,skip=1)
        
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

        mask = mask_index + 0.0   ! ocean 
        do q = 1, size(regs)
            in_reg = point_in_polygon(real(grid%x),real(grid%y),real(regs(q)%x),real(regs(q)%y))
            where (in_reg) mask = mask_index + real(q)*0.1  
        end do 

        return 

    end function get_region_map_greenland

end module regions

