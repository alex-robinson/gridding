
module regions
    ! This module defines regions globally (lat-lon)
    ! and returns them on the coordinates of interest (eg, polar_stereographic)

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    real(4), parameter :: mask_index_north = 1.0  ! North
    real(4), parameter :: mask_index_south = 2.0  ! Antarctica
    real(4), parameter :: mask_index_grl   = 3.0  ! Greenland

    private
    public :: get_region_map_north
    public :: get_region_map_greenland 
    public :: write_regions
contains 

    function get_region_map_north(grid) result(mask)
        ! For a given grid (input), output a mask of the regions
        ! that overlap with it

        implicit none

        type(grid_class), intent(IN) :: grid 
        real(4) :: mask(grid%G%nx,grid%G%ny)

        type(points_class), allocatable :: regs(:)
        type(points_class) :: pts 
        logical :: in_reg(grid%G%nx,grid%G%ny)
        integer :: q 

        ! Allocate the region_type to hold all regions of interest
        allocate(regs(5))

        ! === Define each region ===
        call points_init(regs(1),grid0=grid,name="grl",filename="regions/polygon_grl.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(2),grid0=grid,name="grl",filename="regions/polygon_grl-inner.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(3),grid0=grid,name="grl",filename="regions/polygon_ellesmere.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(4),grid0=grid,name="grl",filename="regions/polygon_svalbard.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(5),grid0=grid,name="grl",filename="regions/polygon_iceland.txt",latlon=.TRUE.,skip=1)
        
        mask = mask_index_north + 0.0   ! ocean 

        do q = 1, size(regs)
            in_reg = point_in_polygon(real(grid%x),real(grid%y),real(regs(q)%x),real(regs(q)%y))
            where (in_reg) mask = mask_index_grl + real(q)*0.1  
        end do 

        return 

    end function get_region_map_north

    function get_region_map_south(grid) result(mask)
        ! For a given grid (input), output a mask of the regions
        ! that overlap with it

        implicit none

        type(grid_class), intent(IN) :: grid 
        real(4) :: mask(grid%G%nx,grid%G%ny)

        type(points_class), allocatable :: regs(:)
        type(points_class) :: pts 
        logical :: in_reg(grid%G%nx,grid%G%ny)
        integer :: q 

        mask = mask_index_south + 0.0   ! ocean 

        return 

    end function get_region_map_south

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

        ! Allocate the region_type to hold all regions of interest
        allocate(regs(5))

        ! === Define each region ===
        call points_init(regs(1),grid0=grid,name="grl",filename="regions/polygon_grl.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(2),grid0=grid,name="grl",filename="regions/polygon_grl-inner.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(3),grid0=grid,name="grl",filename="regions/polygon_ellesmere.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(4),grid0=grid,name="grl",filename="regions/polygon_svalbard.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(5),grid0=grid,name="grl",filename="regions/polygon_iceland.txt",latlon=.TRUE.,skip=1)
        
        mask = mask_index_grl + 0.0   ! ocean 
        do q = 1, size(regs)
            in_reg = point_in_polygon(real(grid%x),real(grid%y),real(regs(q)%x),real(regs(q)%y))
            where (in_reg) mask = mask_index_grl + real(q)*0.1  
        end do 

        return 

    end function get_region_map_greenland

    subroutine write_regions(outfldr,grid,domain)

        implicit none 

        character(len=*), intent(IN) :: outfldr
        type(grid_class), intent(IN) :: grid 
        character(len=*), intent(IN) :: domain 

        character(len=512)   :: filename 
        real(4), allocatable :: mask(:,:) 

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_REGIONS.nc"

        ! Allocate the mask 
        call grid_allocate(grid,mask)    
        
        select case(trim(domain)) 

            case("North") 
                mask = get_region_map_north(grid)

            case("Antarctica") 
                mask = get_region_map_south(grid)
                
            case("Greenland")
                mask = get_region_map_greenland(grid)

            case DEFAULT 
                write(*,*) "regions:: error: domain not recognized: "//trim(domain)
                stop 

        end select

        ! Write the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        call nc_write(filename,"mask",mask,dim1="xc",dim2="yc", &  
                      long_name="Regional mask")

        return 

    end subroutine write_regions 

end module regions

