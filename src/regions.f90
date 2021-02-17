
module regions
    ! This module defines regions globally (lat-lon)
    ! and returns them on the coordinates of interest (eg, polar_stereographic)

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    real(4), parameter :: mask_index_north = 1.0  ! North
    real(4), parameter :: mask_index_south = 2.0  ! Antarctica

    private
    public :: get_region_map_north
    public :: get_region_map_south
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
        logical :: in_continent(grid%G%nx,grid%G%ny)
        real(4), allocatable :: reg_vals(:) 
        integer :: q, n_regions, k
        
        ! Allocate the region_type to hold all regions of interest
        n_regions = 10 
        allocate(regs(n_regions))
        allocate(reg_vals(n_regions))

        ! === Define each region ===

        ! == Continental regions ==
        call points_init(regs(1), grid0=grid,name="reg1",filename="regions/polygons/polygon_laurentide.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(2), grid0=grid,name="reg2",filename="regions/polygons/polygon_eis.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(3), grid0=grid,name="reg3",filename="regions/polygons/polygon_grl_and_ellesmere.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(4), grid0=grid,name="reg4",filename="regions/polygons/polygon_asia.txt",latlon=.TRUE.,skip=1)
        
        ! == Sub-regions ==
        call points_init(regs(5), grid0=grid,name="reg5", filename="regions/polygons/polygon_ellesmere.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(6), grid0=grid,name="reg6", filename="regions/polygons/polygon_barents-kara.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(7), grid0=grid,name="reg7", filename="regions/polygons/polygon_britain.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(8), grid0=grid,name="reg8", filename="regions/polygons/polygon_svalbard.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(9), grid0=grid,name="reg9", filename="regions/polygons/polygon_iceland.txt",latlon=.TRUE.,skip=1)
        call points_init(regs(10),grid0=grid,name="reg10",filename="regions/polygons/polygon_hudson.txt",latlon=.TRUE.,skip=1)
        
        ! == Continental regions ==
        reg_vals(1)  = 1.0
        reg_vals(2)  = 2.0
        reg_vals(3)  = 3.0
        reg_vals(4)  = 4.0

        ! == Sub-regions ==
        reg_vals(5)  = 1.1
        reg_vals(6)  = 2.1
        reg_vals(7)  = 2.2
        reg_vals(8)  = 2.3
        reg_vals(9)  = 3.1
        reg_vals(10) = 1.2 
        
        mask = mask_index_north + 0.0   ! ocean 

        do q = 1, n_regions

            ! Determine the continent we are working with 
            if (reg_vals(q) .eq. floor(reg_vals(q))) then 
                ! This region is a continent itself, every point is allowed 
                ! and will be determined soley by in_reg later. 
                in_continent = .TRUE. 
            else 
                ! This region is a subregion of a continent. Determine
                ! the extent of its continent. 

                k = minloc(abs(reg_vals - floor(reg_vals(q))),1)
                in_continent = point_in_polygon(real(grid%x),real(grid%y),real(regs(k)%x),real(regs(k)%y))
                
            end if 

            ! Check which grid points are contained in the current region defined 
            ! by a polygon - return a mask in_reg(nx,ny) of true/false
            in_reg = point_in_polygon(real(grid%x),real(grid%y),real(regs(q)%x),real(regs(q)%y))
            
            ! Assume that subregions have the same integer value as the continental 
            ! region that contains them. Don't allow subregions to have points outside 
            ! of their continental region 
            where (in_continent .and. in_reg) mask = mask_index_north + reg_vals(q)*0.1  
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
        real(4), allocatable :: reg_vals(:) 
        integer :: q, n_regions 
        
        ! Allocate the region_type to hold all regions of interest
        n_regions = 2 
        allocate(regs(n_regions))
        allocate(reg_vals(n_regions))

        ! === Define each region ===

        ! == Continental regions ==
        call points_init(regs(1),grid0=grid,name="reg1",filename="regions/polygons/polygon_antarctica.txt",latlon=.TRUE.,skip=1)

        ! == Sub-regions ==
        call points_init(regs(2),grid0=grid,name="reg5",filename="regions/polygons/polygon_antarctica_inner.txt",latlon=.TRUE.,skip=1)
        
        ! ajr, to do!!
        !call points_init(regs(2),grid0=grid,name="reg5",filename="regions/polygons/polygon_wais.txt",latlon=.TRUE.,skip=1)
        !call points_init(regs(3),grid0=grid,name="reg6",filename="regions/polygons/polygon_eais.txt",latlon=.TRUE.,skip=1)

        ! == Continental regions ==
        reg_vals(1) = 1.0

        ! == Sub-regions ==
        reg_vals(2) = 1.1

        ! ajr, to do!!
        !reg_vals(3) = 1.2
        !reg_vals(4) = 1.3
        
        mask = mask_index_south + 0.0   ! ocean 

        do q = 1, n_regions
            in_reg = point_in_polygon(real(grid%x),real(grid%y),real(regs(q)%x),real(regs(q)%y))
            where (in_reg) mask = mask_index_south + reg_vals(q)*0.1  
        end do 

        return 

    end function get_region_map_south

    subroutine write_regions(outfldr,grid,domain)

        implicit none 

        character(len=*), intent(IN) :: outfldr
        type(grid_class), intent(IN) :: grid 
        character(len=*), intent(IN) :: domain 

        character(len=512)   :: filename 
        real(4), allocatable :: mask(:,:) 
        real(4), allocatable :: mask2(:,:) 

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_REGIONS.nc"

        ! Allocate the mask 
        call grid_allocate(grid,mask)    
        
        select case(trim(domain)) 

            case("North","Laurentide","Eurasia","Greenland") 
                mask = get_region_map_north(grid)

            case("Antarctica") 
                mask = get_region_map_south(grid)
                
!             case("Greenland")
!                 mask = get_region_map_greenland(grid)
            
            case("Global") 
                mask  = get_region_map_north(grid)
                mask2 = get_region_map_south(grid)
                
                where(mask .eq. 0.0) mask = mask2 
                
            case DEFAULT 
                write(*,*) "regions:: error: domain not recognized: "//trim(domain)
                write(*,*) "Setting all values of regions to zero."
                mask = 0.0  

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

