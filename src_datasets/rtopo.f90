module rtopo 

    use coord
    use ncio
    use gridding_datasets
    use generic 
    
    implicit none 

    private 
    public :: rtopo_latlon_to_grid 
!     public :: rtopo_to_grid 

contains 

    subroutine rtopo_latlon_to_grid(outfldr,grid,domain)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       RTOPO DATA
        !       * These data files are enormous (3.5 Gb). Instead of
        !       * using mapping routines, that may crash due to memory
        !       * requirements, the first step will be a quick and dirty
        !       * nearest neighbor solution to a high resolution, local 
        !       * grid. 
        !
        ! =========================================================
        ! PATH: /p/projects/megarun/greenrise/datasets/rtopo-2/
        ! londim = 43201
        ! latdim = 21601
        ! lon, lat 
        ! 
        ! RTopo-2.0.1_30sec_bedrock_topography.nc : bedrock_topography
        ! long_name = "ocean bathymetry; surface topography of continents; bedrock topography under grounded or floating ice"
        ! RTopo-2.0.1_30sec_aux.nc : amask
        ! long_name = "ice ocean rock mask"
        ! RTopo-2.0.1_30sec_ice_base_topography.nc : ice_base_topography
        ! long_name = "ice base topography for the Antarctic and Greenland ice sheets / ice shelves (ice draft for ice shelves and floating glaciers; zero in absence of ice)"
        ! RTopo-2.0.1_30sec_surface_elevation.nc : surface_elevation
        ! long_name = "upper ice surface height for the Antarctic and Greenland ice sheets / ice shelves (bedrock topography for ice-free continent; zero for ocean)"


        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid
        double precision :: lat_lim
        character(len=1024) :: path
        character(len=512)  :: filename_in, filename  
        character(len=1024) :: desc, ref 
        character(len=1024) :: var_name, var_name_out, long_name, units  
        
        type inp_type 
            real(4), allocatable :: lon(:), lat(:)
            real(4), allocatable :: var(:,:)
        end type 

        type(inp_type)     :: inp
        integer :: nx0, ny0 

        type neighbs_type 
            integer, allocatable :: i(:,:), j(:,:) 
        end type 
        type(neighbs_type) :: nbs 

        real(4), allocatable :: var(:,:) 

        ! ### Define input information #####
        path = "/p/projects/megarun/greenrise/datasets/rtopo-2"
        filename_in = trim(path)//"/"//"RTopo-2.0.1_30sec_bedrock_topography.nc"

        ! Allocate input field dimensions
        nx0 = nc_size(filename_in,"londim")
        ny0 = nc_size(filename_in,"latdim")
        allocate(inp%lon(nx0),inp%lat(ny0))
        allocate(inp%var(nx0,ny0))

        write(*,*) "nx0, ny0: ", nx0, ny0 

        ! Load lat/lon 
        call nc_read(filename_in,"lon",inp%lon)
        call nc_read(filename_in,"lat",inp%lat)
        
        write(*,*) "range(lon): ", minval(inp%lon), maxval(inp%lon)
        write(*,*) "range(lat): ", minval(inp%lat), maxval(inp%lat)
        
        desc    = "RTOPO-2.0.1 present-day Earth topography data"
        ref     = "Timmermann et al.: A consistent data set of Antarctic &
                  &ice sheet topography, cavity geometry, and global bathymetry, &
                  &Earth Syst. Sci. Data, 2, 261-273, doi:10.5194/essd-2-261-2010, 2010.\n &
                  &https://doi.pangaea.de/10.1594/PANGAEA.741917"

        ! ### Define output information #####
        
        ! Allocate output variable
        call grid_allocate(grid,var)
        call grid_allocate(grid,nbs%i)
        call grid_allocate(grid,nbs%j)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_TOPO-RTOPO-2.0.1.nc"

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ### Process each variable, interpolate to nearest neighbor on grid #####

        ! 1. Bedrock topography
        var_name     = "bedrock_topography"
        var_name_out = "z_bed"
        long_name    = "ocean bathymetry; surface topography of continents; &
                       &bedrock topography under grounded or floating ice"
        units        = "m" 

        path = "/p/projects/megarun/greenrise/datasets/rtopo-2"
        filename_in = trim(path)//"/"//"RTopo-2.0.1_30sec_bedrock_topography.nc"
        call nc_read(filename_in,var_name,inp%var,missing_value=real(mv))

        write(*,*) "input range(var):  ", minval(inp%var), maxval(inp%var)

        ! Interpolate to output grid 
        call nearest_to_grid(zout=var,grid=grid,x=inp%lon,y=inp%lat,z=inp%var,latlon=.TRUE., &
                             max_dist=10e3,lat_lim=0.1)

        write(*,*) "output range(var): ", minval(var), maxval(var)
        
        ! Write output variable to output file
        call nc_write(filename,var_name_out,real(var),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,var_name_out,"units",units)
        call nc_write_attr(filename,var_name_out,"long_name",long_name)
        call nc_write_attr(filename,var_name_out,"coordinates","lat2D lon2D")
            

        return 

    end subroutine rtopo_latlon_to_grid

    
!     subroutine rtopo_to_grid(outfldr,grid,domain)
!         ! Convert the variables to the desired grid format and write to file
!         ! =========================================================
!         !
!         !       RTOPO DATA
!         !       * Given rtopo data on a regional grid, perform 
!         !       * normal interpolation using mapping methods
!         !
!         ! =========================================================
!         ! PATH: /p/projects/megarun/greenrise/datasets/rtopo-2/
!         ! londim = 43201
!         ! latdim = 21601
!         ! lon, lat 
!         ! 
!         ! RTopo-2.0.1_30sec_bedrock_topography.nc : bedrock_topography
!         ! long_name = "ocean bathymetry; surface topography of continents; bedrock topography under grounded or floating ice"
!         ! RTopo-2.0.1_30sec_aux.nc : amask
!         ! long_name = "ice ocean rock mask"
!         ! RTopo-2.0.1_30sec_ice_base_topography.nc : ice_base_topography
!         ! long_name = "ice base topography for the Antarctic and Greenland ice sheets / ice shelves (ice draft for ice shelves and floating glaciers; zero in absence of ice)"
!         ! RTopo-2.0.1_30sec_surface_elevation.nc : surface_elevation
!         ! long_name = "upper ice surface height for the Antarctic and Greenland ice sheets / ice shelves (bedrock topography for ice-free continent; zero for ocean)"


!         implicit none 

!         character(len=*) :: domain, outfldr 
!         type(grid_class) :: grid 
!         integer :: max_neighbors 
!         double precision :: lat_lim, grad_lim 
!         character(len=512) :: filename 
!         character(len=1024) :: desc, ref 


!         ! Initialize output file 
!         file_out = trim(outfldr)//trim(grid%name)//trim(file_out_suffix)
!         call grid_write(mgrid%grid(q),file_out,xnm="xc",ynm="yc",create=.TRUE.)
            

!         varnames = ["zs  ","zb  ","H   ","mask"]

!         do k = 1, size(varnames)

!             varname = varnames(k)
!             write(*,*) k, trim(varname)

!             ! Load test data
!             call nc_read(file_input,varname,var)
!             target_val = calc_grid_total(grid0%G%x,grid0%G%y,var,xlim=xlim,ylim=ylim)

!             if (trim(varname) .ne. "mask") then
!                 call map_field_conservative_map1(map%map,varname,var,var1)

!                 current_val = calc_grid_total(grid%G%x,grid%G%y,var1,xlim=xlim,ylim=ylim)
!                 err_percent = 100.d0 * (current_val-target_val) / target_val
!                 write(*,"(a,3g12.4)") "mass comparison (hi, con, % diff): ", &
!                         target_val, current_val, err_percent                  

!             else 
!                 var1 = interp_nearest(x=grid0%G%x,y=grid0%G%y,z=var, &
!                                       xout=grid%G%x,yout=grid%G%y)

!             end if 

!             ! Write to file 
!             call nc_write(file_out,varname,var1,dim1="xc",dim2="yc")
 
!         end do 

!         return 

!     end subroutine rtopo_to_grid


end module rtopo 
