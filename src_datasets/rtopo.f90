module rtopo 

    use coord
    use ncio
    use gridding_datasets
    use generic 
    use control 

    implicit none 

    private 
    public :: rtopo_latlon_to_grid 
    public :: rtopo_to_grid 

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
        ! RTopo-2.0.1_30sec_surface_elevation.nc : surface_elevation
        ! long_name = "upper ice surface height for the Antarctic and Greenland ice sheets / ice shelves (bedrock topography for ice-free continent; zero for ocean)"
        ! RTopo-2.0.1_30sec_ice_base_topography.nc : ice_base_topography
        ! long_name = "ice base topography for the Antarctic and Greenland ice sheets / ice shelves (ice draft for ice shelves and floating glaciers; zero in absence of ice)"
        ! RTopo-2.0.1_30sec_aux.nc : amask
        ! long_name = "ice ocean rock mask"
        

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
            integer :: i0, i1, j0, j1, ni, nj    ! Range of subset of data of interest
        end type 

        type(inp_type)     :: inp

        integer :: nx0, ny0 

        type neighbs_type 
            integer, allocatable :: ii(:,:), jj(:,:) 
        end type 
        type(neighbs_type) :: nbs 

        real(4), allocatable :: var(:,:)
        integer :: k 

        ! ### Define input information #####
        path = "/p/projects/megarun/greenrise/datasets/rtopo-2"
        filename_in = trim(path)//"/"//"RTopo-2.0.1_30sec_bedrock_topography.nc"

        ! Allocate input field dimensions
        nx0 = nc_size(filename_in,"londim")
        ny0 = nc_size(filename_in,"latdim")
        allocate(inp%lon(nx0),inp%lat(ny0))
        write(*,*) "nx0, ny0: ", nx0, ny0 

        ! Load lat/lon 
        call nc_read(filename_in,"lon",inp%lon)
        call nc_read(filename_in,"lat",inp%lat)
        
        write(*,*) "range(lon): ", minval(inp%lon), maxval(inp%lon)
        write(*,*) "range(lat): ", minval(inp%lat), maxval(inp%lat)
        
        ! Determine subset of interest

        ! Take all longitudes for now 
        inp%i0 = 1
        inp%i1 = size(inp%lon)

        ! # Longitude range
        inp%i0 = size(inp%lon)
        do k = 1, size(inp%lon)
            if (inp%lon(k) .ge. minval(grid%lon)-1.d0) exit   
        end do
        inp%i0 = k 
        inp%i1 = inp%i0 
        do k = inp%i0, size(inp%lon)
            if (inp%lon(k) .ge. maxval(grid%lon)+1.d0) exit   
        end do
        inp%i1 = k

        ! # Latitude range
        inp%j0 = size(inp%lat)
        do k = 1, size(inp%lat)
            if (inp%lat(k) .ge. minval(grid%lat)-1.d0) exit   
        end do
        inp%j0 = k 
        inp%j1 = inp%j0 
        do k = inp%j0, size(inp%lat)
            if (inp%lat(k) .ge. maxval(grid%lat)+1.d0) exit   
        end do
        inp%j1 = k

        inp%ni = inp%i1-inp%i0+1
        inp%nj = inp%j1-inp%j0+1 

        ! Allocate input array to the size of the reduced matrix 
        allocate(inp%var(inp%ni,inp%nj))
        
        write(*,*) "range(lon): ", minval(inp%lon(inp%i0:inp%i1)), maxval(inp%lon(inp%i0:inp%i1))
        write(*,*) "range(lat): ", minval(inp%lat(inp%j0:inp%j1)), maxval(inp%lat(inp%j0:inp%j1))
        write(*,*) "ni, nj:     ", inp%ni, inp%nj 

        desc    = "RTOPO-2.0.1 present-day Earth topography data"
        ref     = "Timmermann et al.: A consistent data set of Antarctic &
                  &ice sheet topography, cavity geometry, and global bathymetry, &
                  &Earth Syst. Sci. Data, 2, 261-273, doi:10.5194/essd-2-261-2010, 2010. &
                  &Data download: https://doi.pangaea.de/10.1594/PANGAEA.741917"

        ! ### Define output information #####
        
        ! Allocate output variable
        call grid_allocate(grid,var)

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

        ! ### Find indices of nearest neighbors #####

        ! Allocate neighbor maps
        call grid_allocate(grid,nbs%ii)
        call grid_allocate(grid,nbs%jj)

        ! Search for nearest neighbors
        call find_nearest_grid(nbs%ii,nbs%jj,x=inp%lon(inp%i0:inp%i1),y=inp%lat(inp%j0:inp%j1), &
                               xout=real(grid%lon),yout=real(grid%lat),latlon=.TRUE., &
                               max_dist=2e3,lat_lim=0.05)
            

        ! ### Process each variable, interpolate to nearest neighbor on grid #####

        ! Source path for all input files 
        path = "/p/projects/megarun/greenrise/datasets/rtopo-2"

        ! 1. Bedrock topography ------------------------------------------------
        var_name     = "bedrock_topography"
        var_name_out = "z_bed"
        long_name    = "ocean bathymetry; surface topography of continents; &
                       &bedrock topography under grounded or floating ice"
        units        = "m" 

        filename_in = trim(path)//"/"//"RTopo-2.0.1_30sec_bedrock_topography.nc"
        call nc_read(filename_in,var_name,inp%var,missing_value=real(mv), &
                     start=[inp%i0,inp%j0],count=[inp%ni,inp%nj])

        write(*,*) "input range(var):  ", minval(inp%var), maxval(inp%var)

        ! Interpolate to output grid 
        call nearest_to_grid(zout=var,z=inp%var,ii=nbs%ii,jj=nbs%jj)

        write(*,*) "output range(var): ", minval(var,mask=var.ne.mv), maxval(var,mask=var.ne.mv)
        
        ! Write output variable to output file
        call nc_write(filename,var_name_out,real(var),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,var_name_out,"units",units)
        call nc_write_attr(filename,var_name_out,"long_name",long_name)
        call nc_write_attr(filename,var_name_out,"coordinates","lat2D lon2D")
            
        ! 2. Surface elevation ------------------------------------------------
        var_name     = "surface_elevation"
        var_name_out = "z_srf"
        long_name    = "upper ice surface height for the Antarctic and Greenland ice sheets &
                       &/ ice shelves (bedrock topography for ice-free continent; zero for ocean)"
        units        = "m" 

        filename_in = trim(path)//"/"//"RTopo-2.0.1_30sec_surface_elevation.nc"
        call nc_read(filename_in,var_name,inp%var,missing_value=real(mv), &
                     start=[inp%i0,inp%j0],count=[inp%ni,inp%nj])

        write(*,*) "input range(var):  ", minval(inp%var), maxval(inp%var)

        ! Interpolate to output grid 
        call nearest_to_grid(zout=var,z=inp%var,ii=nbs%ii,jj=nbs%jj)

        write(*,*) "output range(var): ", minval(var,mask=var.ne.mv), maxval(var,mask=var.ne.mv)
        
        ! Write output variable to output file
        call nc_write(filename,var_name_out,real(var),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,var_name_out,"units",units)
        call nc_write_attr(filename,var_name_out,"long_name",long_name)
        call nc_write_attr(filename,var_name_out,"coordinates","lat2D lon2D")
        
        ! 3. Ice base topography ------------------------------------------------
        var_name     = "ice_base_topography"
        var_name_out = "z_ice_base"
        long_name    = "ice base topography for the Antarctic and Greenland ice sheets &
                       &/ ice shelves (ice draft for ice shelves and floating glaciers; &
                       &zero in absence of ice)"
        units        = "m" 

        filename_in = trim(path)//"/"//"RTopo-2.0.1_30sec_ice_base_topography.nc"
        call nc_read(filename_in,var_name,inp%var,missing_value=real(mv), &
                     start=[inp%i0,inp%j0],count=[inp%ni,inp%nj])

        write(*,*) "input range(var):  ", minval(inp%var), maxval(inp%var)

        ! Interpolate to output grid 
        call nearest_to_grid(zout=var,z=inp%var,ii=nbs%ii,jj=nbs%jj)

        write(*,*) "output range(var): ", minval(var,mask=var.ne.mv), maxval(var,mask=var.ne.mv)
        
        ! Write output variable to output file
        call nc_write(filename,var_name_out,real(var),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,var_name_out,"units",units)
        call nc_write_attr(filename,var_name_out,"long_name",long_name)
        call nc_write_attr(filename,var_name_out,"coordinates","lat2D lon2D")
        
        ! 4. amask ------------------------------------------------
        var_name     = "amask"
        var_name_out = "mask"
        long_name    = "ice ocean rock mask"
        units        = "--" 

        filename_in = trim(path)//"/"//"RTopo-2.0.1_30sec_aux.nc"
        call nc_read(filename_in,var_name,inp%var,missing_value=real(mv), &
                     start=[inp%i0,inp%j0],count=[inp%ni,inp%nj])

        write(*,*) "input range(var):  ", minval(inp%var), maxval(inp%var)

        ! Interpolate to output grid 
        call nearest_to_grid(zout=var,z=inp%var,ii=nbs%ii,jj=nbs%jj)

        write(*,*) "output range(var): ", minval(var,mask=var.ne.mv), maxval(var,mask=var.ne.mv)
        
        ! Write output variable to output file
        call nc_write(filename,var_name_out,real(var),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,var_name_out,"units",units)
        call nc_write_attr(filename,var_name_out,"long_name",long_name)
        call nc_write_attr(filename,var_name_out,"coordinates","lat2D lon2D")
        
        return 

    end subroutine rtopo_latlon_to_grid

    
    subroutine rtopo_to_grid(outfldr,grid,domain)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       RTOPO DATA
        !       * Given rtopo data on a regional grid, perform 
        !       * normal interpolation using mapping methods
        !
        ! =========================================================

        implicit none 

        character(len=*), intent(IN) :: outfldr, domain
        type(grid_class), intent(IN) :: grid  

        ! Local variables
        type(map_class) :: map  
        character(len=512)  :: path, filename_in, filename  
        character(len=1024) :: desc, ref 
        character(len=56), allocatable :: varnames(:) 
        character(len=56)  :: varname
        character(len=256) :: units, long_name  
        integer :: k  
        real(8) :: current_val, target_val, err_percent

        type(grid_class) :: grid0 
        real(8), allocatable :: var0(:,:), var(:,:)
        real(8) :: xlim(2), ylim(2) 

        ! ### Input information #####

        ! Define input grid 
        select case(trim(domain))
            case("Greenland")
                call domain_definition(grid0,"GRL-2KM")
            case("Antarctica")
                call domain_definition(grid0,"ANT-1KM")
            case("North")
                call domain_definition(grid0,"NH-1KM")
            case DEFAULT 
                write(*,*) "rtopo_to_grid:: error: domain not recognized: "//trim(domain)
                stop

        end select 

        xlim = [minval(grid0%G%x), maxval(grid0%G%x)]
        ylim = [minval(grid0%G%y), maxval(grid0%G%y)]

        call grid_allocate(grid0,var0) 

        path = "output/"//trim(domain)
        filename_in = trim(path)//"/"//trim(grid0%name)//"/"// &
                        trim(grid0%name)//"_TOPO-RTOPO-2.0.1.nc"

        ! ### Output information ##### 

        call grid_allocate(grid,var) 
        
        ! Initialize the map 
        call map_init(map,grid0,grid,max_neighbors=20,lat_lim=1.d0,fldr="maps",load=.TRUE.)
        
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

        allocate(varnames(4))
        varnames = ["z_bed     ","z_srf     ","z_ice_base","mask      "]

        do k = 1, size(varnames)

            varname = varnames(k)
            write(*,*) k, trim(varname)

            ! Load test data
            call nc_read(filename_in,varname,var0)
            target_val = calc_grid_total(grid0%G%x,grid0%G%y,var,xlim=xlim,ylim=ylim)

            if (trim(varname) .eq. "mask") then
                ! Perform nearest neighbor interpolation

                var = interp_nearest(x=grid0%G%x,y=grid0%G%y,z=var0, &
                                      xout=grid%G%x,yout=grid%G%y)

            else 
                ! Perform conservative interpolation 

                call map_field_conservative_map1(map%map,varname,var0,var,missing_value=mv)

                current_val = calc_grid_total(grid%G%x,grid%G%y,var,xlim=xlim,ylim=ylim)
                err_percent = 100.d0 * (current_val-target_val) / target_val
                write(*,"(a,3g12.4)") "mass comparison (hi, con, % diff): ", &
                        target_val, current_val, err_percent                  
  
            end if 

            ! Write to file 
            call nc_write(filename,varname,var,dim1="xc",dim2="yc",missing_value=mv)
    
            ! Write variable metadata
            call nc_read_attr(filename_in,varname,"units",units)
            call nc_write_attr(filename,  varname,"units",units)
            call nc_read_attr(filename_in,varname,"long_name",long_name)
            call nc_write_attr(filename,  varname,"long_name",long_name)
            
            call nc_write_attr(filename,varname,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine rtopo_to_grid


end module rtopo 
