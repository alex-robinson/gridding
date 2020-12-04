module rtopo 

    use coord
    use ncio
    use gridding_datasets
    use generic 
    use control 

    use gaussian_filter 

    implicit none 

    private 
    public :: rtopo_latlon_to_grid 
    public :: rtopo_to_grid 

contains 
    
    subroutine rtopo_latlon_to_grid_cdo(outfldr,grid,domain)
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
        
        integer :: nx1, ny1 
        real(4), allocatable :: var1(:,:)
        real(4), allocatable :: var2(:,:)
        character(len=56) :: grid_in_name 

        type(map_scrip_class) :: mps 

        ! ### Define input information #####
        path = "/p/projects/megarun/greenrise/datasets/rtopo-2"

        grid_in_name = "RTopo-2.0.1_30sec"

        ! 1. Bedrock topography ======================

        filename_in = trim(path)//"/"//trim(grid_in_name)//"_bedrock_topography.nc"

        nx1 = nc_size(filename_in,"lon")
        ny1 = nc_size(filename_in,"lat")
        

        ! Assume SCRIP map is already written for RTOPO file 
        ! (done manually to avoid generating a huge grid object)
        !call grid_write_cdo_desc_short(grid_in,fldr="maps") 
        
        ! Define output grid in grid description file 
        call grid_write_cdo_desc_short(grid,fldr="maps") 
        
        ! Generate SCRIP interpolation weights 
        call map_scrip_init(mps,grid_in_name,grid%name,fldr="maps",src_nc=filename_in)

        stop 

        desc    = "RTOPO-2.0.1 present-day Earth topography data"
        ref     = "Schaffer, J., Timmermann, R., Arndt, J. E., Kristensen, S. S., Mayer, C., &
                  &Morlighem, M., and Steinhage, D.: A global, high-resolution data set of &
                  &ice sheet topography, cavity geometry, and ocean bathymetry, &
                  &Earth Syst. Sci. Data, 8, 543-557, &
                  *https://doi.org/10.5194/essd-8-543-2016, 2016. &
                  &Data download: https://doi.pangaea.de/10.1594/PANGAEA.856844"

        ! ### Define output information #####
        
        ! Allocate output variable
        call grid_allocate(grid,var1)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_TOPO-RTOPO-2.0.1-new.nc"

        write(*,*) "Output filename: "//trim(filename)
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)


        return 

    end subroutine rtopo_latlon_to_grid_cdo


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
        inp%i1 = min(k,size(inp%lon))

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
        inp%j1 = min(k,size(inp%lat))

        inp%ni = inp%i1-inp%i0+1
        inp%nj = inp%j1-inp%j0+1 

        ! Allocate input array to the size of the reduced matrix 
        allocate(inp%var(inp%ni,inp%nj))
        
        write(*,*) "range(lon): ", minval(inp%lon(inp%i0:inp%i1)), maxval(inp%lon(inp%i0:inp%i1))
        write(*,*) "range(lat): ", minval(inp%lat(inp%j0:inp%j1)), maxval(inp%lat(inp%j0:inp%j1))
        write(*,*) "ni, nj:     ", inp%ni, inp%nj 

        desc    = "RTOPO-2.0.1 present-day Earth topography data"
        ref     = "Schaffer, J., Timmermann, R., Arndt, J. E., Kristensen, S. S., Mayer, C., &
                  &Morlighem, M., and Steinhage, D.: A global, high-resolution data set of &
                  &ice sheet topography, cavity geometry, and ocean bathymetry, &
                  &Earth Syst. Sci. Data, 8, 543-557, &
                  *https://doi.org/10.5194/essd-8-543-2016, 2016. &
                  &Data download: https://doi.pangaea.de/10.1594/PANGAEA.856844"

        ! ### Define output information #####
        
        ! Allocate output variable
        call grid_allocate(grid,var)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_TOPO-RTOPO-2.0.1.nc"

        write(*,*) "Output filename: "//trim(filename)
        
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
                               max_dist=1e3,lat_lim=0.02)
        
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
        real(8), allocatable :: var0(:,:), var0b(:,:), var(:,:)
        real(8) :: xlim(2), ylim(2) 

        real(8), allocatable :: z_base(:,:), z_srf(:,:), z_bed(:,:), H_ice(:,:)   
        integer, allocatable :: mask(:,:) 
        real(8) :: sigma 

        real(8), parameter :: rho_ice =  910.0 
        real(8), parameter :: rho_sw  = 1028.0 
        
        ! ### Input information #####

        ! Define input grid 
        select case(trim(domain))
            case("North")
                call domain_definition(grid0,"NH-4KM")
            case("Laurentide")
                call domain_definition(grid0,"LIS-4KM")
            case("Eurasia")
                call domain_definition(grid0,"EIS-4KM")
            case("Greenland")
                call domain_definition(grid0,"GRL-1KM")
            case("Antarctica")
                call domain_definition(grid0,"ANT-4KM")
            case DEFAULT 
                write(*,*) "rtopo_to_grid:: error: domain not recognized: "//trim(domain)
                stop

        end select 

        xlim = [minval(grid0%G%x), maxval(grid0%G%x)]
        ylim = [minval(grid0%G%y), maxval(grid0%G%y)]

        call grid_allocate(grid0,var0) 
        call grid_allocate(grid0,var0b) 
        
        path = "output/"//trim(domain)
        filename_in = trim(path)//"/"//trim(grid0%name)//"/"// &
                        trim(grid0%name)//"_TOPO-RTOPO-2.0.1.nc"

        desc    = "RTOPO-2.0.1 present-day Earth topography data"
        ref     = "Schaffer, J., Timmermann, R., Arndt, J. E., Kristensen, S. S., Mayer, C., &
                  &Morlighem, M., and Steinhage, D.: A global, high-resolution data set of &
                  &ice sheet topography, cavity geometry, and ocean bathymetry, &
                  &Earth Syst. Sci. Data, 8, 543-557, &
                  *https://doi.org/10.5194/essd-8-543-2016, 2016. &
                  &Data download: https://doi.pangaea.de/10.1594/PANGAEA.856844"

        ! ### Output information ##### 

        call grid_allocate(grid,var) 
        
        ! Initialize the map 
        call map_init(map,grid0,grid,max_neighbors=20,lat_lim=0.5d0,fldr="maps",load=.TRUE.)
        
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

        allocate(varnames(3))
        varnames = ["z_bed     ","z_srf     ","z_ice_base","H_ice     "]

        do k = 1, size(varnames)

            varname = varnames(k)
            write(*,*) k, trim(varname)

            ! Load test data
            if (trim(varname) .ne. "H_ice") then 
                call nc_read(filename_in,varname,var0,missing_value=mv)
            else 
                call nc_read(filename_in,"z_srf",var0,missing_value=mv)
                call nc_read(filename_in,"z_ice_base",var0b,missing_value=mv)
                var0 = max(var0 - var0b,0.0)
                where(var0b .eq. 0.0) var0 = 0.0  
            end if 

            target_val = calc_grid_total(grid0%G%x,grid0%G%y,var0,xlim=xlim,ylim=ylim)

            write(*,*) "range(var_in): ", minval(var0,mask=var0.ne.mv), maxval(var0,mask=var0.ne.mv)

if (.FALSE.) then 
            ! Smooth data on hires grid
            call filter_gaussian(var=var0,sigma=grid%G%dx/2.d0,dx=grid0%G%dx)

            ! Perform nearest neighbor interpolation to get lores output 
            var = mv 
            var = interp_nearest(x=grid0%G%x,y=grid0%G%y,z=var0, &
                                 xout=grid%G%x,yout=grid%G%y)
else 
            ! Perform conservative interpolation 
            var = mv 
            call map_field_conservative_map1(map%map,varname,var0,var,method="mean",missing_value=mv)

end if 

            current_val = calc_grid_total(grid%G%x,grid%G%y,var,xlim=xlim,ylim=ylim)
            err_percent = 100.d0 * (current_val-target_val) / target_val
            write(*,"(a,3g12.4)") "mass comparison (hi, con, % diff): ", &
                    target_val, current_val, err_percent                  

            write(*,*) "range(var_out): ", minval(var,mask=var.ne.mv), maxval(var,mask=var.ne.mv)

            ! Write to file
            call nc_write(filename,varname,real(var),dim1="xc",dim2="yc",missing_value=real(mv))

            ! Write variable metadata
            if (trim(varname) .ne. "H_ice") then
                call nc_read_attr(filename_in,varname,"units",units)
                call nc_write_attr(filename,  varname,"units",units)
                call nc_read_attr(filename_in,varname,"long_name",long_name)
                call nc_write_attr(filename,  varname,"long_name",long_name)
                call nc_write_attr(filename,varname,"coordinates","lat2D lon2D")
            else 
                ! Write variable metadata
                call nc_write_attr(filename,varname,"units","m")
                call nc_write_attr(filename,varname,"long_name","Ice thickness")
                call nc_write_attr(filename,varname,"coordinates","lat2D lon2D")
            end if 
            
        end do 

        ! Calculate additional variables of interest 
        call grid_allocate(grid,z_srf)
        call grid_allocate(grid,z_bed)
        call grid_allocate(grid,z_base)
        call grid_allocate(grid,H_ice)
        call grid_allocate(grid,mask)

        call nc_read(filename,"z_ice_base",z_base)
        call nc_read(filename,"z_srf",z_srf)
        call nc_read(filename,"z_bed",z_bed)
        call nc_read(filename,"H_ice",H_ice)
        
        ! Delete grounded points with small ice thickness 
        ! Set threshold for "floating" to 100m difference with bed, in order to avoid 
        ! artifacts 
        where (H_ice .lt. 20.0 .and. abs((z_srf-H_ice)-z_bed) .lt. 100e0)
            H_ice = 0.0 
            z_srf = max(z_bed,0.0)
        end where 

        ! Delete floating points with small ice thickness 
        ! Set threshold for "floating" to 100m difference with bed, in order to avoid 
        ! artifacts 
        where (H_ice .lt. 100.0 .and. abs((z_srf-H_ice)-z_bed) .ge. 100e0)
            H_ice = 0.0 
            z_srf = 0.0 
        end where 

        ! Now generate the consistent mask 
        where (z_srf .gt. 0.0 .and. H_ice .gt. 0.0 .and. abs((z_srf-H_ice)-z_bed) .lt. 10e0) 
            ! Ice thickness touches bedrock to within 1m, then it is grounded ice 
            mask = 2.0
        else where (z_srf .gt. 0.0 .and. H_ice .gt. 0.0) 
            ! Floating ice 
            mask = 3.0 
        else where (z_srf .gt. 0.0) 
            ! Ice-free land 
            mask = 1.0
        elsewhere
            ! Ice-free ocean 
            mask = 0.0  
        end where 

        ! Rewrite surface elevation to file after modification 
        varname = "z_srf"
        call nc_write(filename,varname,real(z_srf),dim1="xc",dim2="yc",missing_value=real(mv))

        ! Write new variable to file: Ice thickness 
        varname = "H_ice"
        call nc_write(filename,varname,real(H_ice),dim1="xc",dim2="yc",missing_value=real(mv))

        ! Write new variable to file: ice mask  
        varname = "mask"
        call nc_write(filename,varname,mask,dim1="xc",dim2="yc",missing_value=int(mv))

        ! Write variable metadata
        call nc_write_attr(filename,varname,"units","")
        call nc_write_attr(filename,varname, &
                "long_name","Mask (0:ocean, 1:land, 2:grounded ice, 3:floating ice)")
        call nc_write_attr(filename,varname,"coordinates","lat2D lon2D")
        
        ! For bedrock elevation, additionally calculate the standard deviation of the field 
        ! var0 == high resolution field 
        ! var  == destination field 

        call nc_read(filename_in,"z_bed",var0,missing_value=mv)

        var = mv 

        call map_field_conservative_map1(map%map,"z_bed",var0,var,method="stdev",missing_value=mv)

        sigma = grid%G%dx
        mask = 0
        where(var.eq.mv) mask = 1 
        call fill_weighted(var,missing_value=mv)
        call filter_gaussian(var=var,sigma=sigma,dx=grid%G%dx,mask=mask.eq.1)

        call nc_write(filename,"z_bed_sd",real(var),dim1="xc",dim2="yc",missing_value=real(mv))
        
        ! Write variable metadata
        call nc_write_attr(filename,"z_bed_sd","units","m")
        call nc_write_attr(filename,"z_bed_sd","long_name","Bedrock standard deviation")
        call nc_write_attr(filename,"z_bed_sd","coordinates","lat2D lon2D")

        return 

    end subroutine rtopo_to_grid


end module rtopo 
