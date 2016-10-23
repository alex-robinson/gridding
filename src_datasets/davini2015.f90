module davini2015 

    use gridding_datasets
    use coord
    use ncio 
    use gaussian_filter
    
    implicit none 

    private 
    public :: davini2015_to_grid

contains 

    subroutine davini2015_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       Davini et al. (2015, grl) atmospheric DATA
        !
        ! =========================================================
        implicit none 

        character(len=*) :: domain, outfldr, subfldr, path_in 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: sigma, lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:), var(:,:)
        end type 

        type(inp_type)     :: inp
        type(grid_class)   :: grid0
        character(len=256) :: fldr_in, file_in
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, np 

        type(map_class)  :: map
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        ! Define the input filenames
        fldr_in      = trim(path_in)
        file_in      = trim(fldr_in)//"davini_"//trim(domain)//".nc"

        write(*,*) "fldr_in: "//trim(fldr_in)
        write(*,*) "file_in: "//trim(file_in)

        desc    = "EC-EARTH simulation output"
        ref     = "Davini et al.: Impact of Greenland orography on the &
                  &Atlantic Meridional Overturning Circulation, &
                  &Geophysical Research Letters, doi:10.1002/2014GL062668, 2015."

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                              trim(grid%name)//"_davini_"//trim(domain)//".nc"

        ! Load the domain information 
        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")
        np = nx*ny 

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)
        
        ! Define points and input variable field
        call grid_init(grid0,name="davini-latlon",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat )

        ! Define the variables to be mapped 
        allocate(vars(11))
        call def_var_info(vars( 1),trim(file_in),"tas","t2m",units="Kelvin", &
                          long_name="Near-surface temperature (2-m)",method="nng")
        call def_var_info(vars( 2),trim(file_in),"pr","pr",units="mm*d**-1", &
                          long_name="Precipitation",method="nng")
        call def_var_info(vars( 3),trim(file_in),"al","al",units="1", &
                          long_name="Surface albedo",method="nng")
        call def_var_info(vars( 4),trim(file_in),"uas","uas",units="m*s**-1", &
                          long_name="Surface velocity, x-component",method="nng")
        call def_var_info(vars( 5),trim(file_in),"vas","vas",units="m*s**-1", &
                          long_name="Surface velocity, y-component",method="nng")
        call def_var_info(vars( 6),trim(file_in),"uvas","uvas",units="m*s**-1", &
                          long_name="Surface velocity, magnitude",method="nng")
        call def_var_info(vars( 7),trim(file_in),"p700_u","p700_u",units="m*s**-1", &
                          long_name="700 Mb velocity, x-component",method="nng")
        call def_var_info(vars( 8),trim(file_in),"p700_v","p700_v",units="m*s**-1", &
                          long_name="700 Mb velocity, y-component",method="nng")
        call def_var_info(vars( 9),trim(file_in),"p700_uv","p700_uv",units="m*s**-1", &
                          long_name="700 Mb velocity, magnitude",method="nng")
        call def_var_info(vars(10),trim(file_in),"p700_Z","p700_Z",units="m**2*s**-2", &
                          long_name="700 Mb geopotential",method="nng")
        call def_var_info(vars(11),trim(file_in),"p700_z","p700_z",units="m", &
                          long_name="700 Mb geopotential height",method="nng")

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"month",x=1,dx=1,nx=12,units="month")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Load reference topography in order to adjust temps to sea-level temps 
        call nc_read(file_in,"zs",inp%var,missing_value=missing_value)
        write(*,*) "input zs : ", minval(inp%var), maxval(inp%var)

        ! Map zs to new grid
        call map_field(map,"zs",inp%var,outvar,outmask,"nng", &
                          fill=.TRUE.,missing_value=missing_value,sigma=sigma)
        call nc_write(filename,"zs",real(outvar),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,"zs","units","m")
        call nc_write_attr(filename,"zs","long_name","Surface elevation")
        call nc_write_attr(filename,"zs","coordinates","lat2D lon2D")
        
        ! Also generate a land mask
        where(inp%var .lt. 10.d0) inp%var = 0.d0 
        where(inp%var .gt.  0.d0) inp%var = 1.d0 

        ! Map mask to new grid
        call map_field(map,"mask_land",inp%var,outvar,outmask,"nn", &
                          fill=.TRUE.,missing_value=missing_value)
        call nc_write(filename,"mask_land",nint(outvar),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,"mask_land","units","1")
        call nc_write_attr(filename,"mask_land","long_name","Land mask (land=1)")
        call nc_write_attr(filename,"mask_land","coordinates","lat2D lon2D")

        ! ## Map climatological gridded variables ##
        
        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            do k = 1, 12 

                ! Read in current variable
                call nc_read(trim(var_now%filename),var_now%nm_in,inp%var, &
                             start=[1,1,k],count=[nx,ny,1], &
                             missing_value=missing_value)
                where(abs(inp%var) .ge. 1d10) inp%var = missing_value 

                ! Map variable to new grid
                call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
                              fill=.TRUE.,missing_value=missing_value,sigma=sigma)
                
                ! Write output variable to output file
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month", &
                              start=[1,1,k],count=[grid%G%nx,grid%G%ny,1])

            end do 

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
        
        end do 

        return 

    end subroutine davini2015_to_grid

end module davini2015

