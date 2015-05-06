module ECMWF 

    use gridding_datasets
    use coordinates
    use interp_time 
    use ncio 
    
    implicit none 

    private 
    public :: ecmwf_to_grid
    public :: ecmwf_ocn_to_grid 

contains 

    subroutine ecmwf_to_grid(outfldr,grid,sigma,max_neighbors,lat_lim,clim_range)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       ECMWF DATA (ERA-INTERIM 1979-2013)
        !
        ! ========================================================= 

        implicit none 

        character(len=*)    :: outfldr 
        type(grid_class)    :: grid 
        double precision, optional :: sigma
        integer, optional   :: max_neighbors 
        double precision, optional :: lat_lim 
        integer, optional   :: clim_range(2) 
        character(len=512)  :: filename 
        character(len=512)  :: subfldr, cmd 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:), tmp(:)
            double precision, allocatable :: var(:,:)
        end type 

        type(inp_type)     :: inp

        type(grid_class)   :: gECMWF 
        character(len=256) :: file_invariant, file_surface, files_pres(9), file_precip
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:), precip(:)
        double precision, allocatable :: invar(:,:) 
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: nx, ny, nyr, nm, q, k, year, m, i, l 
        integer :: yearf, k0, nk 
        character(len=512) :: filename_clim 
        double precision, allocatable :: var3D(:,:,:), var2D(:,:)

        ! Assign the filenames
        file_invariant = "/data/sicopolis/data/ECMWF/ERA-INTERIM-invariant_historical_mon_197901-201212.nc"
        file_surface   = "/data/sicopolis/data/ECMWF/ERA-INTERIM-surface_historical_mon_197901-201212.nc"
        file_precip    = "/data/sicopolis/data/ECMWF/ERA-INTERIM-sf-tp_historical_mon_197901-201412.nc"
        files_pres(1)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-1000Mb_historical_mon_197901-201212.nc"
        files_pres(2)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-950Mb_historical_mon_197901-201212.nc"
        files_pres(3)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-850Mb_historical_mon_197901-201212.nc"
        files_pres(4)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-750Mb_historical_mon_197901-201212.nc"
        files_pres(5)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-700Mb_historical_mon_197901-201212.nc"
        files_pres(6)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-650Mb_historical_mon_197901-201212.nc"
        files_pres(7)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-600Mb_historical_mon_197901-201212.nc"
        files_pres(8)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-550Mb_historical_mon_197901-201212.nc"
        files_pres(9)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-500Mb_historical_mon_197901-201212.nc"

        desc    = "ERA-Interim dataset"
        ref     = "Dee et al., 2011, BAMS, &
                  &http://onlinelibrary.wiley.com/doi/10.1002/qj.828/abstract"
        
        ! Initialize the grid
        nx = nc_size(file_invariant,"longitude")
        ny = nc_size(file_invariant,"latitude")
        allocate(inp%lon(nx),inp%lat(ny))
        call nc_read(file_invariant,"longitude",inp%lon)
        call nc_read(file_invariant,"latitude", inp%lat)

        ! Flip latitudes because it is reversed in the file
        allocate(inp%tmp(ny))
        inp%tmp = inp%lat
        do i = 1, ny 
            inp%lat(i) = inp%tmp(ny-i+1)
        end do 
        deallocate(inp%tmp)

        call grid_init(gECMWF,name="ECMWF-075",mtype="latlon",units="kilometers",lon180=.TRUE., &
                       x=inp%lon,y=inp%lat)

        ! Add a subfolder to outfldr to hold all of the ECMWF files
        subfldr = "ERA-INT"
        cmd = "mkdir "//trim(outfldr)//"/"//trim(subfldr)
        call system(cmd)

        ! ## First make file for surface fields including invariants ##
        write(filename,"(a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                                trim(grid%name)//"_ERA-INT_197901-201212.nc"

        ! For climatology
        if (present(clim_range)) then  
            k0 = clim_range(1) - 1979 + 1
            nk = clim_range(2) - clim_range(1) + 1 

            ! Add a subfolder to outfldr to hold all of the ECMWF files
            subfldr = "ERA-INT"
            cmd = "mkdir "//trim(outfldr)//"_clim/"//trim(subfldr)
            call system(cmd)

            ! Make filename for the climatologies
            write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(subfldr)//"/"// &
                    trim(grid%name)//"_ERA-INT_",clim_range(1),"-",clim_range(2),".nc"
        end if 

        ! Define the pressure levels to be mapped
        plev = [1000,950,850,750,700,650,600,550,500]

        ! Define the variables to be mapped 
        
        allocate(invariant(1))
        call def_var_info(invariant(1),trim(file_invariant),"z","zs",units="m", &
                          long_name="Surface elevation",conv=1.d0/9.81d0)

        allocate(surf(12))
        call def_var_info(surf( 1),trim(file_surface),"sp", "sp", units="Pa",long_name="Surface pressure")
        call def_var_info(surf( 2),trim(file_surface),"tcw","tcw",units="kg m**-2", &
                          long_name="Total column water content")
        call def_var_info(surf( 3),trim(file_surface),"tclw","tclw",units="kg m**-2", &
                          long_name="Total column liquid water content")
        call def_var_info(surf( 4),trim(file_surface),"tciw","tciw",units="kg m**-2", &
                          long_name="Total column ice water content")
        call def_var_info(surf( 5),trim(file_surface),"p56.162","clw",units="kg m**-2", &
                          long_name="Total column liquid condensation")
        call def_var_info(surf( 6),trim(file_surface),"p57.162","ciw",units="kg m**-2", &
                          long_name="Total column ice condensation")
        call def_var_info(surf( 7),trim(file_surface),"tcc","tcc",units="(0 - 1)", &
                          long_name="Cloud cover fraction")
        call def_var_info(surf( 8),trim(file_surface),"u10","u10",units="m s**-1", &
                          long_name="10-m wind speed, u-component")
        call def_var_info(surf( 9),trim(file_surface),"v10","v10",units="m s**-1", &
                          long_name="10-m wind speed, v-component")
        call def_var_info(surf(10),trim(file_surface),"t2m","t2m",units="K", &
                          long_name="Near-surface temperature (2-m)")
        call def_var_info(surf(11),trim(file_surface),"al", "al", units="(0 - 1)", &
                          long_name="Surface albedo")
        call def_var_info(surf(12),trim(file_surface),"sst","sst",units="K", &
                          long_name="Sea surface temperature")

        allocate(precip(2))
        call def_var_info(precip( 1),trim(file_precip),"sf", "sf", units="kg m**-2 d**-1", &
                          long_name="Snowfall",conv=1d3)
        call def_var_info(precip( 2),trim(file_precip),"tp", "pr", units="kg m**-2 d**-1", &
                          long_name="Precipitation",conv=1d3)
        
        allocate(pres(5))
        call def_var_info(pres( 1),"None","t", "p_t",units="K",         plev="plev",filenames=files_pres, &
                          long_name="Temperature")
        call def_var_info(pres( 2),"None","z", "p_z",units="m**2 s**-2",plev="plev",filenames=files_pres, &
                          long_name="Geopotential height")
        call def_var_info(pres( 3),"None","u", "p_u",units="m s**-1",   plev="plev",filenames=files_pres, &
                          long_name="Wind speed, u-component")
        call def_var_info(pres( 4),"None","v", "p_v",units="m s**-1",   plev="plev",filenames=files_pres, &
                          long_name="Wind speed, v-component")
        call def_var_info(pres( 5),"None","w", "p_w",units="Pa s**-1",  plev="plev",filenames=files_pres, &
                          long_name="Vertical wind speed")

        nyr = 2012-1979+1
        nm  = 12 

        if (present(max_neighbors) .and. present(lat_lim)) then 

            ! Make sure Gaussian filter sigma available 
            if (.not. present(sigma)) then 
                write(*,*) "ecmwf_to_grid:: error:", "sigma must be provided for Gaussian filtering step."
                stop 
            end if 

            ! Allocate the input grid variable
            call grid_allocate(gECMWF,invar)

            ! Initialize mapping
            call map_init(map,gECMWF,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

            ! Initialize output variable arrays
            call grid_allocate(grid,outvar)
            call grid_allocate(grid,outmask)    

            ! Initialize the output file
            call nc_create(filename)
            call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
            call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
            call nc_write_dim(filename,"plev", x=dble(plev),units="hPa")
            call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
            call nc_write_dim(filename,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
            call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
            
            ! Write meta data 
            call nc_write_attr(filename,"Description",desc)
            call nc_write_attr(filename,"Reference",ref)

            ! ## INVARIANT FIELDS ##
            var_now = invariant(1) 
            call nc_read(trim(var_now%filename),var_now%nm_in,invar, &
                         start=[1,1,1],count=[gECMWF%G%nx,gECMWF%G%ny,1],missing_value=mv)
            call flip_lat(invar)
            invar = invar*var_now%conv 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,"nng",sigma=sigma,missing_value=mv)
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")


            ! ## SURFACE FIELDS ##
            do i = 1, size(surf)
                var_now = surf(i)

                q = 0 
                do k = 1, nyr 

                    year = 1978 + k 
                    write(*,*) trim(var_now%nm_in)," :",year

                    do m = 1, nm 
                        q = q+1 
                        call nc_read(trim(var_now%filename),var_now%nm_in,invar,start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1], &
                                     missing_value=mv)
                        call flip_lat(invar)
                        call map_field(map,var_now%nm_in,invar,outvar,outmask,"nng",sigma=sigma,missing_value=mv)
                        call nc_write(filename,var_now%nm_out,real(outvar),  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                      start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1],missing_value=real(mv))
                    end do 
                end do
                
                ! Write variable metadata
                call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
                call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
                
            end do  

            ! ## PRECIP FIELDS ##
            do i = 1, size(precip)
                var_now = precip(i)

                q = 0 
                do k = 1, nyr 

                    year = 1978 + k 
                    write(*,*) trim(var_now%nm_in)," :",year

                    do m = 1, nm 
                        q = q+1 
                        call nc_read(trim(var_now%filename),var_now%nm_in,invar,start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1], &
                                     missing_value=mv)
                        where (invar .ne. mv) invar = invar*var_now%conv 
                        call flip_lat(invar)
                        call map_field(map,var_now%nm_in,invar,outvar,outmask,"nng",sigma=sigma,missing_value=mv)
                        call nc_write(filename,var_now%nm_out,real(outvar),  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                      start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1],missing_value=real(mv))
                    end do 
                end do
                
                ! Write variable metadata
                call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
                call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
                
            end do  


            ! ## PRESSURE FIELDS ##
            do l = 1, size(files_pres)   ! Loop over pressure layers

                ! ## Make one file for each pressure level ##
                if (plev(l) .ge. 1000) then 
                    write(filename,"(a,i4,a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                                               trim(grid%name)//"_ERA-INT-",plev(l),"Mb_197901-201212.nc"
                else
                    write(filename,"(a,i3,a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                                               trim(grid%name)//"_ERA-INT-",plev(l),"Mb_197901-201212.nc"
                end if 

                ! Initialize the output file
                call nc_create(filename)
                call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
                call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
                call nc_write_dim(filename,"plev", x=dble(plev),units="hPa")
                call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
                call nc_write_dim(filename,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
                call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
                
                ! Write meta data 
                call nc_write_attr(filename,"Description",desc)
                call nc_write_attr(filename,"Reference",ref)

                ! ## PRESSURE FIELDS ##
                do i = 1, size(pres)
                    var_now = pres(i) 

                    q = 0 
                    do k = 1, nyr 

                        year = 1978 + k 
                        write(*,*) trim(var_now%nm_in)," :",year

                        do m = 1, nm 
                            q = q+1 
                            call nc_read(trim(var_now%filenames(l)),var_now%nm_in,invar, &
                                         start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1],missing_value=mv)
                            call flip_lat(invar)
                            call map_field(map,var_now%nm_in,invar,outvar,outmask,"nng",sigma=sigma,missing_value=mv)
                            call nc_write(filename,var_now%nm_out,real(outvar), &
                                          dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                          start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1],missing_value=real(mv))

                        end do 
                    end do 

                    ! Write variable metadata
                    call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
                    call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
                    call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
                
                end do 

            end do 

        end if 

        if (present(clim_range)) then 

            ! Create climatology too (month by month)

            call grid_allocate(grid,var2D)
            allocate(var3D(grid%G%nx,grid%G%ny,nk))    
            
            ! Initialize the output file
            call nc_create(filename_clim)
            call nc_write_dim(filename_clim,"xc",   x=grid%G%x,units="kilometers")
            call nc_write_dim(filename_clim,"yc",   x=grid%G%y,units="kilometers")
            call nc_write_dim(filename_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
            call grid_write(grid,filename_clim,xnm="xc",ynm="yc",create=.FALSE.)
            
            ! Write meta data 
            call nc_write_attr(filename_clim,"Description",desc)
            call nc_write_attr(filename_clim,"Reference",ref)

            ! ## INVARIANT FIELDS ##
            var_now = invariant(1) 
            call nc_read(filename,var_now%nm_out,var2D,missing_value=mv)
            call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",missing_value=real(mv))

            ! Write variable metadata
            call nc_write_attr(filename_clim,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename_clim,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename_clim,var_now%nm_out,"coordinates","lat2D lon2D")
            
            do i = 1, size(surf)
                var_now = surf(i)
                do m = 1, nm  
                    call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk], &
                                 missing_value=mv)
                    var2D = time_average(var3D)
                    call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                  start=[1,1,m],count=[grid%G%nx,grid%G%ny,1],missing_value=real(mv))
                end do 

                ! Write variable metadata
                call nc_write_attr(filename_clim,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename_clim,var_now%nm_out,"long_name",var_now%long_name)
                call nc_write_attr(filename_clim,var_now%nm_out,"coordinates","lat2D lon2D")
                
            end do 

            do i = 1, size(precip)
                var_now = precip(i)
                do m = 1, nm  
                    call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk], &
                                 missing_value=mv)
                    var2D = time_average(var3D)
                    call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                  start=[1,1,m],count=[grid%G%nx,grid%G%ny,1],missing_value=real(mv))
                end do 

                ! Write variable metadata
                call nc_write_attr(filename_clim,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename_clim,var_now%nm_out,"long_name",var_now%long_name)
                call nc_write_attr(filename_clim,var_now%nm_out,"coordinates","lat2D lon2D")
                
            end do 

            do l = 1, size(files_pres)   ! Loop over pressure layers

                ! ## Make one file for each pressure level ##
                if (plev(l) .ge. 1000) then 
                    write(filename,"(a,i4,a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                                               trim(grid%name)//"_ERA-INT-",plev(l),"Mb_197901-201212.nc"
                    write(filename_clim,"(a,i4,a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(subfldr)//"/"// &
                                               trim(grid%name)//"_ERA-INT-",plev(l),"Mb_",clim_range(1),"-",clim_range(2),".nc"
                else
                    write(filename,"(a,i3,a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                                               trim(grid%name)//"_ERA-INT-",plev(l),"Mb_197901-201212.nc"
                    write(filename_clim,"(a,i3,a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(subfldr)//"/"// &
                                               trim(grid%name)//"_ERA-INT-",plev(l),"Mb_",clim_range(1),"-",clim_range(2),".nc"
                end if 

                ! Initialize the output file
                call nc_create(filename_clim)
                call nc_write_dim(filename_clim,"xc",   x=grid%G%x,units="kilometers")
                call nc_write_dim(filename_clim,"yc",   x=grid%G%y,units="kilometers")
                call nc_write_dim(filename_clim,"plev", x=dble(plev),units="hPa")
                call nc_write_dim(filename_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
                call grid_write(grid,filename_clim,xnm="xc",ynm="yc",create=.FALSE.)
                
                ! Write meta data 
                call nc_write_attr(filename_clim,"Description",desc)
                call nc_write_attr(filename_clim,"Reference",ref)

                do i = 1, size(pres)
                    var_now = pres(i)

                    do m = 1, nm  
                        call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk], &
                                     missing_value=mv)
                        var2D = time_average(var3D)
                        call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                      start=[1,1,m],count=[grid%G%nx,grid%G%ny,1],missing_value=real(mv))
                    end do 

                    ! Write variable metadata
                    call nc_write_attr(filename_clim,var_now%nm_out,"units",var_now%units_out)
                    call nc_write_attr(filename_clim,var_now%nm_out,"long_name",var_now%long_name)
                    call nc_write_attr(filename_clim,var_now%nm_out,"coordinates","lat2D lon2D")
                
                end do 

            end do 

        end if 

        return 

    end subroutine ecmwf_to_grid


    subroutine ecmwf_ocn_to_grid(outfldr,grid,sigma,max_neighbors,lat_lim,clim_range)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       ECMWF ORAS4 oceanic DATA
        !
        ! =========================================================
        implicit none 

        character(len=*) :: outfldr
        type(grid_class) :: grid 
        integer, optional :: max_neighbors 
        double precision, optional :: sigma, lat_lim 
        integer, optional   :: clim_range(2) 
        character(len=512)  :: filename, subfldr, filename_clim 
        character(len=1024) :: desc, ref, cmd

        type inp_type 
            double precision, allocatable :: lon(:), lat(:), z_ocn(:), depth(:), time(:)
            double precision, allocatable :: var(:,:)
            integer,          allocatable :: mask(:,:)
        end type 

        type(inp_type)     :: inp
        type(grid_class)   :: grid0
        character(len=256) :: fldr_in, file_in, file_in_template 
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, nz, nt 

        type(map_class)  :: map
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), var3D(:,:,:), var2D(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var, d
        integer :: k0, k1, nk
        integer :: nm = 12
        character(len=4) :: year 

        ! Define the input filenames
        fldr_in          = "/data/sicopolis/data/ECMWF_ORAS4/"
        file_in_template = trim(fldr_in)//"{var}_oras4_1m_{year}_grid_1x1.nc"

        desc    = "ECMWF ORAS4 oceanic reanalysis output"
        ref     = "source folder: "//trim(fldr_in)

        ! Add a subfolder to outfldr to hold all of the ECMWF files
        subfldr = "ERA-INT-ORAS4"
        cmd = "mkdir "//trim(outfldr)//"/"//trim(subfldr)
        call system(cmd)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                              trim(grid%name)//"_ERA-INT-ORAS4_195801-201412.nc"

        ! Get correct filename for thetao, 2014
        file_in = trim(file_in_template)
        call replace(file_in,"{var}","thetao")
        call replace(file_in,"{year}","2014")
                    
        ! Load the domain information 
        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")
        nz = nc_size(file_in,"depth")
        allocate(inp%lon(nx),inp%lat(ny),inp%z_ocn(nz),inp%depth(nz))

        call nc_read(file_in,"lon",inp%lon)
        call nc_read(file_in,"lat",inp%lat)
        call nc_read(file_in,"depth",inp%depth)

        ! Make z negative and reverse it (rel to ocean surface)
        do k = 1, nz 
            inp%z_ocn(k) = -inp%depth(nz-k+1)
        end do  

        ! Define time info too 
        nt = 2014-1958+1 
        allocate(inp%time(nt))
        do k = 1, nt 
            inp%time(k) = 1958.d0 + (k-1)
        end do 
        
        ! Define and input grid variable fields
        call grid_init(grid0,name="ORAS4-1DEG",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat )
        call grid_allocate(grid0,inp%var)
        call grid_allocate(grid0,inp%mask)

        ! Define the variables to be mapped 
        allocate(vars(3))
        call def_var_info(vars( 1),trim(file_in_template),"thetao","to",units="degrees Celcius", &
                          long_name="Potential temperature",method="nng")
        call def_var_info(vars( 2),trim(file_in_template),"so","so",units="psu", &
                          long_name="Salinity",method="nng")
        call def_var_info(vars( 3),trim(file_in_template),"mask","mask_ocn",units="1", &
                          long_name="Land-ocean mask (0=land, 1=ocean)",method="nn")

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    

        if (present(sigma) .and. present(max_neighbors) .and. present(lat_lim)) then 

        ! Initialize mappings
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x, units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y, units="kilometers")
        call nc_write_dim(filename,"depth",x=inp%depth,units="meters")
        call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call nc_write_dim(filename,"time", x=inp%time, units="year")

        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ## Map variables ##
        
        ! Loop over each variable (except for mask)
        do q = 1, 2 

            var_now = vars(q)

            ! Map variables for each year and month and depth
            do k = 1, nt
            do m = 1, nm 
            do d = 1, nz 

                write(year,"(i4)") int(inp%time(k)) 
                if (k.eq.1) write(*,*) "YEAR, MONTH = ", year, m 

                ! Get correct filename
                file_in = trim(var_now%filename)
                call replace(file_in,"{var}",trim(var_now%nm_in))
                call replace(file_in,"{year}",trim(year))
                
                ! Read in current variable
                call nc_read(file_in,var_now%nm_in,inp%var,missing_value=mv, &
                             start=[1,1,d,m],count=[nx,ny,1,1])
                where(abs(inp%var) .ge. 1d10) inp%var = mv 

                ! Map the 2D field
                call map_field(map, var_now%nm_in,inp%var,outvar,outmask,"nng", &
                               fill=.TRUE.,missing_value=mv,sigma=sigma)

                ! Clean up infinite values or all missing layers
                ! (eg, for deep bathymetry levels for GRL domain)
                if (trim(var_now%nm_out) .eq. "to") then
                    where(outvar .ne. outvar .or. &
                        count(outvar.eq.mv) .eq. grid%npts) outvar = 1.d0
                end if 
                if (trim(var_now%nm_out) .eq. "so") then
                    where(outvar .ne. outvar .or. &
                        count(outvar.eq.mv) .eq. grid%npts) outvar = 35.d0
                end if 
                
                ! Write output variable to output file
                call nc_write(filename,var_now%nm_out,real(outvar), &
                              dim1="xc",dim2="yc",dim3="depth",dim4="month",dim5="time", &
                              start=[1,1,d,m,k],count=[grid%G%nx,grid%G%ny,1,1,1])

            end do 
            end do 
            end do 

            call nc_write_attr(filename,var_now%nm_out,"units",    var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        end do 

        ! Also map mask
        ! Map variable for each depth level

        var_now = vars(3) 
        ! note: filename taken from last file processed above 

        do d = 1, nz 
 
            ! Read in mask 
            call nc_read(file_in,var_now%nm_in,inp%mask,missing_value=int(mv), &
                         start=[1,1,1,d],count=[nx,ny,1,1])
            where(inp%mask == mv) inp%mask = 0 

            ! Map the 2D variable
            call map_field(map,var_now%nm_in,dble(inp%mask),outvar,outmask,method="nn", &
                          fill=.TRUE.,missing_value=mv)

            ! Write output mask to output file
            call nc_write(filename,var_now%nm_out,int(outvar), &
                          dim1="xc",dim2="yc",dim3="depth", &
                          start=[1,1,d],count=[grid%G%nx,grid%G%ny,1])
        
        end do 

        ! Write variable metadata
        call nc_write_attr(filename,var_now%nm_out,"units",    var_now%units_out)
        call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
        call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        end if 


        if (present(clim_range)) then 

            ! For climatology
            k0 = clim_range(1) - 1958 + 1
            nk = clim_range(2) - clim_range(1) + 1 

            ! Add a subfolder to outfldr to hold all of the ECMWF files
            subfldr = "ERA-INT"
            cmd = "mkdir "//trim(outfldr)//"_clim/"//trim(subfldr)
            call system(cmd)

            ! Make filename for the climatologies
            write(filename_clim,"(a,i4,a1,i4,a3)")  &
                    trim(outfldr)//"_clim/"//trim(subfldr)//"/"// &
                    trim(grid%name)//"_ERA-INT-ORAS4_",clim_range(1),"-",clim_range(2),".nc"


            ! Create climatology too (month by month)

            call grid_allocate(grid,var2D)
            allocate(var3D(grid%G%nx,grid%G%ny,nk))    
            
            ! Initialize the output file
            call nc_create(filename_clim)
            call nc_write_dim(filename_clim,"xc",   x=grid%G%x,units="kilometers")
            call nc_write_dim(filename_clim,"yc",   x=grid%G%y,units="kilometers")
            call nc_write_dim(filename_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
            call nc_write_dim(filename,"depth",x=inp%depth,units="meters")
            call grid_write(grid,filename_clim,xnm="xc",ynm="yc",create=.FALSE.)
            
            ! Write meta data 
            call nc_write_attr(filename_clim,"Description",desc)
            call nc_write_attr(filename_clim,"Reference",ref)

            do i = 1, 2
                var_now = vars(i)
                do d = 1, nz 
                do m = 1, nm  
                    call nc_read(filename,var_now%nm_out,var3D,start=[1,1,d,m,k0],count=[grid%G%nx,grid%G%ny,1,1,nk], &
                                 missing_value=mv)
                    var2D = time_average(var3D)
                    call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="depth",dim4="month", &
                                  start=[1,1,d,m],count=[grid%G%nx,grid%G%ny,1,1],missing_value=real(mv))
                end do 
                end do 

                ! Write variable metadata
                call nc_write_attr(filename_clim,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename_clim,var_now%nm_out,"long_name",var_now%long_name)
                call nc_write_attr(filename_clim,var_now%nm_out,"coordinates","lat2D lon2D")
                
            end do 

            ! Mask too 
            var_now = vars(3)
            do d = 1, nz 
                call nc_read(filename,var_now%nm_out,var2D,start=[1,1,d],count=[grid%G%nx,grid%G%ny,1], &
                                     missing_value=mv)
                call nc_write(filename_clim,var_now%nm_out,int(var2D),dim1="xc",dim2="yc",dim3="depth", &
                              start=[1,1,d],count=[grid%G%nx,grid%G%ny,1],missing_value=int(mv))
            end do 
        end if 


        return 

    end subroutine ecmwf_ocn_to_grid



    subroutine flip_lat(var)

        implicit none 
        double precision :: var(:,:)
        integer :: i, ny 
        double precision :: tmp(size(var,1),size(var,2))

        ny = size(var,2)
        
        tmp = var 

        do i = 1, ny 
            var(:,i) = tmp(:,ny-i+1)
        end do 

        return 

    end subroutine flip_lat 

end module ECMWF 
