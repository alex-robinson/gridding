module RACMO2 

    use gridding_datasets
    use coordinates 
    use ncio 
    
    implicit none 

    private 
    public :: RACMO2rot_to_grid
    
contains 

    subroutine RACMO2rot_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,clim_range)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       RACMO2/ANT (RCM) DATA - RACMO2 data obtained from
        !       Willem van de Berg (per. comm.). Available for 
        !       HadCM3 driven A1B simulation (2000-2199)
        !
        ! =========================================================

        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer, optional :: max_neighbors 
        double precision, optional :: lat_lim 
        integer, optional :: clim_range(2)

        character(len=512) :: filename 
        character(len=512) :: fldr_input, file_suffix1, file_suffix2
        type(points_class) :: pIN1, pIN2, pIN3
        type(var_defs), allocatable :: vars0(:), vars(:)
        double precision, allocatable :: invar(:), lon(:), lat(:)
        integer, allocatable :: invar_int(:) 
        integer :: nx, ny 

        type(map_class)  :: map
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:), outvar_int(:,:)

        integer :: nyr, nm, q, k, year, m, i, l, year0, n_prefix, n_var 
        integer :: yearf, k0, nk 
        character(len=512) :: filename_clim 
        double precision, allocatable :: var3D(:,:,:), var2D(:,:)
        double precision :: conv_tosec 

        ! Define input grid
        if (trim(domain) .eq. "Antarctica-A1B") then 
            
            ! Define the input filenames
            fldr_input     = "tmpdata/"
            file_suffix1   = "_RACMO2_ANT3K55_HadCM3-A1B.nc"
            file_suffix2   = "_RACMO2_ANT3K55_HadCM3-A1B_2000-2199.nc"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_RACMO2-ANT3K55_HadCM3-A1B-monthly_2000-2199.nc"

            year0       = 2000 
            nyr         = 2199-2000+1

            ! For climatology
            if (present(clim_range)) then  
                k0 = clim_range(1) - year0+1
                nk = clim_range(2) - clim_range(1) + 1 

                write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
                    "_RACMO2-ANT3K55_HadCM3-A1B-monthly__",clim_range(1),"-",clim_range(2),".nc"
            end if 

        else if (trim(domain) .eq. "Antarctica-c20") then 

            ! Define the input filenames
            fldr_input     = "tmpdata2/"
            file_suffix1   = "_RACMO2_ANT3K55_HadCM3-c20.nc"
            file_suffix2   = "_RACMO2_ANT3K55_HadCM3-c20_1980-1999.nc"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_RACMO2-ANT3K55_HadCM3-c20-monthly_1980-1999.nc"

            year0       = 1980 
            nyr         = 1999-1980+1

            ! For climatology
            if (present(clim_range)) then  
                k0 = clim_range(1) - year0+1
                nk = clim_range(2) - clim_range(1) + 1 

                write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
                    "_RACMO2-ANT3K55_HadCM3-c20-monthly_",clim_range(1),"-",clim_range(2),".nc"
            end if 




        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define RACMO2 input grids/points ===========
        
        nx = 134 
        ny = 122
        
        if (allocated(lon)) deallocate(lon)
        if (allocated(lat)) deallocate(lat)
        allocate(lon(nx*ny),lat(nx*ny))
        call nc_read(trim(fldr_input)//"Geopotential"//trim(file_suffix1),"g10_lon_1",lon,start=[1,1],count=[nx,ny])
        call nc_read(trim(fldr_input)//"Geopotential"//trim(file_suffix1),"g10_lat_0",lat,start=[1,1],count=[nx,ny])
        call points_init(pIN1,name="ANT3K55",mtype="latlon",units="degrees",lon180=.TRUE., &
                         x=lon,y=lat)

        ! Define the variables to be mapped

        ! ## INVARIANT (2D) FIELDS ## 
        allocate(vars0(3))
        call def_var_info(vars0(1),trim(fldr_input)//"Geopotential"//trim(file_suffix1), &
                          "GP_GDS10_HTGL_ave1h", "zs",  units="m",fill=.TRUE.,conv=1.d0/9.81d0)
        call def_var_info(vars0(2),trim(fldr_input)//"IceMask"//trim(file_suffix1), &
                          "ICE_C_GDS10_HTGL_ave1h", "mask_ice",  units="-",fill=.TRUE.,method="nn")
        call def_var_info(vars0(3),trim(fldr_input)//"LSM"//trim(file_suffix1), &
                          "LAND_GDS10_HTGL_ave1h", "mask_land",  units="-",fill=.TRUE.,method="nn")


        ! ## SURFACE (3D) FIELDS ##

        conv_tosec = 1.d0/(30.d0*86400.d0)

        if (allocated(vars)) deallocate(vars)
        allocate(vars(21))
        call def_var_info(vars(1),trim(fldr_input)//"t2m"//trim(file_suffix2), &
                          "t2m", "t2m",  units="K",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(2),trim(fldr_input)//"clcov"//trim(file_suffix2), &
                          "clcov", "cc",  units="-",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(3),trim(fldr_input)//"evap"//trim(file_suffix2), &
                          "evap", "evap",  units="kg m**-2 d**-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(4),trim(fldr_input)//"precip"//trim(file_suffix2), &
                          "precip", "pr",  units="kg m**-2 d**-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(5),trim(fldr_input)//"q2m"//trim(file_suffix2), &
                          "q2m", "qs",  units="kg kg-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(6),trim(fldr_input)//"rain"//trim(file_suffix2), &
                          "rain", "rf",  units="kg m**-2 d**-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(7),trim(fldr_input)//"refreeze"//trim(file_suffix2), &
                          "refreeze", "rz",  units="kg m**-2 d**-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(8),trim(fldr_input)//"runoff"//trim(file_suffix2), &
                          "runoff", "ru",  units="kg m**-2 d**-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(9),trim(fldr_input)//"smb"//trim(file_suffix2), &
                          "smb", "smb",  units="kg m**-2 d**-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(10),trim(fldr_input)//"snowfall"//trim(file_suffix2), &
                          "snowfall", "sf",  units="kg m**-2 d**-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(11),trim(fldr_input)//"snowmelt"//trim(file_suffix2), &
                          "snowmelt", "me",  units="kg m**-2 d**-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(12),trim(fldr_input)//"sublim"//trim(file_suffix2), &
                          "sublim", "subl",  units="kg m**-2 d**-1",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(13),trim(fldr_input)//"tskin"//trim(file_suffix2), &
                          "tskin", "ts",  units="K",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(14),trim(fldr_input)//"u10m"//trim(file_suffix2), &
                          "u10m", "uas",  units="m s**-1 [rot]",fill=.TRUE.,dimextra=.TRUE.)
        call def_var_info(vars(15),trim(fldr_input)//"v10m"//trim(file_suffix2), &
                          "v10m", "vas",  units="m s**-1 [rot]",fill=.TRUE.,dimextra=.TRUE.)
        
        call def_var_info(vars(16),trim(fldr_input)//"LHF"//trim(file_suffix2), &
                "LHTFL_GDS10_HTGL_acc", "lhf",  units="W m**-2",fill=.TRUE.,conv=conv_tosec)
        call def_var_info(vars(17),trim(fldr_input)//"LWD"//trim(file_suffix2), &
                "VAR_177_GDS10_HTGL_acc", "lwd",  units="W m**-2",fill=.TRUE.,conv=conv_tosec)
        call def_var_info(vars(18),trim(fldr_input)//"LWN"//trim(file_suffix2), &
                "VAR_177_GDS10_HTGL_acc", "lwn",  units="W m**-2",fill=.TRUE.,conv=conv_tosec)
        call def_var_info(vars(19),trim(fldr_input)//"SWD"//trim(file_suffix2), &
                "VAR_176_GDS10_HTGL_acc", "swd",  units="W m**-2",fill=.TRUE.,conv=conv_tosec)
        call def_var_info(vars(20),trim(fldr_input)//"SWN"//trim(file_suffix2), &
                "VAR_176_GDS10_HTGL_acc", "swn",  units="W m**-2",fill=.TRUE.,conv=conv_tosec)
        call def_var_info(vars(21),trim(fldr_input)//"Albedo"//trim(file_suffix2), &
                "ALBDO_GDS10_HTGL_ave1h", "al",  units="-",fill=.TRUE.)
        
        if (trim(domain) .ne. "Antarctica-c20") then
            ! Rename albedo to load empty data, since field is not available for other scenarios
            call def_var_info(vars(21),trim(fldr_input)//"SWN"//trim(file_suffix2), &
                "VAR_176_GDS10_HTGL_acc", "al",  units="-",fill=.TRUE.,conv=0.d0)
        end if 
        
        nm       = 12
        n_var    = size(vars)

        if (present(max_neighbors) .and. present(lat_lim)) then 

            ! Allocate the input grid variable
            call points_allocate(pIN1,invar)
            call points_allocate(pIN1,invar_int)

            ! Initialize mapping
            call map_init(map,pIN1,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

            ! Initialize output variable arrays
            call grid_allocate(grid,outvar)
            call grid_allocate(grid,outmask)    
            call grid_allocate(grid,outvar_int)

            ! Initialize the output file
            call nc_create(filename)
            call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
            call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
            call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
            call nc_write_dim(filename,"time", x=year0,dx=1,nx=nyr,units="years",calendar="360_day")
            call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
            
            ! ## INVARIANT (2D) FIELDS ##
            do i = 1, size(vars0)

                var_now = vars0(i)
                call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value, &
                             start=[1,1],count=[nx,ny])
                outvar = missing_value
                call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,50.d3, &
                               fill=.TRUE.,missing_value=missing_value)
                if (var_now%fill) call fill_mean(outvar,missing_value=missing_value,fill_value=0.d0)
                if (trim(var_now%method) .eq. "nn") then 
                    call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc", &
                                  units=var_now%units_out,missing_value=nint(missing_value))
                else 
                    call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc", &
                                  units=var_now%units_out,missing_value=real(missing_value))
                end if 

            end do 

            ! ## SURFACE (3D) FIELDS ##
            do i = 1, size(vars)
                var_now = vars(i)     
                write(*,*) "=== ",trim(var_now%nm_out)," ==="
     
                do k = 1, nyr 
!                     write(*,*) year0-1+k
                    do m = 1, nm 
                        q = (k-1)*12 + m 

                        if (var_now%dimextra) then 
                            call nc_read(trim(var_now%filename),var_now%nm_in,invar, &
                                     missing_value=missing_value, &
                                     start=[1,1,1,q],count=[nx,ny,1,1])
                        else 
                            call nc_read(trim(var_now%filename),var_now%nm_in,invar,&
                                     missing_value=missing_value, &
                                     start=[1,1,q],count=[nx,ny,1])
                        end if
                        where (invar .ne. missing_value) invar = invar*var_now%conv 
                        outvar = missing_value 
                        call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",100.d3, &
                                       fill=.FALSE.,missing_value=missing_value)
                        call nc_write(filename,var_now%nm_out,real(outvar),units=var_now%units_out, &
                                      dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                      start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])
            
                    end do 
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
        
            ! ## INVARIANT (2D) FIELDS ##
            do i = 1, size(vars0)
                var_now = vars0(i) 
                call nc_read(filename,var_now%nm_out,var2D)
                call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc", &
                              units=var_now%units_out)
            end do 

            ! ## SURFACE (3D) FIELDS ##
            do i = 1, size(vars)
                var_now = vars(i)

                do m = 1, nm  
                    call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk])
                    var2D = time_average(var3D)
                    call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                  units=var_now%units_out,start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
                end do 
            end do 

        end if 

        return 

    end subroutine RACMO2rot_to_grid

end module RACMO2