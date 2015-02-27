module MAR 

    use gridding_datasets
    use coordinates
    use interp2D
    use interp_time 
    use ncio 
    
    implicit none 

    private 
    public :: MARv35_to_grid
!     public :: MARv33_to_grid
!     public :: MARv32_to_grid

contains 

    subroutine MARv35_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,clim_range)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MAR (RCM) DATA - MARv3.5 downloaded from the ftp site:
        !       ftp://ftp.climato.be/fettweis/MARv3.5/Greenland
        !       * Routine expects the -RAW- MAR data on eg, 30 km grid
        !       Note: data is also available on the Bamber et al. (2001) 5km grid
        !       domain="Greenland-ERA": ERA-40 + ERA-Interim combined datasets
        !
        ! =========================================================

        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer, optional :: max_neighbors 
        double precision, optional :: lat_lim 
        integer, optional :: clim_range(2)

        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type(grid_class)   :: gMAR
        character(len=256) :: file_invariant, file_surface, file_prefix(2)
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:)
        integer, allocatable :: invar_int(:,:) 
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:), outvar_int(:,:)

        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 
        integer :: yearf, k0, nk 
        character(len=512) :: filename_clim 
        double precision, allocatable :: var3D(:,:,:), var2D(:,:)

        ! Define input grid
        if (trim(domain) .eq. "Greenland-ERA") then 
            
            ! Define MAR raw grid and input variable field
            call grid_init(gMAR,name="MAR-30KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                           x0=-780.d0,dx=30.d0,nx=51,y0=-1230.d0,dy=30.d0,ny=93, &
                           lambda=-40.d0,phi=71.d0,alpha=7.5d0)

            ! Define the input filenames
            file_invariant = "/data/sicopolis/data/MARv3.5/Greenland/ERA_1958-2013_30km-rawclean/"// &
                             "MARv3.5-ERA-30km-monthly-2013.nc"
            file_surface   = "/data/sicopolis/data/MARv3.5/Greenland/"
            file_prefix(1) = "ERA_1958-2013_30km-rawclean/MARv3.5-ERA-30km-monthly-"
            file_prefix(2) = "ERA_1958-2013_30km-rawclean/MARv3.5-ERA-30km-monthly-"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_MARv3.5-ERA-30km-monthly_195801-201312.nc"
            desc    = "Greenland regional climate simulated by MARv3.5"
            ref     = "Fettweis et al., ftp.climato.be/fettweis/MARv3.5/Greenland"

            year0       = 1958 
            year_switch = 1979   ! Switch scenarios (ERA-40 to ERA-INTERIM)
            nyr         = 2013-1958+1

            ! For climatology
            if (present(clim_range)) then  
                k0 = clim_range(1) - year0+1
                nk = clim_range(2) - clim_range(1) + 1 

                write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
                    "_MARv3.5-ERA-30km-monthly_",clim_range(1),"-",clim_range(2),".nc"
            end if 

        else if (trim(domain) .eq. "Greenland-MIROC5-RCP85") then 

            ! Define MAR (Bamber et al. 2001) grid and input variable field
            call grid_init(gMAR,name="Bamber01-5KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                           x0=-800.d0,dx=5.d0,nx=301,y0=-3400.d0,dy=5.d0,ny=561, &
                           lambda=-39.d0,phi=90.d0,alpha=7.5d0)

            ! Define the input filenames
            file_invariant = "/data/sicopolis/data/MARv3.5/Greenland/MIROC5-histo_1976-2005_30km/"// &
                             "MARv3.5-monthly-MIROC5-histo-1976.nc"
            file_surface   = "/data/sicopolis/data/MARv3.5/Greenland/"
            file_prefix(1) = "MIROC5-histo_1976-2005_30km/MARv3.5-monthly-MIROC5-histo-"
            file_prefix(2) = "MIROC5-rcp85_2006-2100_30km/MARv3.5-monthly-MIROC5-rcp85-"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_MARv3.5-30km-monthly-MIROC5-rcp85_197601-210012.nc"
            desc    = "Greenland regional climate simulated by MARv3.5"
            ref     = "Fettweis et al., ftp.climato.be/fettweis/MARv3.5/Greenland"

            year0       = 1976
            year_switch = 2006   ! Switch scenarios (historical to RCP85)
            nyr         = 2100-1976+1

            ! For climatology
            if (present(clim_range)) then  
                k0 = clim_range(1) - year0+1
                nk = clim_range(2) - clim_range(1) + 1 

                write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
                    "_MARv3.5-30km-monthly-MIROC5-rcp85_",clim_range(1),"-",clim_range(2),".nc"
            end if 

        else if (trim(domain) .eq. "Antarctica") then 

            ! TODO 

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(3))
        call def_var_info(invariant(1),trim(file_invariant),"SH", "zs",  units="m",fill=.TRUE., &
            long_name="Surface elevation")
        call def_var_info(invariant(2),trim(file_invariant),"SRF","mask",units="(0 - 4)", &
            long_name="Land ice mask",method="nn",fill=.TRUE.)
        call def_var_info(invariant(3),trim(file_invariant),"MSK","msk", units="%", &
            long_name="Cell ice coverage",fill=.TRUE.)
        
        allocate(surf(28))
        call def_var_info(surf( 1),trim(file_surface),"SHSN0", "SH0", units="m", &
            long_name="Snow height (layer 0)")  
        call def_var_info(surf( 2),trim(file_surface),"SHSN2", "SH2", units="m", &
            long_name="Snow height (layer 2)")  
        call def_var_info(surf( 3),trim(file_surface),"SHSN3", "SH3", units="m", &
            long_name="Snow height (layer 3)")
        call def_var_info(surf( 4),trim(file_surface),"SMB", "smb", units="mm d**-1", &
            long_name="Surface mass balance") 
        call def_var_info(surf( 5),trim(file_surface),"SU",  "su",  units="mm d**-1", &
            long_name="Sublimation") 
        call def_var_info(surf( 6),trim(file_surface),"ME",  "me",  units="mm d**-1", &
            long_name="Total melt",fill=.TRUE.) 
        call def_var_info(surf( 7),trim(file_surface),"RZ",  "rz",  units="mm d**-1", &
            long_name="Total refreezing",fill=.TRUE.) 
        call def_var_info(surf( 8),trim(file_surface),"SF",  "sf",  units="mm d**-1", &
            long_name="Snowfall",fill=.TRUE.) 
        call def_var_info(surf( 9),trim(file_surface),"RF",  "rf",  units="mm d**-1", &
            long_name="Rainfall",fill=.TRUE.) 
        call def_var_info(surf(10),trim(file_surface),"RU",  "ru",  units="mm d**-1", &
            long_name="Runoff") 
        call def_var_info(surf(11),trim(file_surface),"UU",  "u",   units="m s**-1", &
            long_name="Surface wind, u-component",fill=.TRUE.)
        call def_var_info(surf(12),trim(file_surface),"VV",  "v",   units="m s**-1", &
            long_name="Surface wind, v-component",fill=.TRUE.)
        call def_var_info(surf(13),trim(file_surface),"TT",  "t3m", units="degrees Celcius", &
            long_name="Near-surface temperature (3-m)",fill=.TRUE.)
        call def_var_info(surf(17),trim(file_surface),"TTMIN","t3m_min", units="degrees Celcius", &
            long_name="Min near-surface temperature (3-m)",fill=.TRUE.)
        call def_var_info(surf(18),trim(file_surface),"TTMAX","t3m_max", units="degrees Celcius", &
            long_name="Max near-surface temperature (3-m)",fill=.TRUE.)
        call def_var_info(surf(14),trim(file_surface),"QQ",  "q",   units="g kg**-1", &
            long_name="Near-surface specific humidity",fill=.TRUE.)
        call def_var_info(surf(15),trim(file_surface),"SP",  "sp",  units="hPa", &
            long_name="Surface pressure",fill=.TRUE.)
        call def_var_info(surf(16),trim(file_surface),"RH",  "rh",  units="%", &
            long_name="Relative humidity",fill=.TRUE.)
        call def_var_info(surf(19),trim(file_surface),"UV",  "uv",  units="m s**-1", &
            long_name="Surface wind, magnitude",fill=.TRUE.)
        call def_var_info(surf(20),trim(file_surface),"SWD", "swd", units="W m**-2", &
            long_name="Near-surface radiation, shortwave down",fill=.TRUE.)
        call def_var_info(surf(21),trim(file_surface),"LWD", "lwd", units="W m**-2", &
            long_name="Near-surface radiation, longwave down",fill=.TRUE.)
        call def_var_info(surf(22),trim(file_surface),"LWU", "lwu", units="W m**-2", &
            long_name="Near-surface radiation, longwave up",fill=.TRUE.)
        call def_var_info(surf(23),trim(file_surface),"SHF", "shf", units="W m**-2", &
            long_name="Sensible heat flux",fill=.TRUE.)
        call def_var_info(surf(24),trim(file_surface),"LHF", "lhf", units="W m**-2", &
            long_name="Latent heat flux",fill=.TRUE.)
        call def_var_info(surf(25),trim(file_surface),"AL",  "al",  units="(0 - 1)", &
            long_name="Surface albedo")
        call def_var_info(surf(26),trim(file_surface),"CC",  "cc",  units="(0 - 1)", &
            long_name="Cloud cover fraction",fill=.TRUE.)
        call def_var_info(surf(27),trim(file_surface),"ST",  "ts",  units="degrees Celcius", &
            long_name="Surface temperature")
        call def_var_info(surf(28),trim(file_surface),"PDD", "pdd", units="degrees Celcius", &
            long_name="Positive degree days")
        
        nm       = 12
        n_var    = size(surf)

        if (present(max_neighbors) .and. present(lat_lim)) then 

            ! Allocate the input grid variable
            call grid_allocate(gMAR,invar)
            call grid_allocate(gMAR,invar_int)

            ! Initialize mapping
            call map_init(map,gMAR,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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
            
            ! Write meta data 
            call nc_write_attr(filename,"Description",desc)
            call nc_write_attr(filename,"Reference",ref)

            ! ## INVARIANT FIELDS ##
            do i = 1, size(invariant)

                var_now = invariant(i)
                call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value)
                outvar = missing_value
                call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,50.d3, &
                               fill=.TRUE.,missing_value=missing_value)
                if (var_now%fill) call fill_mean(outvar,missing_value=missing_value,fill_value=0.d0)
                if (trim(var_now%method) .eq. "nn") then 
                    call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc", &
                                  missing_value=nint(missing_value))
                else 
                    call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc", &
                                  missing_value=real(missing_value))
                end if 

                ! Write variable metadata
                call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            
            end do 

            n_prefix = 1 
            do k = 1, nyr 

                year = year0 + (k-1) 
                if (year .ge. year_switch) n_prefix = 2
                write(*,*) "=== ",year," ==="

                ! ## SURFACE FIELDS ##
                do i = 1, n_var

                    var_now = surf(i)     
                    write(var_now%filename,"(a,a,i4,a3)") &
                        trim(adjustl(file_surface)), trim(file_prefix(n_prefix)),year,".nc"
                    
                    do m = 1, nm 
                    
                        call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
                                 start=[1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1])
                        
                        ! Bug fix with input values - make sure missing values are missing
                        where (invar .lt. -9000.d0) invar = missing_value 

                        where (invar .ne. missing_value) invar = invar*var_now%conv 
                        outvar = missing_value 
                        call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",50.d3, &
                                       fill=.TRUE.,missing_value=missing_value)
                        if (var_now%fill) call fill_weighted(outvar,missing_value=missing_value)
                        call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                      start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])
                    
                        ! Write variable metadata
                        call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
                        call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            
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
            
            ! Write meta data 
            call nc_write_attr(filename_clim,"Description",desc)
            call nc_write_attr(filename_clim,"Reference",ref)

            ! ## INVARIANT FIELDS ##
            do i = 1, size(invariant)
                var_now = invariant(i) 
                call nc_read(filename,var_now%nm_out,var2D)
                call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc")
            
                ! Write variable metadata
                call nc_write_attr(filename_clim,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename_clim,var_now%nm_out,"long_name",var_now%long_name)
            
            end do 

            do i = 1, n_var
                var_now = surf(i)

                do m = 1, nm  
                    call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk])
                    var2D = time_average(var3D)
                    call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                  start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
                end do

                ! Write variable metadata
                call nc_write_attr(filename_clim,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename_clim,var_now%nm_out,"long_name",var_now%long_name)
             
            end do 

        end if 

        return 

    end subroutine MARv35_to_grid




! ## AJR: To do below: add long_name attribution to def_var_info calls, 
! ##      write meta data to output files 



!     subroutine MARv33_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,clim_range)
!         ! Convert the variables to the desired grid format and write to file
!         ! =========================================================
!         !
!         !       MAR (RCM) DATA - MARv3.3 downloaded from the ftp site:
!         !       ftp://ftp.climato.be/fettweis/MARv3.3/Greenland
!         !       Data is available on the Bamber et al. (2001) 5km grid
!         !       domain="Greenland-ERA": ERA-40 + ERA-Interim combined datasets
!         !       domain="Greenland-MIROC5-RCP85": MIROC5 histo+rcp85 combined datasets
!         !
!         ! =========================================================

!         implicit none 

!         character(len=*) :: domain, outfldr 
!         type(grid_class) :: grid 
!         integer, optional :: max_neighbors 
!         double precision, optional :: lat_lim 
!         integer, optional :: clim_range(2)

!         character(len=512) :: filename 
!         type(grid_class)   :: gMAR
!         character(len=256) :: file_invariant, file_surface, file_prefix(2)
!         type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
!         double precision, allocatable :: invar(:,:)
!         integer, allocatable :: invar_int(:,:) 
!         integer :: plev(9) 

!         type(map_class)  :: map 
!         type(var_defs) :: var_now 
!         double precision, allocatable :: outvar(:,:)
!         integer, allocatable          :: outmask(:,:), outvar_int(:,:)

!         integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 
!         integer :: yearf, k0, nk 
!         character(len=512) :: filename_clim 
!         double precision, allocatable :: var3D(:,:,:), var2D(:,:)

!         ! Define input grid
!         if (trim(domain) .eq. "Greenland-ERA") then 
            
!             ! Define MAR (Bamber et al. 2001) grid and input variable field
!             call grid_init(gMAR,name="Bamber01-5KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
!                            x0=-800.d0,dx=5.d0,nx=301,y0=-3400.d0,dy=5.d0,ny=561, &
!                            lambda=-39.d0,phi=90.d0,alpha=7.5d0)

!             ! Define the input filenames
!             file_invariant = "/data/sicopolis/data/MARv3.3/Greenland/ERA_1958-2013_15km/"// &
!                              "MARv3.3-15km-monthly-ERA-Interim-2013.nc"
!             file_surface   = "/data/sicopolis/data/MARv3.3/Greenland/"
!             file_prefix(1) = "ERA_1958-2013_15km/MARv3.3-15km-monthly-ERA-Interim-"
!             file_prefix(2) = "ERA_1958-2013_15km/MARv3.3-15km-monthly-ERA-Interim-"

!             ! Define the output filename 
!             write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
!                               "_MARv3.3-15km-monthly-ERA-Interim_195801-201312.nc"

!             year0       = 1958 
!             year_switch = 1979   ! Switch scenarios (ERA-40 to ERA-INTERIM)
!             nyr         = 2013-1958+1

!             ! For climatology
!             if (present(clim_range)) then  
!                 k0 = clim_range(1) - year0+1
!                 nk = clim_range(2) - clim_range(1) + 1 

!                 write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
!                     "_MARv3.3-15km-monthly-ERA-Interim_",clim_range(1),"-",clim_range(2),".nc"
!             end if 

!         else if (trim(domain) .eq. "Greenland-MIROC5-RCP85") then 

!             ! Define MAR (Bamber et al. 2001) grid and input variable field
!             call grid_init(gMAR,name="Bamber01-5KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
!                            x0=-800.d0,dx=5.d0,nx=301,y0=-3400.d0,dy=5.d0,ny=561, &
!                            lambda=-39.d0,phi=90.d0,alpha=7.5d0)

!             ! Define the input filenames
!             file_invariant = "/data/sicopolis/data/MARv3.3/Greenland/MIROC5-histo_1976-2005_30km/"// &
!                              "MARv3.3-monthly-MIROC5-histo-1976.nc"
!             file_surface   = "/data/sicopolis/data/MARv3.3/Greenland/"
!             file_prefix(1) = "MIROC5-histo_1976-2005_30km/MARv3.3-monthly-MIROC5-histo-"
!             file_prefix(2) = "MIROC5-rcp85_2006-2100_30km/MARv3.3-monthly-MIROC5-rcp85-"

!             ! Define the output filename 
!             write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
!                               "_MARv3.3-30km-monthly-MIROC5-rcp85_197601-210012.nc"

!             year0       = 1976
!             year_switch = 2006   ! Switch scenarios (historical to RCP85)
!             nyr         = 2100-1976+1

!             ! For climatology
!             if (present(clim_range)) then  
!                 k0 = clim_range(1) - year0+1
!                 nk = clim_range(2) - clim_range(1) + 1 

!                 write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
!                     "_MARv3.3-30km-monthly-MIROC5-rcp85_",clim_range(1),"-",clim_range(2),".nc"
!             end if 

!         else if (trim(domain) .eq. "Antarctica") then 

!             ! TODO 

!         else

!             write(*,*) "Domain not recognized: ",trim(domain)
!             stop 
!         end if 

!         ! Define the variables to be mapped 
!         allocate(invariant(2))
!         call def_var_info(invariant(1),trim(file_invariant),"MSK_MAR","mask",units="(0 - 2)",method="nn",fill=.TRUE.)
!         call def_var_info(invariant(2),trim(file_invariant),"SRF_MAR","zs",units="m",fill=.TRUE.)

!         allocate(surf(19))
!         call def_var_info(surf( 1),trim(file_surface),"SMB", "smb", units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
!         call def_var_info(surf( 2),trim(file_surface),"RU",  "ru",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
!         call def_var_info(surf( 3),trim(file_surface),"ME",  "me",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
!         call def_var_info(surf( 4),trim(file_surface),"ST",  "ts",  units="degrees Celcius")
!         call def_var_info(surf( 5),trim(file_surface),"TT",  "t3m", units="degrees Celcius",fill=.TRUE.)
!         call def_var_info(surf( 6),trim(file_surface),"SF",  "sf",  units="mm d**-1",conv=12.d0/365.d0,fill=.TRUE.)  ! mm/month => mm/day
!         call def_var_info(surf( 7),trim(file_surface),"RF",  "rf",  units="mm d**-1",conv=12.d0/365.d0,fill=.TRUE.)  ! mm/month => mm/day
!         call def_var_info(surf( 8),trim(file_surface),"SU",  "su",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
!         call def_var_info(surf( 9),trim(file_surface),"AL",  "al",  units="(0 - 1)")
!         call def_var_info(surf(10),trim(file_surface),"SWD", "swd", units="W m**-2",fill=.TRUE.)
!         call def_var_info(surf(11),trim(file_surface),"LWD", "lwd", units="W m**-2",fill=.TRUE.)
!         call def_var_info(surf(12),trim(file_surface),"SHF", "shf", units="W m**-2",fill=.TRUE.)
!         call def_var_info(surf(13),trim(file_surface),"LHF", "lhf", units="W m**-2",fill=.TRUE.)
!         call def_var_info(surf(14),trim(file_surface),"SP",  "sp",  units="hPa",fill=.TRUE.)
!         call def_var_info(surf(15),trim(file_surface),"UU",  "u",   units="m s**-1",fill=.TRUE.)
!         call def_var_info(surf(16),trim(file_surface),"VV",  "v",   units="m s**-1",fill=.TRUE.)
!         call def_var_info(surf(17),trim(file_surface),"QQ",  "q",   units="g kg**-1",fill=.TRUE.)
!         call def_var_info(surf(18),trim(file_surface),"CC",  "cc",  units="(0 - 1)",fill=.TRUE.)
!         call def_var_info(surf(19),trim(file_surface),"SH3", "SH3", units="mm d**-1",conv=1d3*12.d0/365.d0)   ! m/month => mm/day

!         nm       = 12
!         n_var    = size(surf)
!         if (trim(domain) .ne. "Greenland-ERA") n_var = 16   ! Exclude QQ, CC and SH3 if not available
    
!         if (present(max_neighbors) .and. present(lat_lim)) then 

!             ! Allocate the input grid variable
!             call grid_allocate(gMAR,invar)
!             call grid_allocate(gMAR,invar_int)

!             ! Initialize mapping
!             call map_init(map,gMAR,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

!             ! Initialize output variable arrays
!             call grid_allocate(grid,outvar)
!             call grid_allocate(grid,outmask)    
!             call grid_allocate(grid,outvar_int)

!             ! Initialize the output file
!             call nc_create(filename)
!             call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
!             call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
!             call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
!             call nc_write_dim(filename,"time", x=year0,dx=1,nx=nyr,units="years",calendar="360_day")
!             call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
!             ! ## INVARIANT FIELDS ##
!             do i = 1, size(invariant)

!                 var_now = invariant(i)
!                 call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value)
!                 outvar = missing_value
!                 call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,50.d3, &
!                                fill=.TRUE.,missing_value=missing_value)
!                 if (var_now%fill) call fill_mean(outvar,missing_value=missing_value,fill_value=0.d0)
!                 if (trim(var_now%method) .eq. "nn") then 
!                     call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc", &
!                                   units=var_now%units_out,missing_value=nint(missing_value))
!                 else 
!                     call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc", &
!                                   units=var_now%units_out,missing_value=real(missing_value))
!                 end if 

!             end do 

!             n_prefix = 1 
!             do k = 1, nyr 

!                 year = year0 + (k-1) 
!                 if (year .ge. year_switch) n_prefix = 2
!                 write(*,*) "=== ",year," ==="
         
!                 do m = 1, nm 
!                     q = m 
!                     write(*,*) "month ",m

!                     ! ## SURFACE FIELDS ##
!                     do i = 1, n_var
!                         var_now = surf(i)     
!                         write(var_now%filename,"(a,a,i4,a3)") &
!                             trim(adjustl(file_surface)), trim(file_prefix(n_prefix)),year,".nc"
!                         call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
!                                  start=[1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1])
!                         where (invar .ne. missing_value) invar = invar*var_now%conv 
!                         outvar = missing_value 
!                         call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",50.d3, &
!                                        fill=.TRUE.,missing_value=missing_value)
!                         if (var_now%fill) call fill_weighted(outvar,missing_value=missing_value)
!                         call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month",dim4="time", &
!                                       units=var_now%units_out,start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])
!                     end do 

!                 end do 
!             end do 
        
!         end if 

!         if (present(clim_range)) then 

!             ! Create climatology too (month by month)

!             call grid_allocate(grid,var2D)
!             allocate(var3D(grid%G%nx,grid%G%ny,nk))    
            
!             ! Initialize the output file
!             call nc_create(filename_clim)
!             call nc_write_dim(filename_clim,"xc",   x=grid%G%x,units="kilometers")
!             call nc_write_dim(filename_clim,"yc",   x=grid%G%y,units="kilometers")
!             call nc_write_dim(filename_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
!             call grid_write(grid,filename_clim,xnm="xc",ynm="yc",create=.FALSE.)
        
!             ! ## INVARIANT FIELDS ##
!             do i = 1, size(invariant)
!                 var_now = invariant(i) 
!                 call nc_read(filename,var_now%nm_out,var2D)
!                 call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc", &
!                               units=var_now%units_out)
!             end do 

!             do i = 1, n_var
!                 var_now = surf(i)

!                 do m = 1, nm  
!                     call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk])
!                     var2D = time_average(var3D)
!                     call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
!                                   units=var_now%units_out,start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
!                 end do 
!             end do 

!         end if 

!         return 

!     end subroutine MARv33_to_grid

!     subroutine MARv32_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
!         ! Convert the variables to the desired grid format and write to file
!         ! =========================================================
!         !
!         !       MAR (RCM) DATA - MARv3.2 original data passed by Xavier
!         !
!         ! =========================================================

!         implicit none 

!         character(len=*) :: domain, outfldr 
!         type(grid_class) :: grid 
!         integer :: max_neighbors 
!         double precision :: lat_lim 
!         character(len=512) :: filename 

!         type(grid_class)   :: gMAR
!         character(len=256) :: file_invariant, file_surface, file_prefix(2)
!         type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
!         double precision, allocatable :: invar(:,:) 
!         integer :: plev(9) 

!         type(map_class)  :: map 
!         type(var_defs) :: var_now 
!         double precision, allocatable :: outvar(:,:)
!         integer, allocatable          :: outmask(:,:)

!         integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 

!         ! Define ECMWF input grid
!         if (trim(domain) .eq. "Greenland-ERA") then 
            
!             ! Define MAR grid and input variable field
!             call grid_init(gMAR,name="MAR-25KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
!                            x0=-750.d0,dx=25.d0,nx=58,y0=-1200.d0,dy=25.d0,ny=108, &
!                            lambda=-40.d0,phi=71.d0,alpha=7.5d0)

!             ! Define the input filenames
!             file_invariant = "data/MAR/MAR_ERA-INTERIM/MARv3.2_historical_mon_197901-197912.nc"
!             file_surface   = "data/MAR/"
!             file_prefix(1) = "MAR_ERA-INTERIM/MARv3.2_historical_mon_"

!             ! Define the output filename 
!             write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
!                               "_MARv3.2-ERA-INTERIM_197901-201112.nc"

!             year0       = 1979
!             year_switch = 0   ! Switch scenarios (ERA-40 to ERA-INTERIM)
!             nyr         = 2011-1979+1

!         else

!             write(*,*) "Domain not recognized: ",trim(domain)
!             stop 
!         end if 

!         ! Define the variables to be mapped 
!         allocate(invariant(3))
!         call def_var_info(invariant(1),trim(file_invariant),"SRF","mask_srf",units="(1 - 4)",method="nn")
!         call def_var_info(invariant(2),trim(file_invariant),"SOL","mask_sol",units="(0 - 12)",method="nn")
!         call def_var_info(invariant(3),trim(file_invariant),"SH","zs",  units="m")

!         allocate(surf(23))
!         call def_var_info(surf( 1),trim(file_surface),"SMB", "smb", units="mm d**-1",dimextra=.TRUE.)
!         call def_var_info(surf( 2),trim(file_surface),"RU",  "ru",  units="mm d**-1")
!         call def_var_info(surf( 3),trim(file_surface),"ME",  "me",  units="mm d**-1",dimextra=.TRUE.)
!         call def_var_info(surf( 4),trim(file_surface),"RZ",  "rz",  units="mm d**-1",dimextra=.TRUE.)
!         call def_var_info(surf( 5),trim(file_surface),"SF",  "sf",  units="mm d**-1")
!         call def_var_info(surf( 6),trim(file_surface),"RF",  "rf",  units="mm d**-1")
!         call def_var_info(surf( 7),trim(file_surface),"SU",  "su",  units="mm d**-1",dimextra=.TRUE.)
!         call def_var_info(surf( 8),trim(file_surface),"SF",  "sf",  units="mm d**-1")
!         call def_var_info(surf( 9),trim(file_surface),"TT",  "t3m", units="degrees Celcius",dimextra=.TRUE.)
!         call def_var_info(surf(10),trim(file_surface),"QQ",  "Q",   units="g kg**-1",dimextra=.TRUE.)
!         call def_var_info(surf(11),trim(file_surface),"UU",  "u",   units="m s**-1",dimextra=.TRUE.)
!         call def_var_info(surf(12),trim(file_surface),"VV",  "v",   units="m s**-1",dimextra=.TRUE.)
!         call def_var_info(surf(13),trim(file_surface),"SP",  "sp",  units="hPa")
!         call def_var_info(surf(14),trim(file_surface),"SWD", "swd", units="W m**-2")
!         call def_var_info(surf(15),trim(file_surface),"LWD", "lwd", units="W m**-2")
!         call def_var_info(surf(16),trim(file_surface),"LWU", "lwu", units="W m**-2")
!         call def_var_info(surf(17),trim(file_surface),"SHF", "shf", units="W m**-2")
!         call def_var_info(surf(18),trim(file_surface),"LHF", "lhf", units="W m**-2")
!         call def_var_info(surf(19),trim(file_surface),"AL1", "al1", units="(0 - 1)")
!         call def_var_info(surf(20),trim(file_surface),"AL2", "al2", units="(0 - 1)")
!         call def_var_info(surf(21),trim(file_surface),"CC",  "cc",  units="(0 - 1)")
!         call def_var_info(surf(22),trim(file_surface),"STT", "ts",  units="degrees Celcius",dimextra=.TRUE.)
!         call def_var_info(surf(23),trim(file_surface),"SHSN2","Hs", units="m",dimextra=.TRUE.)
    
!         ! Allocate the input grid variable
!         call grid_allocate(gMAR,invar)

!         ! Initialize mapping
!         call map_init(map,gMAR,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

!         ! Initialize output variable arrays
!         call grid_allocate(grid,outvar)
!         call grid_allocate(grid,outmask)    
        
!         ! Initialize the output file
!         call nc_create(filename)
!         call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
!         call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
!         call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
!         call nc_write_dim(filename,"time", x=year0,dx=1,nx=nyr,units="years",calendar="360_day")
!         call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
    
!         ! ## INVARIANT FIELDS ##
!         do i = 1, size(invariant)
!             var_now = invariant(i) 
!             call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value)
!             outvar = missing_value 
!             call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,100.d3, &
!                            fill=.FALSE.,missing_value=missing_value)
!             if (var_now%method .eq. "nn") then 
!                 call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
!             else
!                 call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
!             end if  
!         end do 

!         nm       = 12    
!         n_prefix = 1 
!         n_var    = size(surf)

!         do k = 1, nyr 

!             year = year0 + (k-1) 
!             if (year .ge. year_switch) n_prefix = 1
!             write(*,*) "=== ",year," ==="
     
!             do m = 1, nm 
!                 q = m 
!                 write(*,*) "month ",m

!                 ! ## SURFACE FIELDS ##
!                 do i = 1, n_var
!                     var_now = surf(i)     
!                     write(var_now%filename,"(a,a,i4,a3,i4,a5)")  &
!                         trim(file_surface),trim(file_prefix(n_prefix)),year,"01-",year,"12.nc"
!                     if (var_now%dimextra) then 
!                         call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
!                                       start=[1,1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1,1])
!                     else 
!                         call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
!                                  start=[1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1])
!                     end if
!                     where (invar .ne. missing_value) invar = invar*var_now%conv 
!                     outvar = missing_value 
!                     call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",100.d3, &
!                                    fill=.FALSE.,missing_value=missing_value)
!                     call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month",dim4="time", &
!                                   units=var_now%units_out,start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])
!                 end do 

!             end do 
!         end do 
        
!         return 

!     end subroutine MARv32_to_grid


end module MAR 
