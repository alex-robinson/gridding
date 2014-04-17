
module gridding_datasets

    use coordinates 
    use ncio 
    use vargrid 

    implicit none 

    private
    public :: ecmwf_to_grid, MARv33_to_grid, MARv32_to_grid

contains

    subroutine ecmwf_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file

        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 

        type(grid_class)   :: gECMWF 
        character(len=256) :: file_invariant, file_surface, files_pres(9)
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:) 
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: nyr, nm, q, k, year, m, i, l 

        ! Define ECMWF input grid
        if (trim(domain) .eq. "GRL075") then 
            
            ! Initialize the grid
            call grid_init(gECMWF,name="ECMWF-GRL075",mtype="latlon",units="kilometers",lon180=.TRUE., &
                           x0=-100.d0,dx=0.75d0,nx=161,y0=49.5d0,dy=0.75d0,ny=55)
            
            ! Assign the filenames
            file_invariant = "data/ECMWF/NEW/ERA-INTERIM-GRL-invariant_historical_mon_197901-201212.nc"
            file_surface   = "data/ECMWF/NEW/ERA-INTERIM-GRL-surface_historical_mon_197901-201212.nc"
            files_pres(1)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-1000Mb_historical_mon_197901-201212.nc"
            files_pres(2)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-950Mb_historical_mon_197901-201212.nc"
            files_pres(3)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-850Mb_historical_mon_197901-201212.nc"
            files_pres(4)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-750Mb_historical_mon_197901-201212.nc"
            files_pres(5)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-700Mb_historical_mon_197901-201212.nc"
            files_pres(6)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-650Mb_historical_mon_197901-201212.nc"
            files_pres(7)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-600Mb_historical_mon_197901-201212.nc"
            files_pres(8)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-550Mb_historical_mon_197901-201212.nc"
            files_pres(9)  = "data/ECMWF/NEW/ERA-INTERIM-GRL-500Mb_historical_mon_197901-201212.nc"

        else if (trim(domain) .eq. "ANT075") then 

            ! TODO 

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the pressure levels to be mapped
        plev = [1000,950,850,750,700,650,600,550,500]

        ! Define the variables to be mapped 
        
        allocate(invariant(1))
        call def_var_info(invariant(1),trim(file_invariant),"z","zs",units="m")

        allocate(surf(12))
        call def_var_info(surf( 1),trim(file_surface),"sp", "sp", units="Pa")
        call def_var_info(surf( 2),trim(file_surface),"tcw","tcw",units="kg m**-2")
        call def_var_info(surf( 3),trim(file_surface),"tclw","tclw",units="kg m**-2")
        call def_var_info(surf( 4),trim(file_surface),"tciw","tciw",units="kg m**-2")
        call def_var_info(surf( 5),trim(file_surface),"p56.162","clw",units="kg m**-2")
        call def_var_info(surf( 6),trim(file_surface),"p57.162","ciw",units="kg m**-2")
        call def_var_info(surf( 7),trim(file_surface),"tcc","tcc",units="(0 - 1)")
        call def_var_info(surf( 8),trim(file_surface),"u10","u10",units="m s**-1")
        call def_var_info(surf( 9),trim(file_surface),"v10","v10",units="m s**-1")
        call def_var_info(surf(10),trim(file_surface),"t2m","t2m",units="K")
        call def_var_info(surf(11),trim(file_surface),"al", "al", units="(0 - 1)")
        call def_var_info(surf(12),trim(file_surface),"sst","sst",units="K")

        allocate(pres(5))
        call def_var_info(pres( 1),"None","t", "p_t",units="K",         plev="plev",filenames=files_pres)
        call def_var_info(pres( 2),"None","z", "p_z",units="m**2 s**-2",plev="plev",filenames=files_pres)
        call def_var_info(pres( 3),"None","u", "p_u",units="m s**-1",   plev="plev",filenames=files_pres)
        call def_var_info(pres( 4),"None","v", "p_v",units="m s**-1",   plev="plev",filenames=files_pres)
        call def_var_info(pres( 5),"None","w", "p_w",units="Pa s**-1",  plev="plev",filenames=files_pres)

        ! Allocate the input grid variable
        call grid_allocate(gECMWF,invar)

        ! Initialize mapping
        call map_init(map,gECMWF,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    

        ! ## First make file for surface fields including invariants ##
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_ERA-INTERIM"// &
                             "_197901-201212.nc"

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"plev", x=dble(plev),units="hPa")
        call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call nc_write_dim(filename,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
    
        ! ## INVARIANT FIELDS ##
        var_now = invariant(1) 
        call nc_read(var_now%filename,var_now%nm_in,invar)
        call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",400.d3,missing_value=missing_value)
        call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",units=var_now%units_out)

        nyr = 2012-1979+1
        nm  = 12 

        q = 0 
        do k = 1, nyr 

            year = 1978 + k 
            write(*,*) "=== ",year," ==="

            do m = 1, nm 
                q = q+1 

                write(*,*) "month ",m

                ! ## SURFACE FIELDS ##
                do i = 1, size(surf)
                    var_now = surf(i) 
                    call nc_read(var_now%filename,var_now%nm_in,invar,start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1])
                    call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",400.d3,missing_value=missing_value)
                    call nc_write(filename,var_now%nm_out,real(outvar),  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                  units=var_now%units_out,start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])
                end do 
            end do
        end do  

        do l = 1, size(files_pres)   ! Loop over pressure layers

            ! ## Make one file for each pressure level ##
            if (plev(l) .ge. 1000) then 
                write(filename,"(a,i4,a)") trim(outfldr)//"/"//trim(grid%name)//"_ERA-INTERIM-", &
                                           plev(l),"Mb_197901-201212.nc"
            else
                write(filename,"(a,i3,a)") trim(outfldr)//"/"//trim(grid%name)//"_ERA-INTERIM-", &
                                           plev(l),"Mb_197901-201212.nc"
            end if 

            ! Initialize the output file
            call nc_create(filename)
            call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
            call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
            call nc_write_dim(filename,"plev", x=dble(plev),units="hPa")
            call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
            call nc_write_dim(filename,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
            call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
            q = 0 
            do k = 1, nyr 

                year = 1978 + k 
                write(*,*) "=== ",year," ==="

                do m = 1, nm 
                    q = q+1 

                    write(*,*) "month ",m

                    ! ## PRESSURE FIELDS ##
                    do i = 1, size(pres)
                        var_now = pres(i) 

                        call nc_read(var_now%filenames(l),var_now%nm_in,invar,start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1])
                        call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",400.d3,missing_value=missing_value)
                        call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                      units=var_now%units_out,start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])

                    end do 
                end do 
            end do 

        end do 

        return 

    end subroutine ecmwf_to_grid

    subroutine MARv33_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file

        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 

        type(grid_class)   :: gMAR
        character(len=256) :: file_invariant, file_surface, file_prefix(2)
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:) 
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 

        ! Define ECMWF input grid
        if (trim(domain) .eq. "Greenland-ERA") then 
            
            ! Define MAR (Bamber et al. 2001) grid and input variable field
            call grid_init(gMAR,name="Bamber01-5KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                           x0=-800.d0,dx=5.d0,nx=301,y0=-3400.d0,dy=5.d0,ny=561, &
                           lambda=-39.d0,phi=90.d0,alpha=7.5d0)

            ! Define the input filenames
            file_invariant = "/data/sicopolis/data/MARv3.3/Greenland/ERA_1958-2013_15km/"// &
                             "MARv3.3-15km-monthly-ERA-Interim-2013.nc"
            file_surface   = "/data/sicopolis/data/MARv3.3/Greenland/"
            file_prefix(1) = "ERA_1958-2013_15km/MARv3.3-15km-monthly-"
            file_prefix(2) = "ERA_1958-2013_15km/MARv3.3-15km-monthly-"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_MARv3.3-15km-monthly-ERA-Interim_195801-201312.nc"

            year0       = 1958 
            year_switch = 1979   ! Switch scenarios (ERA-40 to ERA-INTERIM)
            nyr         = 2013-1958+1

        else if (trim(domain) .eq. "Greenland-MIROC5-RCP85") then 

            ! Define MAR (Bamber et al. 2001) grid and input variable field
            call grid_init(gMAR,name="Bamber01-5KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                           x0=-800.d0,dx=5.d0,nx=301,y0=-3400.d0,dy=5.d0,ny=561, &
                           lambda=-39.d0,phi=90.d0,alpha=7.5d0)

            ! Define the input filenames
            file_invariant = "/data/sicopolis/data/MARv3.3/Greenland/MIROC5-histo_1976-2005_30km/"// &
                             "MARv3.3-monthly-MIROC5-histo-1976.nc"
            file_surface   = "/data/sicopolis/data/MARv3.3/Greenland/"
            file_prefix(1) = "MIROC5-histo_1976-2005_30km/MARv3.3-monthly-MIROC5-histo-"
            file_prefix(2) = "MIROC5-rcp85_2006-2100_30km/MARv3.3-monthly-MIROC5-rcp85-"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_MARv3.3-30km-monthly-MIROC5-rcp85_197601-210012.nc"

            year0       = 1976
            year_switch = 2006   ! Switch scenarios (historical to RCP85)
            nyr         = 2100-1976+1

        else if (trim(domain) .eq. "Antarctica") then 

            ! TODO 

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(2))
        call def_var_info(invariant(1),trim(file_invariant),"MSK_MAR","mask",units="(0 - 2)",method="nn")
        call def_var_info(invariant(2),trim(file_invariant),"SRF_MAR","zs",units="m")

        allocate(surf(17))
        call def_var_info(surf( 1),trim(file_surface),"SMB", "smb", units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
        call def_var_info(surf( 2),trim(file_surface),"RU",  "ru",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
        call def_var_info(surf( 3),trim(file_surface),"ME",  "me",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
        call def_var_info(surf( 4),trim(file_surface),"ST",  "ts",  units="degrees Celcius")
        call def_var_info(surf( 5),trim(file_surface),"TT",  "t3m", units="degrees Celcius")
        call def_var_info(surf( 6),trim(file_surface),"SF",  "sf",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
        call def_var_info(surf( 7),trim(file_surface),"RF",  "rf",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
        call def_var_info(surf( 8),trim(file_surface),"SU",  "su",  units="mm d**-1",conv=12.d0/365.d0)  ! mm/month => mm/day
        call def_var_info(surf( 9),trim(file_surface),"AL",  "al",  units="(0 - 1)")
        call def_var_info(surf(10),trim(file_surface),"SWD", "swd", units="W m**-2")
        call def_var_info(surf(11),trim(file_surface),"LWD", "lwd", units="W m**-2")
        call def_var_info(surf(12),trim(file_surface),"SHF", "shf", units="W m**-2")
        call def_var_info(surf(13),trim(file_surface),"LHF", "lhf", units="W m**-2")
        call def_var_info(surf(14),trim(file_surface),"SP",  "sp",  units="hPa")
        call def_var_info(surf(15),trim(file_surface),"UU",  "u",   units="m s**-1")
        call def_var_info(surf(16),trim(file_surface),"VV",  "v",   units="m s**-1")
        call def_var_info(surf(17),trim(file_surface),"SH3", "SH3", units="mm d**-1",conv=1d3*12.d0/365.d0)   ! m/month => mm/day

        ! Allocate the input grid variable
        call grid_allocate(gMAR,invar)

        ! Initialize mapping
        call map_init(map,gMAR,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call nc_write_dim(filename,"time", x=year0,dx=1,nx=nyr,units="years",calendar="360_day")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
    
        ! ## INVARIANT FIELDS ##
        do i = 1, size(invariant)
            var_now = invariant(i) 
            call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value)
            outvar = missing_value 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,100.d3, &
                           fill=.FALSE.,missing_value=missing_value)
            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            end if  
        end do 

        nm       = 12    
        n_prefix = 1 
        n_var    = size(surf)
        if (trim(domain) .ne. "Greenland-ERA") n_var = size(surf)-1   ! Exclude SH3 if it doesn't exist

        do k = 1, nyr 

            year = year0 + (k-1) 
            if (year .ge. year_switch) n_prefix = 2
            write(*,*) "=== ",year," ==="
     
            do m = 1, nm 
                q = m 
                write(*,*) "month ",m

                ! ## SURFACE FIELDS ##
                do i = 1, n_var
                    var_now = surf(i)     
                    write(var_now%filename,"(a,a,i4,a3)") trim(file_surface),trim(file_prefix(n_prefix)),year,".nc"
                    call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value, &
                             start=[1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1])
                    where (invar .ne. missing_value) invar = invar*var_now%conv 
                    outvar = missing_value 
                    call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",100.d3, &
                                   fill=.FALSE.,missing_value=missing_value)
                    call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                  units=var_now%units_out,start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])
                end do 

            end do 
        end do 
        
        return 

    end subroutine MARv33_to_grid


    subroutine MARv32_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file

        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 

        type(grid_class)   :: gMAR
        character(len=256) :: file_invariant, file_surface, file_prefix(2)
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:) 
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 

        ! Define ECMWF input grid
        if (trim(domain) .eq. "Greenland-ERA") then 
            
            ! Define MAR grid and input variable field
            call grid_init(gMAR,name="MAR-25KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                           x0=-750.d0,dx=25.d0,nx=58,y0=-1200.d0,dy=25.d0,ny=108, &
                           lambda=-40.d0,phi=71.d0,alpha=7.5d0)

            ! Define the input filenames
            file_invariant = "data/MAR/MAR_ERA-INTERIM/MARv3.2_historical_mon_197901-197912.nc"
            file_surface   = "data/MAR/"
            file_prefix(1) = "MAR_ERA-INTERIM/MARv3.2_historical_mon_"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_MARv3.2-ERA-INTERIM_197901-201112.nc"

            year0       = 1979
            year_switch = 0   ! Switch scenarios (ERA-40 to ERA-INTERIM)
            nyr         = 2011-1979+1

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(3))
        call def_var_info(invariant(1),trim(file_invariant),"SRF","mask_srf",units="(1 - 4)",method="nn")
        call def_var_info(invariant(2),trim(file_invariant),"SOL","mask_sol",units="(0 - 12)",method="nn")
        call def_var_info(invariant(3),trim(file_invariant),"SH","zs",  units="m")

        allocate(surf(23))
        call def_var_info(surf( 1),trim(file_surface),"SMB", "smb", units="mm d**-1",dimextra=.TRUE.)
        call def_var_info(surf( 2),trim(file_surface),"RU",  "ru",  units="mm d**-1")
        call def_var_info(surf( 3),trim(file_surface),"ME",  "me",  units="mm d**-1",dimextra=.TRUE.)
        call def_var_info(surf( 4),trim(file_surface),"RZ",  "rz",  units="mm d**-1",dimextra=.TRUE.)
        call def_var_info(surf( 5),trim(file_surface),"SF",  "sf",  units="mm d**-1")
        call def_var_info(surf( 6),trim(file_surface),"RF",  "rf",  units="mm d**-1")
        call def_var_info(surf( 7),trim(file_surface),"SU",  "su",  units="mm d**-1",dimextra=.TRUE.)
        call def_var_info(surf( 8),trim(file_surface),"SF",  "sf",  units="mm d**-1")
        call def_var_info(surf( 9),trim(file_surface),"TT",  "t3m", units="degrees Celcius",dimextra=.TRUE.)
        call def_var_info(surf(10),trim(file_surface),"QQ",  "Q",   units="g kg**-1",dimextra=.TRUE.)
        call def_var_info(surf(11),trim(file_surface),"UU",  "u",   units="m s**-1",dimextra=.TRUE.)
        call def_var_info(surf(12),trim(file_surface),"VV",  "v",   units="m s**-1",dimextra=.TRUE.)
        call def_var_info(surf(13),trim(file_surface),"SP",  "sp",  units="hPa")
        call def_var_info(surf(14),trim(file_surface),"SWD", "swd", units="W m**-2")
        call def_var_info(surf(15),trim(file_surface),"LWD", "lwd", units="W m**-2")
        call def_var_info(surf(16),trim(file_surface),"LWU", "lwu", units="W m**-2")
        call def_var_info(surf(17),trim(file_surface),"SHF", "shf", units="W m**-2")
        call def_var_info(surf(18),trim(file_surface),"LHF", "lhf", units="W m**-2")
        call def_var_info(surf(19),trim(file_surface),"AL1", "al1", units="(0 - 1)")
        call def_var_info(surf(20),trim(file_surface),"AL2", "al2", units="(0 - 1)")
        call def_var_info(surf(21),trim(file_surface),"CC",  "cc",  units="(0 - 1)")
        call def_var_info(surf(22),trim(file_surface),"STT", "ts",  units="degrees Celcius",dimextra=.TRUE.)
        call def_var_info(surf(23),trim(file_surface),"SHSN2","Hs", units="m",dimextra=.TRUE.)
    
        ! Allocate the input grid variable
        call grid_allocate(gMAR,invar)

        ! Initialize mapping
        call map_init(map,gMAR,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call nc_write_dim(filename,"time", x=year0,dx=1,nx=nyr,units="years",calendar="360_day")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
    
        ! ## INVARIANT FIELDS ##
        do i = 1, size(invariant)
            var_now = invariant(i) 
            call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value)
            outvar = missing_value 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,100.d3, &
                           fill=.FALSE.,missing_value=missing_value)
            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            end if  
        end do 

        nm       = 12    
        n_prefix = 1 
        n_var    = size(surf)

        do k = 1, nyr 

            year = year0 + (k-1) 
            if (year .ge. year_switch) n_prefix = 1
            write(*,*) "=== ",year," ==="
     
            do m = 1, nm 
                q = m 
                write(*,*) "month ",m

                ! ## SURFACE FIELDS ##
                do i = 1, n_var
                    var_now = surf(i)     
                    write(var_now%filename,"(a,a,i4,a3,i4,a5)")  &
                        trim(file_surface),trim(file_prefix(n_prefix)),year,"01-",year,"12.nc"
                    if (var_now%dimextra) then 
                        call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value, &
                                      start=[1,1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1,1])
                    else 
                        call nc_read(var_now%filename,var_now%nm_in,invar,missing_value=missing_value, &
                                 start=[1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1])
                    end if
                    where (invar .ne. missing_value) invar = invar*var_now%conv 
                    outvar = missing_value 
                    call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",100.d3, &
                                   fill=.FALSE.,missing_value=missing_value)
                    call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                  units=var_now%units_out,start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])
                end do 

            end do 
        end do 
        
        return 

    end subroutine MARv32_to_grid


end module gridding_datasets
