
module gridding_datasets

    use coordinates 
    use interp_time
    use ncio 
    
    implicit none 

    double precision, parameter :: missing_value = -9999.d0
    
    type var_defs
        character(len=256) :: filename, filenames(20)
        character(len=256) :: nm_in, nm_out  
        character(len=256) :: units_in, units_out 
        character(len=256) :: method
        logical :: mask, dimextra
        character(len=256) :: plev
        double precision   :: conv 
        logical            :: fill 
    end type 

    private
    public :: Bamber13_to_grid, ecmwf_to_grid
    public :: MARv33_to_grid, MARv32_to_grid
    public :: CERES_to_grid 

contains

    subroutine Bamber13_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       TOPO DATA
        !
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 

        type(grid_class)   :: gTOPO
        character(len=256) :: file_invariant, file_surface, file_prefix(2)
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:) 
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 

        ! Define ECMWF input grid
        if (trim(domain) .eq. "Greenland") then 
            
            ! Define topography (Bamber et al. 2013) grid and input variable field
            call grid_init(gTOPO,name="TOPO-10KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-1300.d0,dx=10.d0,nx=251,y0=-3500.d0,dy=10.d0,ny=301, &
                   lambda=-39.d0,phi=90.d0,alpha=7.5d0)

            ! Define the input filenames
            file_invariant = "data/Greenland/Greenland_bedrock_topography_V3.nc"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(4))
        call def_var_info(invariant(1),trim(file_invariant),"BedrockElevation","zb",units="m")
        call def_var_info(invariant(2),trim(file_invariant),"SurfaceElevation","zs",units="m")
        call def_var_info(invariant(3),trim(file_invariant),"IceThickness",    "H", units="m")
        call def_var_info(invariant(4),trim(file_invariant),"LandMask",      "mask",units="(0 - 4",method="nn")

        ! Allocate the input grid variable
        call grid_allocate(gTOPO,invar)
        
        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp(2501,3001))

        ! Initialize mapping
        call map_init(map,gTOPO,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
    
        ! ## INVARIANT FIELDS ##
        do i = 1, size(invariant)
            var_now = invariant(i) 
            call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=missing_value)
            call thin(invar,tmp,by=10)
            if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
                where( invar .eq. missing_value ) invar = 0.d0 
            end if
            if (trim(var_now%nm_out) .eq. "zb") then 
                call fill(invar,missing_value=missing_value,fill_value=-1000.d0)
            end if 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,20.d3, &
                          fill=.TRUE.,missing_value=missing_value)
            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            end if 
        end do 

        return 

    end subroutine Bamber13_to_grid

    subroutine ecmwf_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,clim_range)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       ECMWF DATA (ERA-INTERIM 1979-2013)
        !
        ! ========================================================= 

        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer, optional :: max_neighbors 
        double precision, optional :: lat_lim 
        character(len=512) :: filename 
        integer, optional  :: clim_range(2) 

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
        integer :: yearf, k0, nk 
        character(len=512) :: filename_clim 
        double precision, allocatable :: var3D(:,:,:), var2D(:,:)

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

            ! ## First make file for surface fields including invariants ##
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_ERA-INTERIM_197901-201212.nc"

            ! For climatology
            if (present(clim_range)) then  
                k0 = clim_range(1) - 1979 + 1
                nk = clim_range(2) - clim_range(1) + 1 

                write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
                    "_ERA-INTERIM_",clim_range(1),"-",clim_range(2),".nc"
            end if 

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

        if (present(max_neighbors) .and. present(lat_lim)) then 

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
        
            ! ## INVARIANT FIELDS ##
            var_now = invariant(1) 
            call nc_read(trim(var_now%filename),var_now%nm_in,invar)
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
                        call nc_read(trim(var_now%filename),var_now%nm_in,invar,start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1])
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

                            call nc_read(trim(var_now%filenames(l)),var_now%nm_in,invar, &
                                         start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1])
                            call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",400.d3, &
                                           missing_value=missing_value)
                            call nc_write(filename,var_now%nm_out,real(outvar), &
                                          dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                          units=var_now%units_out,start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])

                        end do 
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
        
            ! ## INVARIANT FIELDS ##
            var_now = invariant(1) 
            call nc_read(filename,var_now%nm_out,var2D)
            call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc", &
                          units=var_now%units_out)

            do i = 1, size(surf)
                var_now = surf(i)
                write(*,*) "Processing ",trim(var_now%nm_out)
                do m = 1, nm  
                    call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk])
                    var2D = time_average(var3D)
                    call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                  units=var_now%units_out,start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
                end do 
            end do 

            do l = 1, size(files_pres)   ! Loop over pressure layers

                ! ## Make one file for each pressure level ##
                if (plev(l) .ge. 1000) then 
                    write(filename,"(a,i4,a)") trim(outfldr)//"/"//trim(grid%name)//"_ERA-INTERIM-", &
                                               plev(l),"Mb_197901-201212.nc"
                    write(filename_clim,"(a,i4,a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)//"_ERA-INTERIM-", &
                                               plev(l),"Mb_",clim_range(1),"-",clim_range(2),".nc"
                else
                    write(filename,"(a,i3,a)") trim(outfldr)//"/"//trim(grid%name)//"_ERA-INTERIM-", &
                                               plev(l),"Mb_197901-201212.nc"
                    write(filename_clim,"(a,i3,a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)//"_ERA-INTERIM-", &
                                               plev(l),"Mb_",clim_range(1),"-",clim_range(2),".nc"
                end if 

                ! Initialize the output file
                call nc_create(filename_clim)
                call nc_write_dim(filename_clim,"xc",   x=grid%G%x,units="kilometers")
                call nc_write_dim(filename_clim,"yc",   x=grid%G%y,units="kilometers")
                call nc_write_dim(filename_clim,"plev", x=dble(plev),units="hPa")
                call nc_write_dim(filename_clim,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
                call nc_write_dim(filename_clim,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
                call grid_write(grid,filename_clim,xnm="xc",ynm="yc",create=.FALSE.)
            
                do i = 1, size(pres)
                    var_now = pres(i)

                    do m = 1, nm  
                        call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk])
                        var2D = time_average(var3D)
                        call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                      units=var_now%units_out,start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
                    end do 
                end do 

            end do 

        end if 

        return 

    end subroutine ecmwf_to_grid

    subroutine CERES_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CERES DATA
        !
        ! =========================================================
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 

        type(grid_class)   :: gCERES
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
        if (trim(domain) .eq. "Global") then 
            
            ! Define MAR grid and input variable field
            call grid_init(gCERES,name="CERES-1deg",mtype="latlon",units="degrees",lon180=.FALSE., &
                   x0=0.5d0,dx=1d0,nx=360,y0=-90.d0,dy=1d0,ny=180 )

            ! Define the input filenames
            file_surface = "data/CERES/CERES_EBAF-TOA_Ed2.8_Subset_CLIM01-CLIM12.nc"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_CERES_2001-2013.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(surf(10))
        call def_var_info(surf( 1),trim(file_surface),"toa_sw_all_clim","toa_sw_all",units="W m**-2",method="radius")
        call def_var_info(surf( 2),trim(file_surface),"toa_sw_clr_clim","toa_sw_clr",units="W m**-2",method="radius")
        call def_var_info(surf( 3),trim(file_surface),"toa_lw_all_clim","toa_lw_all",units="W m**-2",method="radius")
        call def_var_info(surf( 4),trim(file_surface),"toa_lw_clr_clim","toa_lw_clr",units="W m**-2",method="radius")
        call def_var_info(surf( 5),trim(file_surface),"toa_net_all_clim","toa_net_all",units="W m**-2",method="radius")
        call def_var_info(surf( 6),trim(file_surface),"toa_net_clr_clim","toa_net_clr",units="W m**-2",method="radius")
        call def_var_info(surf( 7),trim(file_surface),"toa_cre_sw_clim","toa_cre_sw",units="W m**-2",method="radius")
        call def_var_info(surf( 8),trim(file_surface),"toa_cre_lw_clim","toa_cre_lw",units="W m**-2",method="radius")
        call def_var_info(surf( 9),trim(file_surface),"toa_cre_net_clim","toa_cre_net",units="W m**-2",method="radius")
        call def_var_info(surf(10),trim(file_surface),"solar_clim","solar",units="W m**-2",method="radius")

        ! Allocate the input grid variable
        call grid_allocate(gCERES,invar)

        ! Initialize mapping
        call map_init(map,gCERES,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
    
        ! Loop over months and map climatological gridded variables
        nm = 12 
        do m = 1, nm 

            write(*,*) "month ",m

            do i = 1, size(surf)
                var_now = surf(i)
                call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
                             start=[1,1,m],count=[gCERES%G%nx,gCERES%G%ny,1])
                call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method, &
                              fill=.TRUE.,missing_value=missing_value)
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month", &
                              units=var_now%units_out,start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
            end do 
        end do 

        return 

    end subroutine CERES_to_grid

    subroutine MARv33_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,clim_range)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MAR (RCM) DATA - MARv3.3 downloaded from the ftp site:
        !       ftp://ftp.climato.be/fettweis/MARv3.3/Greenland
        !       Data is available on the Bamber et al. (2001) 5km grid
        !       domain="Greenland-ERA": ERA-40 + ERA-Interim combined datasets
        !       domain="Greenland-MIROC5-RCP85": MIROC5 histo+rcp85 combined datasets
        !
        ! =========================================================

        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer, optional :: max_neighbors 
        double precision, optional :: lat_lim 
        integer, optional :: clim_range(2)

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
        integer :: yearf, k0, nk 
        character(len=512) :: filename_clim 
        double precision, allocatable :: var3D(:,:,:), var2D(:,:)

        ! Define input grid
        if (trim(domain) .eq. "Greenland-ERA") then 
            
            ! Define MAR (Bamber et al. 2001) grid and input variable field
            call grid_init(gMAR,name="Bamber01-5KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                           x0=-800.d0,dx=5.d0,nx=301,y0=-3400.d0,dy=5.d0,ny=561, &
                           lambda=-39.d0,phi=90.d0,alpha=7.5d0)

            ! Define the input filenames
            file_invariant = "/data/sicopolis/data/MARv3.3/Greenland/ERA_1958-2013_15km/"// &
                             "MARv3.3-15km-monthly-ERA-Interim-2013.nc"
            file_surface   = "/data/sicopolis/data/MARv3.3/Greenland/"
            file_prefix(1) = "ERA_1958-2013_15km/MARv3.3-15km-monthly-ERA-Interim-"
            file_prefix(2) = "ERA_1958-2013_15km/MARv3.3-15km-monthly-ERA-Interim-"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_MARv3.3-15km-monthly-ERA-Interim_195801-201312.nc"

            year0       = 1958 
            year_switch = 1979   ! Switch scenarios (ERA-40 to ERA-INTERIM)
            nyr         = 2013-1958+1

            ! For climatology
            if (present(clim_range)) then  
                k0 = clim_range(1) - year0+1
                nk = clim_range(2) - clim_range(1) + 1 

                write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
                    "_MARv3.3-15km-monthly-ERA-Interim_",clim_range(1),"-",clim_range(2),".nc"
            end if 

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

            ! For climatology
            if (present(clim_range)) then  
                k0 = clim_range(1) - year0+1
                nk = clim_range(2) - clim_range(1) + 1 

                write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
                    "_MARv3.3-15km-monthly-MIROC5-rcp85_",clim_range(1),"-",clim_range(2),".nc"
            end if 

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

        allocate(surf(19))
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
        call def_var_info(surf(17),trim(file_surface),"QQ",  "q",   units="g kg**-1")
        call def_var_info(surf(18),trim(file_surface),"CC",  "cc",  units="(0 - 1)")
        call def_var_info(surf(19),trim(file_surface),"SH3", "SH3", units="mm d**-1",conv=1d3*12.d0/365.d0)   ! m/month => mm/day

        nm       = 12
        n_var    = size(surf)
        if (trim(domain) .ne. "Greenland-ERA") n_var = 16   ! Exclude QQ, CC and SH3 if not available
    
        if (present(max_neighbors) .and. present(lat_lim)) then 

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
                call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value)
                outvar = missing_value 
                call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,100.d3, &
                               fill=.FALSE.,missing_value=missing_value)
                if (var_now%method .eq. "nn") then 
                    call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
                else
                    call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
                end if  
            end do 

            n_prefix = 1 
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
                        write(var_now%filename,"(a,a,i4,a3)") &
                            trim(adjustl(file_surface)), trim(file_prefix(n_prefix)),year,".nc"
                        call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
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
        
            ! ## INVARIANT FIELDS ##
            do i = 1, size(invariant)
                var_now = invariant(i) 
                call nc_read(filename,var_now%nm_out,var2D)
                call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc", &
                              units=var_now%units_out)
            end do 

            do i = 1, n_var
                var_now = surf(i)

                do m = 1, nm  
                    call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk])
                    var2D = time_average(var3D)
                    call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                  units=var_now%units_out,start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
                end do 
            end do 

        end if 

        return 

    end subroutine MARv33_to_grid

    subroutine MARv32_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MAR (RCM) DATA - MARv3.2 original data passed by Xavier
        !
        ! =========================================================

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
            call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value)
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
                        call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
                                      start=[1,1,1,q],count=[gMAR%G%nx,gMAR%G%ny,1,1])
                    else 
                        call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
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



    !##############################################
    !
    ! General subroutines related to the module
    !
    !##############################################



    ! Define some variable info for later manipulation
    subroutine def_var_info(var,filename,nm_in,nm_out,units,method,mask,dimextra,conv,plev,filenames)
        implicit none 

        type(var_defs) :: var 
        character(len=*) :: filename,nm_in,nm_out,units
        character(len=*), optional :: method 
        logical, optional :: mask, dimextra
        character(len=*), optional :: plev 
        character(len=*), optional :: filenames(:)
        double precision, optional :: conv 

        var%filename  = trim(filename)
        var%nm_in     = trim(nm_in)
        var%nm_out    = trim(nm_out)
        var%units_in  = trim(units)
        var%units_out = trim(units)

        var%method = "shepard"
        if (present(method)) var%method = trim(method)

        var%mask = .FALSE. 
        if (present(mask)) var%mask = mask 

        var%dimextra = .FALSE.
        if (present(dimextra)) var%dimextra = dimextra 

        var%plev = "None"
        if (present(plev)) var%plev = trim(plev)

        var%filenames(:) = "None"
        if (present(filenames)) var%filenames = filenames

        var%conv = 1.d0 
        if (present(conv)) var%conv = conv 

        return 

    end subroutine def_var_info

    ! Extract a thinner version of an input array
    ! (new array should be a multiple of input array)
    subroutine thin(var1,var,by)
        implicit none

        double precision, dimension(:,:) :: var, var1 
        integer :: by 
        integer :: i,j, nx, ny 
        integer :: i1, j1

        nx = size(var,1)
        ny = size(var,2) 

        var1 = missing_value 

        i1 = 0
        do i = 1, nx, by 
            i1 = i1+1 
            j1 = 0 
            do j = 1, ny, by  
                j1 = j1 + 1 
                var1(i1,j1) = var(i,j)
            end do 
        end do 

        return
    end subroutine thin 

    ! Fill in missing values of an array with neighbor averages
    ! or with a specified fill_value
    subroutine fill(var,missing_value,fill_value)
        implicit none 
        double precision, dimension(:,:) :: var 
        double precision :: missing_value 
        double precision, optional :: fill_value

        integer :: q, nx, ny, i, j 
        integer, parameter :: qmax = 50 ! Iterations 

        double precision, dimension (3,3) :: neighb, weight
        double precision :: wtot, mval 
        double precision, dimension(:,:), allocatable :: filled
        nx = size(var,1)
        ny = size(var,2) 

        allocate(filled(nx,ny))

        if (present(fill_value)) then
            where(var .eq. missing_value) var = fill_value 
        end if 

        do q = 1, qmax 

            filled = missing_value 

            do i = 2, nx-1 
                do j = 2, ny-1 
                    neighb = var(i-1:i+1,j-1:j+1)

                    weight = 0.d0 
                    where (neighb .ne. missing_value) weight = 1.d0
                    wtot = sum(weight)

                    if (wtot .gt. 0.d0) then 
                        mval = sum(neighb*weight)/wtot
                        where (neighb .eq. missing_value) neighb = mval 
                    end if 

                    filled(i-1:i+1,j-1:j+1) = neighb 

                end do 
            end do 

            where(filled .ne. missing_value) var = filled 

            write(*,*) q," : Missing values: ", count(var .eq. missing_value)
            if ( count(var .eq. missing_value) .eq. 0 ) exit 
        end do 

        return
    end subroutine fill 

end module gridding_datasets

