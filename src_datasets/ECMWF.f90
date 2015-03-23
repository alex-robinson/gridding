module ECMWF 

    use gridding_datasets
    use coordinates
    use interp_time 
    use ncio 
    
    implicit none 

    private 
    public :: ecmwf_to_grid
    
contains 

    subroutine ecmwf_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,clim_range)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       ECMWF DATA (ERA-INTERIM 1979-2013)
        !
        ! ========================================================= 

        implicit none 

        character(len=*)    :: domain, outfldr 
        type(grid_class)    :: grid 
        integer, optional   :: max_neighbors 
        double precision, optional :: lat_lim 
        integer, optional   :: clim_range(2) 
        character(len=512)  :: filename 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:)
            double precision, allocatable :: var(:,:)
        end type 

        type(inp_type)     :: inp

        type(grid_class)   :: gECMWF 
        character(len=256) :: file_invariant, file_surface, files_pres(9)
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
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

!         ! Define ECMWF input grid
!         if (trim(domain) .eq. "GRL075") then 
            
!             ! Initialize the grid
!             call grid_init(gECMWF,name="ECMWF-GRL075",mtype="latlon",units="kilometers",lon180=.TRUE., &
!                            x0=-100.d0,dx=0.75d0,nx=161,y0=49.5d0,dy=0.75d0,ny=55)
            
!             ! Assign the filenames
!             file_invariant = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-invariant_historical_mon_197901-201212.nc"
!             file_surface   = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-surface_historical_mon_197901-201212.nc"
!             files_pres(1)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-1000Mb_historical_mon_197901-201212.nc"
!             files_pres(2)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-950Mb_historical_mon_197901-201212.nc"
!             files_pres(3)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-850Mb_historical_mon_197901-201212.nc"
!             files_pres(4)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-750Mb_historical_mon_197901-201212.nc"
!             files_pres(5)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-700Mb_historical_mon_197901-201212.nc"
!             files_pres(6)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-650Mb_historical_mon_197901-201212.nc"
!             files_pres(7)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-600Mb_historical_mon_197901-201212.nc"
!             files_pres(8)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-550Mb_historical_mon_197901-201212.nc"
!             files_pres(9)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-GRL-500Mb_historical_mon_197901-201212.nc"

!             desc    = "Regional Greenland subset of the ERA-Interim dataset"
!             ref     = "Dee et al., 2011, BAMS, &
!                       &http://onlinelibrary.wiley.com/doi/10.1002/qj.828/abstract"

!         else if (trim(domain) .eq. "ANT075") then 

!             ! Initialize the grid
!             call grid_init(gECMWF,name="ECMWF-ANT075",mtype="latlon",units="kilometers",lon180=.TRUE., &
!                            x0=-180.d0,dx=0.75d0,nx=480,y0=-90.d0,dy=0.75d0,ny=67)
            
!             ! Assign the filenames
!             file_invariant = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-invariant_historical_mon_197901-201212.nc"
!             file_surface   = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-surface_historical_mon_197901-201212.nc"
!             files_pres(1)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-1000Mb_historical_mon_197901-201212.nc"
!             files_pres(2)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-950Mb_historical_mon_197901-201212.nc"
!             files_pres(3)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-850Mb_historical_mon_197901-201212.nc"
!             files_pres(4)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-750Mb_historical_mon_197901-201212.nc"
!             files_pres(5)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-700Mb_historical_mon_197901-201212.nc"
!             files_pres(6)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-650Mb_historical_mon_197901-201212.nc"
!             files_pres(7)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-600Mb_historical_mon_197901-201212.nc"
!             files_pres(8)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-550Mb_historical_mon_197901-201212.nc"
!             files_pres(9)  = "/data/sicopolis/data/ECMWF/ERA-INTERIM-ANT-500Mb_historical_mon_197901-201212.nc"

!             desc    = "Regional Antarctica subset of the ERA-Interim dataset"
!             ref     = "Dee et al., 2011, BAMS, &
!                       &http://onlinelibrary.wiley.com/doi/10.1002/qj.828/abstract"

!     else ! GLOBAL DOMAIN

            ! Assign the filenames
            file_invariant = "/data/sicopolis/data/ECMWF/ERA-INTERIM-invariant_historical_mon_197901-201212.nc"
            file_surface   = "/data/sicopolis/data/ECMWF/ERA-INTERIM-surface_historical_mon_197901-201212.nc"
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
            call nc_read(file_invariant,"latitude",inp%lon)
            call nc_read(file_invariant,"latitude",inp%lat)
            inp%lat = inp%lat(ny:1)

            call grid_init(gECMWF,name="ECMWF-075",mtype="latlon",units="kilometers",lon180=.TRUE., &
                           x=inp%lon,y=inp%lat)
            
            write(*,*) "GRID CHECK: "
            write(*,*) gECMWF%G%x(1), gECMWF%G%x(nx)
            write(*,*) gECMWF%G%y(1), gECMWF%G%y(ny)
              
!         end if 

        ! ## First make file for surface fields including invariants ##
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_ERA-INTERIM_197901-201212.nc"

        ! For climatology
        if (present(clim_range)) then  
            k0 = clim_range(1) - 1979 + 1
            nk = clim_range(2) - clim_range(1) + 1 

            write(filename_clim,"(a,i4,a1,i4,a3)") trim(outfldr)//"_clim/"//trim(grid%name)// &
                "_ERA-INTERIM_",clim_range(1),"-",clim_range(2),".nc"
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
                         start=[1,1,1],count=[gECMWF%G%nx,gECMWF%G%ny,1])
            call flip_lat(invar)
            invar = invar*var_now%conv 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,"nn",sigma=20.d0,missing_value=missing_value)
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",units=var_now%units_out)






            stop 






            ! ## SURFACE FIELDS ##
            do i = 1, size(surf)
                var_now = surf(i)

                q = 0 
                do k = 1, nyr 

                    year = 1978 + k 
                    write(*,*) trim(var_now%nm_in)," :",year

                    do m = 1, nm 
                        q = q+1 
                        call nc_read(trim(var_now%filename),var_now%nm_in,invar,start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1])
                        call flip_lat(invar)
                        call map_field(map,var_now%nm_in,invar,outvar,outmask,"nn",sigma=20.d0,missing_value=missing_value)
                        call nc_write(filename,var_now%nm_out,real(outvar),  dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                      units=var_now%units_out,start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])
                    end do 
                end do
                
                ! Write variable metadata
                call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
!                 call nc_write_attr(filename,var_now%nm_out,"grid_mapping",trim(grid%mtype))
                call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
                
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
                                         start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1])
                            call flip_lat(invar)
                            call map_field(map,var_now%nm_in,invar,outvar,outmask,"nn",sigma=20.d0, &
                                           missing_value=missing_value)
                            call nc_write(filename,var_now%nm_out,real(outvar), &
                                          dim1="xc",dim2="yc",dim3="month",dim4="time", &
                                          units=var_now%units_out,start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])

                        end do 
                    end do 

                    ! Write variable metadata
                    call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
                    call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
!                     call nc_write_attr(filename,var_now%nm_out,"grid_mapping",trim(grid%mtype))
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
            call nc_read(filename,var_now%nm_out,var2D)
            call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc", &
                          units=var_now%units_out)

            do i = 1, size(surf)
                var_now = surf(i)
                do m = 1, nm  
                    call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk])
                    var2D = time_average(var3D)
                    call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                  units=var_now%units_out,start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
                end do 

                ! Write variable metadata
                call nc_write_attr(filename_clim,var_now%nm_out,"units",var_now%units_out)
                call nc_write_attr(filename_clim,var_now%nm_out,"long_name",var_now%long_name)
!                 call nc_write_attr(filename_clim,var_now%nm_out,"grid_mapping",trim(grid%mtype))
                call nc_write_attr(filename_clim,var_now%nm_out,"coordinates","lat2D lon2D")
                
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
                call grid_write(grid,filename_clim,xnm="xc",ynm="yc",create=.FALSE.)
                
                ! Write meta data 
                call nc_write_attr(filename_clim,"Description",desc)
                call nc_write_attr(filename_clim,"Reference",ref)

                do i = 1, size(pres)
                    var_now = pres(i)

                    do m = 1, nm  
                        call nc_read(filename,var_now%nm_out,var3D,start=[1,1,m,k0],count=[grid%G%nx,grid%G%ny,1,nk])
                        var2D = time_average(var3D)
                        call nc_write(filename_clim,var_now%nm_out,real(var2D),dim1="xc",dim2="yc",dim3="month", &
                                      units=var_now%units_out,start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
                    end do 

                    ! Write variable metadata
                    call nc_write_attr(filename_clim,var_now%nm_out,"units",var_now%units_out)
                    call nc_write_attr(filename_clim,var_now%nm_out,"long_name",var_now%long_name)
!                     call nc_write_attr(filename_clim,var_now%nm_out,"grid_mapping",trim(grid%mtype))
                    call nc_write_attr(filename_clim,var_now%nm_out,"coordinates","lat2D lon2D")
                
                end do 

            end do 

        end if 

        return 

    end subroutine ecmwf_to_grid


    subroutine flip_lat(var)

        implicit none 
        double precision :: var(:,:)
        integer :: ny 

        ny = size(var,2)
        var = var(:,ny:1)

        return 

    end subroutine flip_lat 

end module ECMWF 
