
module gridding_datasets

    use coordinates 
    use ncio 
    use vargrid 

    implicit none 

    private
    public :: ecmwf_to_grid

contains

    subroutine ecmwf_to_grid(filename,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file

        implicit none 

        character(len=*) :: domain, filename 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 

        type(grid_class)   :: gECMWF 
        character(len=256) :: file_invariant, file_surface, files_pres(9)
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:) 

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

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"plev", x=[1000.d0,950.d0,850.d0,750.d0,700.d0,650.d0,600.d0,550.d0,500.d0],units="hPa")
        call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call nc_write_dim(filename,"time", x=1979,dx=1,nx=34,units="years",calendar="360_day")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
    
        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    

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
            write(*,*) 
            write(*,*) "=== ",year," ==="
            write(*,*)

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

                ! ## PRESSURE FIELDS ##
                do i = 1, size(pres)
                    var_now = pres(i) 

                    do l = 1, size(files_pres)   ! Loop over pressure layers
                        call nc_read(var_now%filenames(l),var_now%nm_in,invar,start=[1,1,q],count=[gECMWF%G%nx,gECMWF%G%ny,1])
                        call map_field(map,var_now%nm_in,invar,outvar,outmask,"shepard",400.d3,missing_value=missing_value)
                        call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="plev",dim4="month", &
                                      dim5="time",units=var_now%units_out,start=[1,1,l,m,k],count=[grid%G%nx,grid%G%ny,1,1,1])
                    end do 

                end do 
            end do 
        end do 


        return 

    end subroutine ecmwf_to_grid




end module gridding_datasets

