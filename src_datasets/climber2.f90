module climber2 

    use gridding_datasets
    use coordinates 
    use interp2D 
    use gaussian_filter
    use ncio 
    
    implicit none 

    private 
    public :: climber2_atm_to_grid

    ! ocean directory: /home/jalvarez/outputs_clim3/ocean_marisa

contains 

    subroutine climber2_atm_to_grid(outfldr,subfldr,grid,sim,path_in,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CLIMBER2 atmospheric DATA
        !
        !   For Greenland, eg: 
        !   m=1:12 (months jan-dec), n=1 (atlantic sector),
        !   i=16:17 (lats 60-70,70-80), ntyp=2 (temperature)
        ! =========================================================
        implicit none 

        character(len=*) :: sim, outfldr, subfldr, path_in 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: sigma, lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:), var(:,:)
            double precision, allocatable :: var0(:,:,:)
            double precision, allocatable :: time(:)
            double precision, allocatable :: zs(:,:) 
            double precision :: lapse = 6.5d-3 
        end type 

        type(inp_type)     :: inp
        type(grid_class)   :: grid0
        character(len=256) :: fldr_in, file_in
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, np, nt

        type(map_class)  :: map
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), outzs(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        double precision, parameter :: sigma_climber2 = 1000.d0 

        ! Define the input filenames
        fldr_in      = trim(path_in)
        file_in      = trim(fldr_in)//"a2g_c_small.nc"

        desc    = "CLIMBER-2 simulation output"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(subfldr)//"/"// &
                              trim(grid%name)//"_"//trim(sim)//".nc"

        ! Load the domain information 
        nx =  7
        ny = 18
        np = nx*ny 

        nt = nc_size(file_in,"time")

        allocate(inp%time(nt))
        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))
        allocate(inp%var0(nx,ny,12))
        allocate(inp%zs(nx,ny))

        call nc_read(file_in,"time",inp%time)
        call nc_read(file_in,"lon",inp%lon,start=[1],count=[7])   ! 8-11 are zonal averages
        call nc_read(file_in,"lat",inp%lat)
        
        ! Longitude in file are sectors (1-7), specify actual longitudes here
        inp%lon(1) = -40.d0
        do i = 2, 5 
            inp%lon(i) = inp%lon(1) + (i-1)*51.4d0
        end do 
        inp%lon(7) = inp%lon(1) - (1)*51.4d0 
        inp%lon(6) = inp%lon(1) - (2)*51.4d0 

        ! Shift longitude so that it goes from negative to positive (easier for netcdf viewing)
        inp%lon = cshift(inp%lon,-2)

        write(*,"(a,20f10.2)") "lon: ", inp%lon 
        write(*,"(a,20f10.2)") "lat: ", inp%lat 

        ! Define climber2 points and input variable field
        call grid_init(grid0,name="climber2-atmos",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat )
        call grid_write(grid0,"output/"//trim(grid0%name)//".nc",xnm="xc",ynm="yc",create=.TRUE.) 

        ! Define the variables to be mapped 
        allocate(vars(2))
        call def_var_info(vars( 1),trim(file_in),"ts","t2m_sl",units="degrees Celcius", &
                          long_name="Near-surface temperature at sea level",method="nn")
        call def_var_info(vars( 2),trim(file_in),"prc","pr",units="mm*d**-1", &
                          long_name="Precipitation",method="nn")

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        call grid_allocate(grid,outzs) 
        
        ! Initialize mappings
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"month",x=1,nx=12,units="1")
        call nc_write_dim(filename,"time", x=inp%time,units="kiloyears")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Load reference topography in order to adjust temps to sea-level temps 
        call nc_read(file_in,"horo",inp%zs,missing_value=mv,start=[1,1,1,1],count=[nx,ny,1,1])
        inp%zs = cshift(inp%zs,-2,dim=1)

        ! Map zs to new grid
        call map_field(map,"zs",inp%zs,outzs,outmask,"nn", &
                          fill=.TRUE.,missing_value=mv,sigma=sigma)
        call climber2_smooth(outzs,sigma=sigma_climber2,dx=grid%G%dx,mask=outzs .ne. mv)
        call nc_write(filename,"zs",real(outzs),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,"zs","units","m")
        call nc_write_attr(filename,"zs","long_name","Surface elevation")
        call nc_write_attr(filename,"zs","coordinates","lat2D lon2D")
        
        ! Also generate a land mask
        inp%var = 0.d0
        where(inp%zs .gt. 1.d0) inp%var = 1.d0 

        ! Map mask to new grid
        call map_field(map,"mask_land",inp%var,outvar,outmask,"nn", &
                          fill=.TRUE.,missing_value=mv,sigma=sigma)
        call nc_write(filename,"mask_land",nint(outvar),dim1="xc",dim2="yc")

        ! Write variable metadata
        call nc_write_attr(filename,"mask_land","units","1")
        call nc_write_attr(filename,"mask_land","long_name","Land mask (land=1)")
        call nc_write_attr(filename,"mask_land","coordinates","lat2D lon2D")

        ! ## Map climatological gridded variables ##
        
        ! Load present day temperature field 
        var_now = vars(1)
        call nc_read(trim(var_now%filename),var_now%nm_in,inp%var0,missing_value=mv, &
                             start=[1,1,1,nt],count=[nx,ny,12,1])

        do m = 1, 12 
            inp%var0 = cshift(inp%var0,-2,dim=1)
            inp%var0(:,:,m) = inp%var0(:,:,m) + inp%lapse*inp%zs 
        end do 

        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            do k = 1, nt 
            do m = 1, 12 

                ! Read in current variable
                call nc_read(trim(var_now%filename),var_now%nm_in,inp%var,missing_value=mv, &
                             start=[1,1,m,k],count=[nx,ny,1,1])
                inp%var = cshift(inp%var,-2,dim=1)

                ! Scale to sea-level temperature for interpolation
                if (trim(var_now%nm_out) .eq. "t2m_sl") &
                    inp%var = inp%var + inp%lapse*inp%zs 

                ! Map variable to new grid
                call map_field(map,var_now%nm_in,inp%var,outvar,outmask,"nn", &
                              fill=.TRUE.,missing_value=mv,sigma=sigma)
                call climber2_smooth(outvar,sigma=sigma_climber2,dx=grid%G%dx,mask=outvar .ne. mv)

                ! Write output variable to output file
                call nc_write(filename,var_now%nm_out,real(outvar), &
                              dim1="xc",dim2="yc",dim3="month",dim4="time", &
                              start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])

!                 ! Also map anomalies
!                 if (trim(var_now%nm_out) .eq. "t2m_sl") then 
!                     inp%var = inp%var - inp%var0(:,:,m)

!                     ! Map variable to new grid
!                     call map_field(map,var_now%nm_in,inp%var,outvar,outmask,"nn", &
!                                   fill=.TRUE.,missing_value=mv,sigma=sigma)
!                     call climber2_smooth(outvar,sigma=sigma_climber2,dx=grid%G%dx,mask=outvar .ne. mv)

!                     ! Write output variable to output file
!                     call nc_write(filename,"dt2m_sl",real(outvar), &
!                                   dim1="xc",dim2="yc",dim3="month",dim4="time", &
!                                   start=[1,1,m,k],count=[grid%G%nx,grid%G%ny,1,1])

!                 end if 

            end do     
            end do

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine climber2_atm_to_grid

    subroutine climber2_smooth(var2D,sigma,dx,mask)
        ! Perform gaussian smoothing manually because
        ! climber2 needs a very large sigma (ie multiple iterations)

        implicit none 

        double precision, intent(INOUT) :: var2D(:,:)
        logical          :: mask(:,:)
        double precision :: sigma, dx 

        ! Local variables
        double precision :: tmp2D(size(var2D,1),size(var2D,2))
        double precision, parameter :: sigma_max = 250.d0 
        integer :: nloop, i 

        if (sigma .lt. sigma_max) then 
            write(*,*) "climber2_smooth:: error: desired sigma is too small for the climber2 &
                       &atmospheric grid. Try again."
            stop 
        end if 

        nloop = (sigma / sigma_max)**2 

        write(*,*) "climber2_smooth nloop = ", nloop

        do i = 1, nloop 
            tmp2D = var2D
            call filter_gaussian(input=tmp2D,output=var2D,sigma=sigma_max,dx=dx,mask=mask)
        end do 

        return 

    end subroutine climber2_smooth
    
end module climber2 

