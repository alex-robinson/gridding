module CERES 

    use gridding_datasets
    use coordinates 
    use ncio 
    
    implicit none 

    private 
    public :: CERES_to_grid
    
contains 

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


end module CERES 

