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
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid_in
        character(len=256) :: file_in
        type(var_defs), allocatable :: vars(:)
        double precision, allocatable :: invar(:,:) 
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        ! Define ECMWF input grid
        if (trim(domain) .eq. "Global") then 
            
            ! Define MAR grid and input variable field
            call grid_init(grid_in,name="CERES-1deg",mtype="latlon",units="degrees",lon180=.FALSE., &
                   x0=0.5d0,dx=1d0,nx=360,y0=-90.d0,dy=1d0,ny=180 )

            ! Define the input filenames
            file_in = "/data/sicopolis/data/CERES/CERES_EBAF-TOA_Ed2.8_Subset_CLIM01-CLIM12.nc"
            desc    = "Clouds and the Earth's radiant energy system (CERES) EBAF-TOA products"
            ref     = "NASA's Earth Observing System, &
                      &http://ceres.larc.nasa.gov/products.php?product=EBAF-TOA"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_CERES_2001-2013.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(vars(10))
        call def_var_info(vars( 1),trim(file_in),"toa_sw_all_clim","toa_sw_all",units="W m**-2", &
                          long_name="TOA shortwave radiation (all)",method="nng")
        call def_var_info(vars( 2),trim(file_in),"toa_sw_clr_clim","toa_sw_clr",units="W m**-2", &
                          long_name="TOA shortwave radiation (clear)",method="nng")
        call def_var_info(vars( 3),trim(file_in),"toa_lw_all_clim","toa_lw_all",units="W m**-2", &
                          long_name="TOA longwave radiation (all)",method="nng")
        call def_var_info(vars( 4),trim(file_in),"toa_lw_clr_clim","toa_lw_clr",units="W m**-2", &
                          long_name="TOA longwave radiation (clear)",method="nng")
        call def_var_info(vars( 5),trim(file_in),"toa_net_all_clim","toa_net_all",units="W m**-2", &
                          long_name="TOA net radiation (all)",method="nng")
        call def_var_info(vars( 6),trim(file_in),"toa_net_clr_clim","toa_net_clr",units="W m**-2", &
                          long_name="TOA net radiation (clear)",method="nng")
        call def_var_info(vars( 7),trim(file_in),"toa_cre_sw_clim","toa_cre_sw",units="W m**-2", &
                          long_name="TOA cre shortwave clim",method="nng")
        call def_var_info(vars( 8),trim(file_in),"toa_cre_lw_clim","toa_cre_lw",units="W m**-2", &
                          long_name="TOA cre longwave clim",method="nng")
        call def_var_info(vars( 9),trim(file_in),"toa_cre_net_clim","toa_cre_net",units="W m**-2", &
                          long_name="TOA cre net clim",method="nng")
        call def_var_info(vars(10),trim(file_in),"solar_clim","solar",units="W m**-2", &
                          long_name="Solar clim",method="nng")

        ! Allocate the input grid variable
        call grid_allocate(grid_in,invar)

        ! Initialize mapping
        call map_init(map,grid_in,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ## Map climatological gridded variables ##
        
        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            ! Loop over months
            do m = 1, 12 
                call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
                             start=[1,1,m],count=[grid_in%G%nx,grid_in%G%ny,1])
                call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method, &
                              fill=.TRUE.,sigma=20.d0,missing_value=missing_value)
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month", &
                              start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
            end do 

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            !             call nc_write_attr(filename,var_now%nm_out,"grid_mapping",trim(grid%mtype))
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine CERES_to_grid


end module CERES 

