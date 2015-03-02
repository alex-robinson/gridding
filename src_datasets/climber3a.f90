module climber3a 

    use gridding_datasets
    use coordinates 
    use ncio 
    
    implicit none 

    private 
    public :: climber3a_to_grid
    
contains 

    subroutine climber3a_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       CLIMBER3-alpha DATA
        !
        ! =========================================================
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type inp_type 
            double precision, allocatable :: lon(:), lat(:), var(:,:)
        end type 

        type(inp_type)     :: inp
        type(grid_class)   :: gTOPO
        character(len=256) :: fldr_in, file_in_topo, file_in 
        type(var_defs), allocatable :: vars(:)
        integer :: nx, ny, np 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var 

        ! Define the input filenames
        fldr_in      = "/data/sicopolis/data/CLIMBER3a/brutos/"
        file_in_topo = trim(fldr_in)//"horo.cdf"
        file_in      = trim(fldr_in)//trim(domain)//".cdf"

        desc    = "CLIMBER-3alpha simulation output"
        ref     = "source folder: "//trim(fldr_in)

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_"//trim(domain)//".nc"

        ! Load the domain information 
        nx = nc_size(file_in_topo,"XT_I")
        ny = nc_size(file_in_topo,"YT_J")
        np = nx*ny 

        allocate(inp%lon(nx),inp%lat(ny),inp%var(nx,ny))

        call nc_read(file_in_topo,"XT_I",inp%lon)
        call nc_read(file_in_topo,"YT_J",inp%lat)
        
        ! Define CLIMBER3a points and input variable field
        call grid_init(gTOPO,name="climber3a-atmos",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat )

        
        stop 

        ! Define the variables to be mapped 
        allocate(vars(10))
        call def_var_info(vars( 1),trim(file_in),"toa_sw_all_clim","toa_sw_all",units="W m**-2", &
                          long_name="TOA shortwave radiation (all)",method="radius")
        call def_var_info(vars( 2),trim(file_in),"toa_sw_clr_clim","toa_sw_clr",units="W m**-2", &
                          long_name="TOA shortwave radiation (clear)",method="radius")
        call def_var_info(vars( 3),trim(file_in),"toa_lw_all_clim","toa_lw_all",units="W m**-2", &
                          long_name="TOA longwave radiation (all)",method="radius")
        call def_var_info(vars( 4),trim(file_in),"toa_lw_clr_clim","toa_lw_clr",units="W m**-2", &
                          long_name="TOA longwave radiation (clear)",method="radius")
        call def_var_info(vars( 5),trim(file_in),"toa_net_all_clim","toa_net_all",units="W m**-2", &
                          long_name="TOA net radiation (all)",method="radius")
        call def_var_info(vars( 6),trim(file_in),"toa_net_clr_clim","toa_net_clr",units="W m**-2", &
                          long_name="TOA net radiation (clear)",method="radius")
        call def_var_info(vars( 7),trim(file_in),"toa_cre_sw_clim","toa_cre_sw",units="W m**-2", &
                          long_name="TOA cre shortwave clim",method="radius")
        call def_var_info(vars( 8),trim(file_in),"toa_cre_lw_clim","toa_cre_lw",units="W m**-2", &
                          long_name="TOA cre longwave clim",method="radius")
        call def_var_info(vars( 9),trim(file_in),"toa_cre_net_clim","toa_cre_net",units="W m**-2", &
                          long_name="TOA cre net clim",method="radius")
        call def_var_info(vars(10),trim(file_in),"solar_clim","solar",units="W m**-2", &
                          long_name="Solar clim",method="radius")

        ! Initialize mapping
        call map_init(map,gTOPO,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
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
!                 call nc_read(trim(var_now%filename),var_now%nm_in,invar,missing_value=missing_value, &
!                              start=[1,1,m],count=[grid_in%G%nx,grid_in%G%ny,1])
!                 call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method, &
!                               fill=.TRUE.,missing_value=missing_value)
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="month", &
                              start=[1,1,m],count=[grid%G%nx,grid%G%ny,1])
            end do 

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine climber3a_to_grid


end module climber3a 

