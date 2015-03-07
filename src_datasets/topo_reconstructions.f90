module topo_reconstructions 

    use gridding_datasets
    use coordinates 
    use ncio 
    
    implicit none 

    private 
    public :: ICE6GC_to_grid
    
contains 

    subroutine ICE6GC_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       ICE-6G_C DATA
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
            double precision, allocatable :: var(:,:)
        end type 

        type(inp_type)     :: inp
        integer :: nx, ny 
        type(grid_class)   :: grid0
        character(len=256) :: fldr_in, file_in, file_prefix, times(48)
        type(var_defs), allocatable :: vars(:)
        double precision :: time 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var, kt 

        ! Define the input filenames
        fldr_in         = "/data/sicopolis/data/ICE-6G_C/"
        file_in         = trim(fldr_in)//"I6_C.VM5a_1deg.0.nc"
        file_prefix     = trim(fldr_in)//"I6_C.VM5a_1deg."
        times = ["0   ","0.5 ","1   ","1.5 ","2   ","2.5 ","3   ","3.5 ","4   ","4.5 ", &
                 "5   ","5.5 ","6   ","6.5 ","7   ","7.5 ","8   ","8.5 ","9   ","9.5 ", &
                 "10  ","10.5","11  ","11.5","12  ","12.5","13  ","13.5","14  ","14.5", &
                 "15  ","15.5","16  ","16.5","17  ","17.5","18  ","18.5","19  ","19.5", &
                 "20  ","20.5","21  ","22  ","23  ","24  ","25  ","26  "]

        desc    = "ICE-6G_C reconstructed paleo Earth topography"
        ref     = "Argus, D.F., Peltier, W.R., Drummond, R. and Moore, A.W.: &
                  &The Antarctica component of postglacial rebound model &
                  &ICE-6G_C (VM5a) based upon GPS positioning, exposure &
                  &age dating of ice thicknesses, and relative sea level &
                  &histories. Geophys. J. Int., 198, 537-563, 2014. \n\n &
                  &Peltier, W. R., Argus, D. F. and Drummond, R.: &
                  &Space geodesy constrains ice age terminal deglaciation: &
                  &The global ICE-6G_C (VM5a) model, J. Geophys. Res. Solid Earth, &
                  &120, 450–487, doi:10.1002/2014JB011176, 2015. \n\n &
                  &http://www.atmosp.physics.utoronto.ca/~peltier/data.php"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_TOPO-ICE-6G_C.nc"

        ! Define the input grid
        call grid_init(grid0,name="GLOBAL-1DEG",mtype="latlon",units="degrees",lon180=.TRUE., &
                       x0=0.5d0,dx=1.d0,nx=360,y0=-89.5d0,dy=1.d0,ny=180)
        
        ! Allocate the input array
        call grid_allocate(grid0,inp%var)

        ! Define the variables to be mapped 
        allocate(vars(5))
        call def_var_info(vars( 1),trim(file_prefix),"sftlf","sftlf",units="%", &
                          long_name="Land area fraction",method="quadrant")
        call def_var_info(vars( 2),trim(file_prefix),"sftgif","sftgif",units="%", &
                          long_name="Ice area fraction",method="quadrant")
        call def_var_info(vars( 3),trim(file_prefix),"Topo_Diff","dz",units="m", &
                          long_name="Topography difference from present",method="quadrant")
        call def_var_info(vars( 4),trim(file_prefix),"Topo","z",units="m", &
                          long_name="Topography (Point-value altitude)",method="quadrant")
        call def_var_info(vars( 5),trim(file_prefix),"orog","zs",units="m", &
                          long_name="Orography (Point-value surface altitude)",method="quadrant")

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)     

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"time",x=0.d0,units="kiloyears",unlimited=.TRUE.)
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)
            kt = size(times)-k+1
            
            do k = 1, size(times)
                file_in = trim(var_now%filename)//trim(times(kt))//".nc"
                read(times(kt),*) time
                time = -time   ! negative time 

                ! Read in current variable
                call nc_read(file_in,var_now%nm_in,inp%var,missing_value=missing_value)
                where(abs(inp%var) .ge. 1d8) inp%var = missing_value 

                ! Map variable to new grid
                call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
                              fill=.TRUE.,missing_value=missing_value)

                ! Write output variable to output file
                call nc_write(filename,"time",time,dim1="time",start=[k],count=[1])
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="time", &
                              start=[1,1,k],count=[grid%G%nx,grid%G%ny,1])

            end do 

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine ICE6GC_to_grid


end module topo_reconstructions 

