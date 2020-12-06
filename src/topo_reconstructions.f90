module topo_reconstructions 

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: ICE6GC_to_grid
    public :: ICE5G_to_grid
    public :: LGMsimpson_to_grid 
    public :: huy3_to_grid
    public :: dated1_to_grid


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
                          long_name="Land area fraction",method="nng")
        call def_var_info(vars( 2),trim(file_prefix),"sftgif","sftgif",units="%", &
                          long_name="Ice area fraction",method="nng")
        call def_var_info(vars( 3),trim(file_prefix),"Topo_Diff","dz",units="m", &
                          long_name="Topography difference from present",method="nng")
        call def_var_info(vars( 4),trim(file_prefix),"Topo","z",units="m", &
                          long_name="Topography (Point-value altitude)",method="nng")
        call def_var_info(vars( 5),trim(file_prefix),"orog","zs",units="m", &
                          long_name="Orography (Point-value surface altitude)",method="nng")

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)     

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="km")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="km")
        call nc_write_dim(filename,"time",x=0.d0,units="kiloyears",unlimited=.TRUE.)
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            do k = 1, size(times)
                kt = size(times)-k+1
                file_in = trim(var_now%filename)//trim(times(kt))//".nc"
                read(times(kt),*) time
                time = -time   ! negative time 

                ! Read in current variable
                call nc_read(file_in,var_now%nm_in,inp%var,missing_value=mv)
                where(abs(inp%var) .ge. 1d8) inp%var = mv 

                ! Map variable to new grid
                call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
                              fill=.TRUE.,sigma=40.d0,missing_value=mv)

                ! Write output variable to output file
                call nc_write(filename,"time",time,dim1="time",start=[k],count=[1])
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="time", &
                              start=[1,1,k],count=[grid%G%nx,grid%G%ny,1],missing_value=real(mv))

            end do 

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine ICE6GC_to_grid

    subroutine ICE5G_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       ICE-5G DATA
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
        character(len=256) :: fldr_in, file_in, file_prefix, file_suffix
        type(var_defs), allocatable :: vars(:)
        double precision :: time 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        integer :: q, k, m, i, l, n_var
        character(len=10) :: time_char 
        real(4) :: times(39) 

        ! Define the input filenames
        fldr_in         = "/data/sicopolis/data/ICE-5G/ice5g_v1.2_0-21k_1deg/"
        file_in         = trim(fldr_in)//"ice5g_v1.2_00.0k_1deg.nc"
        file_prefix     = trim(fldr_in)//"ice5g_v1.2_"
        file_suffix     = "k_1deg.nc"

        times = [-21.0,-20.0,-19.0,-18.0,-17.0,-16.5,-16.0,-15.5,-15.0,-14.5,-14.0,-13.5,-13.0, &
                 -12.5,-12.0,-11.5,-11.0,-10.5,-10.0,-9.5,-9.0,-8.5,-8.0,-7.5,-7.0,-6.5,-6.0, &
                 -5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]

        desc    = "ICE-5G (VM2) model 90Km Lithosphere"
        ref     = "Peltier, D.: GLOBAL GLACIAL ISOSTASY AND THE SURFACE &
                  &OF THE ICE-AGE EARTH: The ICE-5G (VM2) Model and GRACE, &
                  &Ann. Rev. Earth Plan. Sci., 32, 111-149, &
                  &doi:10.1146/annurev.earth.32.082503.144359, 2004  \n\n &
                  &http://www.atmosp.physics.utoronto.ca/~peltier/data.php"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_TOPO-ICE-5G.nc"

        ! Define the input grid
        call grid_init(grid0,name="ICE5G-1DEG",mtype="latlon",units="degrees",lon180=.TRUE., &
                       x0=0.0d0,dx=1.d0,nx=360,y0=-89.5d0,dy=1.d0,ny=180)
        
        ! Allocate the input array
        call grid_allocate(grid0,inp%var)

        ! Define the variables to be mapped 
        allocate(vars(3))
        call def_var_info(vars( 1),trim(file_prefix),"sftgif","sftgif",units="%", &
                          long_name="Ice area fraction",method="nng")
        call def_var_info(vars( 2),trim(file_prefix),"orog","zs",units="m", &
                          long_name="Orography (surface altitude)",method="nng")
        call def_var_info(vars( 3),trim(file_prefix),"sftgit","H",units="m", &
                          long_name="Ice sheet thickness",method="nng")

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)     

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="km")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="km")
        call nc_write_dim(filename,"time",x=0.d0,units="kiloyears",unlimited=.TRUE.)
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Loop over variables
        do i = 1, size(vars)
            var_now = vars(i)

            do k = 1, size(times)
                time = times(k)

                if (abs(time).lt.10.0) then 
                    write(time_char,"(i1,f3.1)") 0, abs(time) 
                else 
                    write(time_char,"(f4.1)") abs(time) 
                end if 

                file_in = trim(file_prefix)//trim(time_char)//trim(file_suffix)

                ! Read in current variable
                call nc_read(file_in,var_now%nm_in,inp%var,missing_value=mv)
                where(abs(inp%var) .ge. 1d8) inp%var = mv 

                ! Map variable to new grid
                call map_field(map,var_now%nm_in,inp%var,outvar,outmask,var_now%method, &
                              fill=.TRUE.,sigma=80.d0,missing_value=mv)

                ! Write output variable to output file
                call nc_write(filename,"time",time,dim1="time",start=[k],count=[1])
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",dim3="time", &
                              start=[1,1,k],count=[grid%G%nx,grid%G%ny,1],missing_value=real(mv))

            end do 

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine ICE5G_to_grid

    subroutine LGMsimpson_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       SIMPSON GLACIAL MASK DATA
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
            double precision, allocatable :: lon(:), lat(:)
            double precision, allocatable :: var(:)
        end type 

        type(inp_type)     :: inp
        integer            :: nx, ny, np 
        type(points_class) :: points0
        character(len=256) :: fldr_in, file_in, file_in_grid

        type(map_class)  :: map 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        ! Define the input filenames
        fldr_in         = "/data/sicopolis/data/Simpson2009_GlacialMask//"
        file_in_grid    = trim(fldr_in)//"GISM20.cdf"
        file_in         = trim(fldr_in)//"glacmask.20.cdf"

        desc    = "Reconstructed Greenland ice sheet extent at the LGM"
        ref     = "Simpson, M. J. R., Milne, G. A., Huybrechts, P. and Long, A. J.: &
                  &Calibrating a glaciological model of the Greenland ice sheet &
                  &from the Last Glacial Maximum to present-day using field &
                  &observations of relative sea level and ice extent, &
                  &Quat. Sci. Rev., 28(17-18), 1631–1657, &
                  &doi:10.1016/j.quascirev.2009.03.004, 2009."

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_TOPO-LGM-S09.nc"

        ! Define the input data 
        nx = 83
        ny = 141 
        np = nx*ny 

        allocate(inp%lon(np),inp%lat(np),inp%var(np))

        ! Define the input points
        call nc_read(file_in_grid,"lon",inp%lon,start=[1,1],count=[nx,ny])
        call nc_read(file_in_grid,"lat",inp%lat,start=[1,1],count=[nx,ny])
        call points_init(points0,name="GISMp-20KM",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat) !,latlon=.TRUE.)

        ! Initialize mapping
        call map_init(map,points0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)     

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="km")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="km")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Read in current variable
        call nc_read(file_in,"zm",inp%var,start=[1,1],count=[nx,ny],missing_value=mv)
        
        ! Map variable to new grid
        call map_field(map,"mask",inp%var,outvar,outmask,"nn",fill=.TRUE.,missing_value=mv)

        ! Cut the mask if overlaps with border points 
!         where (grid%border) outvar = 0.d0 
        
        ! Write output variable to output file
        call nc_write(filename,"mask_lgm",int(outvar),dim1="xc",dim2="yc",missing_value=int(mv))

        ! Write variable metadata
        call nc_write_attr(filename,"mask_lgm","units","1")
        call nc_write_attr(filename,"mask_lgm","long_name", &
                    "Mask of ice sheet at extent at LGM (21 ka BP)")
        call nc_write_attr(filename,"mask_lgm","coordinates","lat2D lon2D")

        return 

    end subroutine LGMsimpson_to_grid

    subroutine huy3_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       Huy3 GLACIAL MASK DATA
        !       from Lecavalier et al. (2014)
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
            double precision, allocatable :: lon(:), lat(:)
            double precision, allocatable :: var(:)
        end type 

        type(inp_type)     :: inp
        integer            :: np 
        type(points_class) :: points0
        character(len=256) :: fldr_in, prefix, suffix 
        character(len=512) :: file_in 
        real(4)            :: times(40)
        character(len=4)   :: times_str(40), year_str

        type(map_class)  :: map 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: q 
        
        ! Define the input filenames
        !fldr_in         = "/data/sicopolis/data/Greenland/huy3_extent/"
        fldr_in         = "data/Greenland/huy3_extent/"
        prefix          = trim(fldr_in)//"ice_mask."
        suffix          = ".xyz"

        times = [21.0,20.0,19.0,18.0,17.0,16.5,16.0,15.5,15.0,14.5,14.0,13.5, &
                 13.0,12.5,12.0,11.5,11.0,10.5,10.0,9.5,9.0,8.5,8.0,7.5,7.0,6.5, &
                 6.0,5.5,5.0,4.5,4.0,3.5,3.0,2.5,2.0,1.5,1.0,0.5,0.1,0.0]
        times_str = ["21  ","20  ","19  ","18  ","17  ","16p5","16  ","15p5","15  ","14p5","14  ","13p5", &
                 "13  ","12p5","12  ","11p5","11  ","10p5","10  ","9p5 ","9   ","8p5 ","8   ","7p5 ","7   ","6p5 ", &
                 "6   ","5p5 ","5   ","4p5 ","4   ","3p5 ","3   ","2p5 ","2   ","1p5 ","1   ","0p5 ","0p1 ","0   "]

        desc    = "Reconstructed Greenland ice sheet extent during the last deglaciation"
        ref     = "Lecavalier, B. S., Milne, G. a., Simpson, M. J. R., Wake, L., Huybrechts, P., &
                  &Tarasov, L., Kjeldsen, K. K., Funder, S., Long, A. J., Woodroffe, S., &
                  &Dyke, A. S., et al.: A model of Greenland ice sheet deglaciation &
                  &constrained by observations of relative sea level and ice extent, &
                  &Quat. Sci. Rev., 102, 54–84, doi:10.1016/j.quascirev.2014.07.018, 2014."

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_TOPO-LGM-L14.nc"

        ! Define the input data 
        np = 11703

        allocate(inp%lon(np),inp%lat(np),inp%var(np))

        ! Define the input points
        file_in = trim(prefix)//"21"//trim(suffix)
        inp%lon = read_vector(file_in,n=np,col=1,skip=0)
        inp%lat = read_vector(file_in,n=np,col=2,skip=0)
        call points_init(points0,name="huy3",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat)

        ! Initialize mapping
        call map_init(map,points0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)     

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="km")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="km")
        call nc_write_dim(filename,"time", x=times,units="ka BP",unlimited=.TRUE.)
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        do q = 1, size(times) 

            ! Determine current filename 
            file_in = trim(prefix)//trim(times_str(q))//trim(suffix)

            ! Read in current variable
            inp%var = read_vector(file_in,n=np,col=3,skip=0)

            ! Map variable to new grid
            call map_field(map,"mask",inp%var,outvar,outmask,"nn",fill=.TRUE.,missing_value=mv)

            ! Write output variable to output file
            call nc_write(filename,"mask",int(outvar),dim1="xc",dim2="yc",dim3="time",missing_value=int(mv), &
                          start=[1,1,q],count=[grid%G%nx,grid%G%ny,1])

        end do 

        ! Write variable metadata
        call nc_write_attr(filename,"mask","units","1")
        call nc_write_attr(filename,"mask","long_name", &
                 "Mask of ice sheet at extent through the last deglaciation (21 ka BP to present)")
        call nc_write_attr(filename,"mask","coordinates","lat2D lon2D")

        return 

    end subroutine huy3_to_grid

    subroutine dated1_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       DATED-1 GLACIAL MASK DATA
        !       from Hughes et al. (2016)
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
            double precision, allocatable :: lon(:), lat(:)
            double precision, allocatable :: var(:)
        end type 

        type(inp_type)     :: inp
        integer            :: np 
        type(points_class) :: points0
        character(len=256) :: fldr_in, prefix, suffix 
        character(len=512) :: file_in 
        real(4)            :: times(40)
        character(len=4)   :: times_str(40), year_str

        type(map_class)  :: map 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: i, j, q 
        logical :: in_basin 

        ! Define the input filenames
!         fldr_in         = "/data/sicopolis/data/North/Hughes/"
!        fldr_in         = "/Users/robinson/wrk/mypapers/paleo_indexForcing_ruben/wrk/"
        fldr_in         = "/home/robinson/models/EURICE/gridding/data/"
        file_in         = trim(fldr_in)//"coord_hughes.txt"

        desc    = "Reconstructed Eurasian ice sheet extent during the last deglaciation"
        ref     = "Hughes, Anna L. C., Gyllencreutz, Richard, Lohne, Øystein S., Mangerud, Jan &
                  &Svendsen, John Inge: The last Eurasian ice sheets – a chronological database &
                  &and time-slice reconstruction, DATED-1, BOREAS, doi:10.1111/bor.12142, 2016."

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_TOPO-LGM-H16.nc"

        ! Define the input data 
        np = 13127

        allocate(inp%lon(np),inp%lat(np),inp%var(np))

        ! Define the input points
        inp%lon = read_vector(file_in,n=np,col=1,skip=1)
        inp%lat = read_vector(file_in,n=np,col=2,skip=1)

        ! Define input points for mapping
        call points_init(points0,grid,name="Hughes2016",x=inp%lon,y=inp%lat,latlon=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)     

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="km")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="km")
        call nc_write_dim(filename,"time", x=21.0,units="ka BP",unlimited=.TRUE.)
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Initially set all points to zero 
        outvar = 0.0 

        ! Loop over grid points and check point in polygon
        do j = 1, grid%G%ny 
        do i = 1, grid%G%nx 

            in_basin = point_in_polygon(real(grid%x(i,j)), real(grid%y(i,j)), &
                                        real(points0%x),   real(points0%y)) 
            
            ! If basin was found, save it
            if (in_basin) outvar(i,j) = 1.0 

        end do 
        end do 

        ! Write output variable to output file
        call nc_write(filename,"mask",int(outvar),dim1="xc",dim2="yc",missing_value=int(mv), &
                      start=[1,1],count=[grid%G%nx,grid%G%ny])

        ! Write variable metadata
        call nc_write_attr(filename,"mask","units","1")
        call nc_write_attr(filename,"mask","long_name", &
                 "Mask of ice sheet at extent at the last deglaciation (21 ka BP)")
        call nc_write_attr(filename,"mask","coordinates","lat2D lon2D")

        return 

    end subroutine dated1_to_grid

end module topo_reconstructions 

