module GreenlandVelocity 

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: grlvelj18_to_grid
    public :: grlvelj10_to_grid
    
contains 

    subroutine grlvelj18_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GREENLAND VELOCITY DATA (Joughin et al, 2018)
        !
        !       Data downloaded from here: http://dx.doi.org/10.5067/QUA5Q9SVMSJG
        !       Original format: .tif 250m x 250m resolution
        !       Aggregated to 2km x 2km by Ilaria Tabone in R using
        !       the included script `scripts/joughin2018.r 
        !       
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
        character(len=256) :: fldr_in, file_in

        type(map_class)  :: map
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: ux(:,:), uy(:,:)

        character(len=512) :: filename_H_ice 
        double precision, allocatable :: H_ice(:,:)
        
        ! Define the input filenames
!         fldr_in         = "/data/sicopolis/data/Greenland/"
        fldr_in         = "data/Greenland/Joughin2018_vel/"
        file_in         = trim(fldr_in)//"greenland_vel_mosaic250_v1_5km.txt"

        desc    = "Satellite based reconstruction of Greenland surface velocity"
        ref     = "Joughin, I., Smith, B. and Howat, I.: A complete map of Greenland &
                  &ice velocity derived from satellite data collected over 20 years. &
                  &Journal of Glaciology, 64(243), 1-11, doi:10.1017/jog.2017.73, 2018"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_VEL-J18.nc"

        write(filename_H_ice,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_TOPO-RTOPO-2.0.1.nc"

        ! Define the input data 
        !nx = 301
        !ny = 561
        np = 164346 !(5km x 5 km) nrow for each column of the .txt file
!         np = 1026480  ! (2km x 2km) 

        allocate(inp%lon(np),inp%lat(np),inp%var(np))

        ! Define the input points
        !call nc_read(file_in,"lon",inp%lon,start=[1,1,1],count=[nx,ny,1])
        !call nc_read(file_in,"lat",inp%lat,start=[1,1,1],count=[nx,ny,1])
        inp%lon = read_vector(file_in,n=np,col=3,skip=0)
        inp%lat = read_vector(file_in,n=np,col=4,skip=0)        
 
        call points_init(points0,name="Joughin2018-5km",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat)


        ! Initialize mapping
        call map_init(map,points0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)
        call grid_allocate(grid,H_ice)


        ! First load H_ice field for filtering where velocity should be zero 
        ! (assumes topography is already available on output domain)
        call nc_read(filename_H_ice,"H_ice",H_ice)

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="km")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="km")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)

        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! Read in current variable, map it, and write it
        !call nc_read(file_in,"surfvelx",inp%var,start=[1,1,1],count=[nx,ny,1],missing_value=mv)
        inp%var = read_vector(file_in,n=np,col=1,skip=0)   ! vx
        where(inp%var .eq. 0.0) inp%var = mv
        outvar = mv 
        call map_field(map,"vx",inp%var,outvar,outmask,"radius",fill=.FALSE.,missing_value=mv,radius=grid%G%dx)
        where(H_ice .eq. 0.0) outvar = mv 
        call nc_write(filename,"ux_srf",outvar,dim1="xc",dim2="yc",missing_value=mv)

        ! Write variable metadata
        call nc_write_attr(filename,"ux_srf","units","m/a")
        call nc_write_attr(filename,"ux_srf","long_name", &
                    "Surface velocity, x-component")
        call nc_write_attr(filename,"ux_srf","coordinates","lat2D lon2D")

        ! Read in current variable, map it, and write it
        !call nc_read(file_in,"surfvely",inp%var,start=[1,1,1],count=[nx,ny,1],missing_value=mv)
        inp%var = read_vector(file_in,n=np,col=2,skip=0)   ! vy
        where(inp%var .eq. 0.0) inp%var = mv
        outvar = mv 
        call map_field(map,"vx",inp%var,outvar,outmask,"radius",fill=.FALSE.,missing_value=mv,radius=grid%G%dx)
        where(H_ice .eq. 0.0) outvar = mv 
        call nc_write(filename,"uy_srf",outvar,dim1="xc",dim2="yc",missing_value=mv)

        ! Write variable metadata
        call nc_write_attr(filename,"uy_srf","units","m/a")
        call nc_write_attr(filename,"uy_srf","long_name", &
                    "Surface velocity, y-component")
        call nc_write_attr(filename,"uy_srf","coordinates","lat2D lon2D")

        call grid_allocate(grid,ux)
        call grid_allocate(grid,uy)

        call nc_read(filename,"ux_srf",ux)
        call nc_read(filename,"uy_srf",uy)

        outvar = mv
        where(ux .ne. mv .and. uy .ne. mv) outvar = sqrt(ux**2 + uy**2)

        call nc_write(filename,"uxy_srf",outvar,dim1="xc",dim2="yc",missing_value=mv)

        ! Write variable metadata
        call nc_write_attr(filename,"uxy_srf","units","m/a")
        call nc_write_attr(filename,"uxy_srf","long_name", &
                    "Surface velocity, magnitude")
        call nc_write_attr(filename,"uxy_srf","coordinates","lat2D lon2D")

        return

    end subroutine grlvelj18_to_grid

    subroutine grlvelj10_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GREENLAND VELOCITY DATA (Joughin et al, 2010)
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
        character(len=256) :: fldr_in, file_in

        type(map_class)  :: map 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: ux(:,:), uy(:,:)

        ! Define the input filenames
!         fldr_in         = "/data/sicopolis/data/Greenland/"
        fldr_in         = "data/Greenland/"
        file_in         = trim(fldr_in)//"Greenland_5km_v1.1.nc"

        desc    = "Satellite based reconstruction of Greenland surface velocity &
                  &(source data processed from Searise product Greenland_5km_v1.1.nc)"
        ref     = "Joughin, I., Smith, B. E., Howat, I. M., Scambos, T. and Moon, T.: &
                  &Greenland flow variability from ice-sheet-wide velocity mapping, &
                  &56(197), 415–430, 2010."

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_VEL-J10.nc"

        ! Define the input data 
        nx = 301
        ny = 561
        np = 168861  ! (301x561)

        allocate(inp%lon(np),inp%lat(np),inp%var(np))

        ! Define the input points
        call nc_read(file_in,"lon",inp%lon,start=[1,1,1],count=[nx,ny,1])
        call nc_read(file_in,"lat",inp%lat,start=[1,1,1],count=[nx,ny,1])
        call points_init(points0,name="Searise-5KM",mtype="latlon",units="degrees", &
                         lon180=.TRUE.,x=inp%lon,y=inp%lat)

        ! Initialize mapping
        call map_init(map,points0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

        ! Read in current variable, map it, and write it
        call nc_read(file_in,"surfvelx",inp%var,start=[1,1,1],count=[nx,ny,1],missing_value=mv)
        where(inp%var .eq. 0.0) inp%var = mv 
        call map_field(map,"surfvelx",inp%var,outvar,outmask,"radius",fill=.TRUE.,missing_value=mv,radius=grid%G%dx)
        call nc_write(filename,"Ux_srf",outvar,dim1="xc",dim2="yc",missing_value=mv)

        ! Write variable metadata
        call nc_write_attr(filename,"Ux_srf","units","m/a")
        call nc_write_attr(filename,"Ux_srf","long_name", &
                    "Surface velocity, x-component")
        call nc_write_attr(filename,"Ux_srf","coordinates","lat2D lon2D")

        call nc_read(file_in,"surfvely",inp%var,start=[1,1,1],count=[nx,ny,1],missing_value=mv)
        where(inp%var .eq. 0.0) inp%var = mv 
        call map_field(map,"surfvely",inp%var,outvar,outmask,"radius",fill=.TRUE.,missing_value=mv,radius=grid%G%dx)
        call nc_write(filename,"Uy_srf",outvar,dim1="xc",dim2="yc",missing_value=mv)

        ! Write variable metadata
        call nc_write_attr(filename,"Uy_srf","units","m/a")
        call nc_write_attr(filename,"Uy_srf","long_name", &
                    "Surface velocity, y-component")
        call nc_write_attr(filename,"Uy_srf","coordinates","lat2D lon2D")

        call grid_allocate(grid,ux)     
        call grid_allocate(grid,uy)     

        call nc_read(filename,"Ux_srf",ux)
        call nc_read(filename,"Uy_srf",uy)
        
        outvar = mv 
        where(ux .ne. mv .and. uy .ne. mv) outvar = sqrt(ux**2 + uy**2)

        call nc_write(filename,"Umag_srf",outvar,dim1="xc",dim2="yc",missing_value=mv)

        ! Write variable metadata
        call nc_write_attr(filename,"Umag_srf","units","m/a")
        call nc_write_attr(filename,"Umag_srf","long_name", &
                    "Surface velocity, magnitude")
        call nc_write_attr(filename,"Umag_srf","coordinates","lat2D lon2D")

        return 

    end subroutine grlvelj10_to_grid

end module GreenlandVelocity 

