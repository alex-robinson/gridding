module GreenlandVelocity 

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: grlvelj17_to_grid
    public :: grlvelj10_to_grid
    
contains 

    subroutine grlvelj17_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GREENLAND VELOCITY DATA (Joughin et al, 2017)
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
        fldr_in         = "data/Joughin2017_vel/"
        file_in         = trim(fldr_in)//"greenland_vel_mosaic250_vx_v1.txt"

        desc    = "Satellite based reconstruction of Greenland surface velocity &
                  &(source data processed from Searise product Greenland_5km_v1.1.nc)cHANGE!!!" 
        ref     = "Joughin, I., Smith, B. E., Howat, I. M., Scambos, T. and Moon, T.: &
                  &Greenland flow variability from ice-sheet-wide velocity mapping, &
                  &56(197), 415–430, 2017. (CHANGE!!)"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_VEL-J17.nc"

        ! Define the input data 
        !nx = 301
        !ny = 561
        np = 6032 !(5km x 5 km) nrow for each column of the .txt file

        allocate(inp%lon(np),inp%lat(np),inp%var(np))

        ! Define the input points
        !call nc_read(file_in,"lon",inp%lon,start=[1,1,1],count=[nx,ny,1])
        !call nc_read(file_in,"lat",inp%lat,start=[1,1,1],count=[nx,ny,1])
        inp%lon = read_vector(file_in,n=np,col=3,skip=0)
        inp%lat = read_vector(file_in,n=np,col=4,skip=0)        
 
        call points_init(points0,name="Joughin2017-5km",mtype="latlon",units="degrees", &
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
        !call nc_read(file_in,"surfvelx",inp%var,start=[1,1,1],count=[nx,ny,1],missing_value=mv)
        inp%var = read_vector(file_in,n=np,col=1,skip=0)   ! vx
        where(inp%var .eq. 0.0) inp%var = mv
        call map_field(map,"vx",inp%var,outvar,outmask,"radius",fill=.TRUE.,missing_value=mv,radius=grid%G%dx)
        call nc_write(filename,"Ux_srf",outvar,dim1="xc",dim2="yc",missing_value=mv)

        ! Write variable metadata
        call nc_write_attr(filename,"Ux_srf","units","m/a")
        call nc_write_attr(filename,"Ux_srf","long_name", &
                    "Surface velocity, x-component")
        call nc_write_attr(filename,"Ux_srf","coordinates","lat2D lon2D")

        ! Read in current variable, map it, and write it
        !call nc_read(file_in,"surfvely",inp%var,start=[1,1,1],count=[nx,ny,1],missing_value=mv)
        inp%var = read_vector(file_in,n=np,col=2,skip=0)   ! vy
        where(inp%var .eq. 0.0) inp%var = mv
        call map_field(map,"vx",inp%var,outvar,outmask,"radius",fill=.TRUE.,missing_value=mv,radius=grid%G%dx)
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

    end subroutine grlvelj17_to_grid

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

