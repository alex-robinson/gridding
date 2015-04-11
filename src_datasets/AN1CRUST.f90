module AN1CRUST

    use gridding_datasets
    use coordinates
    use interp2D 
    use ncio 
    
    implicit none 

    private 
    public :: An15litho_to_grid 


contains 

    subroutine An15litho_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       AN1-CRUST DATA
        !       http://www.seismolab.org/model/antarctica/lithosphere/index.html
        ! =========================================================
        
        implicit none 

        character(len=*)    :: domain, outfldr 
        type(grid_class)    :: grid 
        integer             :: max_neighbors 
        double precision    :: lat_lim 
        character(len=512)  :: filename 
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid0
        character(len=256) :: file_in
        double precision, allocatable :: lon(:), lat(:)
        double precision, allocatable :: invar(:,:) 

        type(map_class)  :: map 
        double precision, allocatable :: outvar(:,:), tmp(:,:)
        integer, allocatable          :: outmask(:,:)
        
        integer :: nx, ny
        double precision, parameter :: sigma = 30   ! kilometers

        ! Define the input filenames
        file_in = "/data/sicopolis/data/Antarctica/AN1-CRUST/AN1-CRUST.grd"
        desc    = "AN1-CRUST: A Moho depth map"
        ref     = "Meijian, A., Wiens, D., Yue, Z., Mei, F., Nyblade, A. A., &
                  &Kanao, M., Yuansheng, L., Maggi, A., Lévêque, J., &
                  &S-velocity Model and Inferred Moho Topography beneath the Antarctic Plate from Rayleigh Waves, &
                  &J. Geophys. Res., 120(1),359–383, doi:10.1002/2014JB011332, 2015. \n&
                  &http://www.seismolab.org/model/antarctica/lithosphere/index.html"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                          "_LITHO-AN15.nc"

        ! Define input grid
        nx = nc_size(file_in,"lon")
        ny = nc_size(file_in,"lat")

        allocate(lon(nx),lat(ny))
        call nc_read(file_in,"lon",lon)
        call nc_read(file_in,"lat",lat)

        call grid_init(grid0,name="AN1-CRUST",mtype="latlon",units="degrees",lon180=.TRUE.,x=lon,y=lat)

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)

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

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! ## PROCESS FIELDS ##

        call nc_read(file_in,"z",invar,missing_value=mv)
        where (abs(invar) .gt. 1d6) invar = mv 

        write(*,*) "invar: ", minval(invar), maxval(invar)

        call map_field(map,"z",invar,outvar,outmask,method="nng",sigma=sigma,fill=.TRUE.,missing_value=mv)
        where (abs(outvar) .gt. 1d6) outvar = mv 

        write(*,*) "outvar: ", minval(outvar), maxval(outvar)
        
        call fill_weighted(outvar,missing_value=mv)

        call nc_write(filename,"H_litho",real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
            
        ! Write variable metadata
        call nc_write_attr(filename,"H_litho","units","km")
        call nc_write_attr(filename,"H_litho","long_name","Lithospheric thickness")
        call nc_write_attr(filename,"H_litho","coordinates","lat2D lon2D")

        return 

    end subroutine An15litho_to_grid

end module AN1CRUST 
