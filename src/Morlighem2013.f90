module Morlighem2013

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: morlighem2013_taub_to_grid

contains

    subroutine morlighem2013_taub_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       Morlighem et al (2013) Antarctica inversion data 
        !         - Basal stress 
        !         - 3D velocity field 
        !         - f_shr 
        !
        ! =========================================================
        
        implicit none 

        character(len=*), intent(IN) :: outfldr 
        type(grid_class), intent(IN) :: grid 
        character(len=*), intent(IN) :: domain 
        integer,          intent(IN) :: max_neighbors 
        double precision, intent(IN) :: lat_lim


        ! Local variables 
        character(len=512) :: filename, file_input  
        character(len=1024) :: desc, ref 

        type(points_class) :: pts0 

        type points_vector_type 
            double precision, allocatable :: lon(:), lat(:)
            double precision, allocatable :: taub(:), taud(:), f_vbvs(:)
        end type 

        type(points_vector_type) :: inp
        integer :: i, np 

        type(map_class)  :: map 

        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        if (trim(domain) .ne. "Antarctica") then 
            write(*,*) "Error: Morlighem 2013 dataset is only valid for domain=Antarctica."
            stop 
        end if 

        ! Define the input data information
        file_input = "/data/sicopolis/data/Antarctica/Morlighem2013/Morlighem2013_data.txt"
        desc = "Antarctic basal stress inversion."
        ref  = "Morlighem et al.: Inversion of basal friction in Antarctica using exact and &
               &incomplete adjoints of a higher-order model, JGR: Earth Surface, 118, 1746–1753, &
               &doi:10.1002/jgrf.20125, 2013"

        ! Number of data points 
        np = 32933

        allocate(inp%lon(np),inp%lat(np),inp%taub(np),inp%taud(np),inp%f_vbvs(np))

        ! File format: lat, lon, basin 
        open(2,file=trim(file_input),status="old")
        do i = 1, np 
            read(2,*) inp%lat(i), inp%lon(i), inp%taub(i), inp%taud(i), inp%f_vbvs(i) 
        end do 
        close(2)

        write(*,*) "lon:    ", minval(inp%lon),    maxval(inp%lon)
        write(*,*) "lat:    ", minval(inp%lat),    maxval(inp%lat)
        write(*,*) "taub:   ", minval(inp%taub),   maxval(inp%taub)
        write(*,*) "taud:   ", minval(inp%taud),   maxval(inp%taud)
        write(*,*) "f_vbvs: ", minval(inp%f_vbvs), maxval(inp%f_vbvs)

        ! Define input points for mapping
        call points_init(pts0,name="M13-taub",mtype="latlon",units="degrees",x=inp%lon,y=inp%lat,lon180=.TRUE.)
        
        ! Initialize mapping
        call map_init(map,pts0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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


        return 

    end subroutine morlighem2013_taub_to_grid

end module Morlighem2013
