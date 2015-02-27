module sediments 

    use gridding_datasets
    use coordinates
    use interp2D 
    use ncio 
    
    implicit none 

    private 
    public :: sedLaske_to_grid
    
contains 

    subroutine sedLaske_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       SEDIMENT DATA - CRUST5.1
        !       http://igppweb.ucsd.edu/~gabi/sediment.html
        !
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type inpts_type 
            double precision, allocatable :: lon(:), lat(:), var(:)
        end type 

        type(inpts_type)     :: inp
        type(points_class)   :: pTOPO
        character(len=256)   :: file_in

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: i, np 

        
        ! Define input points from global data
        
        ! Define the input filename
!         file_in = "data/Sediments/CRUST5.1-L97/sedmap.gmt"
        file_in = "/data/sicopolis/data/Sediments/CRUST5.1-L97/sedmap.gmt"

        desc    = "Global sediment map updated from CRUST5.1"
        ref     = "Laske, G. and Masters, G.: &
                  &A Global Digital Map of Sediment Thickness, &
                  &EOS Trans. AGU, 78, F483, 1997. \n &
                  &http://igppweb.ucsd.edu/~gabi/sediment.html"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                          "_SED-L97.nc"

        np = 180*360  ! Number of data points (1x1 deg global grid)

        allocate(inp%lon(np),inp%lat(np),inp%var(np))

        ! File format: lon, lat, sed_thickness
        open(2,file=trim(file_in),status="old")
        do i = 1, np 
            read(2,*) inp%lon(i), inp%lat(i), inp%var(i) 
        end do 
        close(2)

        write(*,*) "lon: ",minval(inp%lon),maxval(inp%lon)
        write(*,*) "lat: ",minval(inp%lat),maxval(inp%lat)
        write(*,*) "var: ",minval(inp%var),maxval(inp%var)

        ! Define input points for mapping
        call points_init(pTOPO,name="L97-1DEG",mtype="latlon",units="degrees",x=inp%lon,y=inp%lat,lon180=.TRUE.)

        ! Initialize mapping
        call map_init(map,pTOPO,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

        ! ## MAP FIELD ##
        call map_field(map,"z_sed",inp%var,outvar,outmask,"quadrant", &
                       fill=.TRUE.,missing_value=missing_value)

        write(*,*) "Range outvar: ",minval(outvar), maxval(outvar)
        write(*,*) "Range outvar: ",minval(outvar,outvar.ne.missing_value), &
                                    maxval(outvar,outvar.ne.missing_value)
        
        ! Convert units [km => m]
        where(outvar .ne. missing_value) outvar = outvar *1d3 

        ! Write field to output file 
        call nc_write(filename,"z_sed",real(outvar),dim1="xc",dim2="yc", &
                      missing_value=real(missing_value))

        ! Write variable metadata
        call nc_write_attr(filename,"z_sed","units","m")
        call nc_write_attr(filename,"z_sed","long_name","Sediment thickness")
!         call nc_write_attr(filename,"z_sed","grid_mapping",trim(grid%mtype))
        call nc_write_attr(filename,"z_sed","coordinates","lat2D lon2D")
            
        return 

    end subroutine sedLaske_to_grid

end module sediments 
