module GeothermalHeatFlux 

    use gridding_datasets
    use coordinates
    use interp2D 
    use ncio 
    
    implicit none 

    private 
    public :: ghfMaule_to_grid
    
contains 

    subroutine ghfMaule_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GEOTHERMAL HEAT FLUX DATA - Maule et al., 2005, 2009
        !       http://core2.gsfc.nasa.gov/research/purucker/heatflux_updates.html
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
        character(len=256)   :: file_in, pname

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: i, np 

        
        ! Define input points from global data
        
        if (trim(domain) .eq. "Greenland") then 

            ! Define the input filename
!             file_in = "data/GeothermalHeatFlux/mct_mf7_foxmaule05_npolar.txt"
            file_in = "/data/sicopolis/data/GeothermalHeatFlux/mct_mf7_foxmaule05_npolar.txt"
            
            pname = "Maul05-npolar"
            np = 1165  ! Number of data points

        else if (trim(domain) .eq. "Antarctica") then

            ! Define the input filename
!             file_in = "data/GeothermalHeatFlux/heatflux_mf7_foxmaule05.txt"
            file_in = "/data/sicopolis/data/GeothermalHeatFlux/heatflux_mf7_foxmaule05.txt"
            
            pname = "Maule05"
            np = 656  ! Number of data points
            
        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        desc    = "Geothermal heat flux"
        ref     = "Maule, C. F., Purucker, M. E., Olsen, N., & Mosegaard, K.: &
                  &Heat flux anomalies in Antarctica revealed by satellite magnetic data, &
                  &Science, 309(5733), 464-467. \n &
                  &http://core2.gsfc.nasa.gov/research/purucker/heatflux_updates.html"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                          "_GHF-M05.nc"

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
        call points_init(pTOPO,name=trim(pname),mtype="latlon",units="degrees",x=inp%lon,y=inp%lat,lon180=.TRUE.)
        
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
        call map_field(map,"ghf",inp%var,outvar,outmask,"quadrant", &
                       fill=.TRUE.,missing_value=missing_value)

        write(*,*) "Range outvar: ",minval(outvar,outvar.ne.missing_value), &
                                    maxval(outvar,outvar.ne.missing_value)
        
        ! Fill any missing values (Antarctica)
        call fill_weighted(outvar,missing_value=missing_value)
        
        ! Write field to output file 
        call nc_write(filename,"ghf",real(outvar),dim1="xc",dim2="yc", &
                      missing_value=real(missing_value))

        ! Write variable metadata
        call nc_write_attr(filename,"ghf","units","mW m**-2")
        call nc_write_attr(filename,"ghf","long_name","Geothermal heat flux")
!         call nc_write_attr(filename,"ghf","grid_mapping",trim(grid%mtype))
        call nc_write_attr(filename,"ghf","coordinates","lat2D lon2D")
            
        return 

    end subroutine ghfMaule_to_grid

end module GeothermalHeatFlux 
