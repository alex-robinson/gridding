module GeothermalHeatFlux 

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: ghfMartos_to_grid
    public :: ghfMaule_to_grid
    public :: ghfDavies_to_grid
    public :: ghfShapiro_to_grid

contains 

    subroutine ghfMartos_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GEOTHERMAL HEAT FLUX DATA - Martos et al., 2017 (Antarctica), 2018 (Greenland)
        !       https://doi.pangaea.de/10.1594/PANGAEA.882503
        !       https://doi.pangaea.de/10.1594/PANGAEA.892973
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
        type(points_class)   :: pts_in

        character(len=56)    :: pname, tmpstr(3)  
        character(len=512)   :: fldr_data 
        character(len=512)   :: files(4), varnames(4), units(4), long_names(4)
        character(len=512)   :: file_now, vnm_now, units_now, long_name_now

        type(map_class)  :: map 
        type(var_defs)   :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: i, q, np 

        ! Define input points from global data on lon/lat coordinates 

        varnames(1) = "ghf"
        varnames(2) = "ghf_err"
        varnames(3) = "depth_curie"
        varnames(4) = "depth_curie_err"
        
        units(1)    = "mW m**-2"
        units(2)    = "mW m**-2"
        units(3)    = "km"
        units(4)    = "km"
        
        long_names(1) = "Geothermal heat flux"
        long_names(2) = "Geothermal heat flux uncertainty"
        long_names(3) = "Curie depth"
        long_names(4) = "Curie depth uncertainty"
        
        if (trim(domain) .eq. "Greenland") then 

            fldr_data = "/data/sicopolis/data/GeothermalHeatFlux/Martos2018/"

            files(1) = trim(fldr_data)//"Martos2018_Greenland_ghf.txt"
            files(2) = trim(fldr_data)//"Martos2018_Greenland_ghf_err.txt"
            files(3) = trim(fldr_data)//"Martos2018_Greenland_depth_curie.txt"
            files(4) = trim(fldr_data)//"Martos2018_Greenland_depth_curie_err.txt"
            
            ! Name of set of points 
            pname = "Martos2018-GRL"

            ! Number of data points 
            np = 10716

            desc    = "Geothermal heat flux"
            ref     = "Martos, Yasmina M; Jordan, Tom A; Catalan, Manuel; Jordan, Thomas M; &
                      &Bamber, Jonathan L; Vaughan, David G (2018): Geothermal heat flux &
                      &reveals the Iceland hotspot track underneath Greenland. Geophysical &
                      & Research Letters, 45(16), 8214-8222, https://doi.org/10.1029/2018GL078289 &
                      &\n https://doi.pangaea.de/10.1594/PANGAEA.892973"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_GHF-M18.nc"

        else if (trim(domain) .eq. "Antarctica") then

            fldr_data = "/data/sicopolis/data/GeothermalHeatFlux/Martos2017/"

            files(1) = trim(fldr_data)//"Martos2017_Antarctica_ghf.txt"
            files(2) = trim(fldr_data)//"Martos2017_Antarctica_ghf_err.txt"
            files(3) = trim(fldr_data)//"Martos2017_Antarctica_depth_curie.txt"
            files(4) = trim(fldr_data)//"Martos2017_Antarctica_depth_curie_err.txt"
            
            ! Name of set of points 
            pname = "Martos2017-ANT"
            
            ! Number of data points 
            np = 59329
            
            desc    = "Geothermal heat flux"
            ref     = "Martos, Yasmina M; Catalan, Manuel; Jordan, Tom A; Golynsky, Alexander V; &
                      &Golynsky, Dmitry A; Eagles, Graeme; Vaughan, David G (2017): Heat flux &
                      &distribution of Antarctica unveiled. Geophysical Research Letters, &
                      &44(22), 11417-11426, https://doi.org/10.1002/2017GL075609 &
                      &\n https://doi.org/10.1594/PANGAEA.882503"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_GHF-M17.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 

        end if 

        ! Allocate input data points 
        allocate(inp%lon(np),inp%lat(np),inp%var(np))

        ! Load data from first file to get lon/lat coordinates of points 
        ! File format: lon, lat, ghf with header
        open(2,file=trim(files(1)),status="old")
        read(2,*) tmpstr(1), tmpstr(2), tmpstr(3) 
        do i = 1, np 
            read(2,*) inp%lon(i), inp%lat(i), inp%var(i) 
        end do 
        close(2)

        write(*,*) "lon: ",minval(inp%lon),maxval(inp%lon)
        write(*,*) "lat: ",minval(inp%lat),maxval(inp%lat)
        write(*,*) "var: ",minval(inp%var),maxval(inp%var)

        ! Define input points for mapping
        call points_init(pts_in,name=trim(pname),mtype="latlon",units="degrees",x=inp%lon,y=inp%lat,lon180=.TRUE.)
        
        ! Initialize mapping
        call map_init(map,pts_in,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

        ! Loop over variables and map each one to our grid 
        do q = 1, size(files) 

            file_now      = files(q)
            vnm_now       = varnames(q)
            units_now     = units(q) 
            long_name_now = long_names(q) 

            ! Read in current data 
            ! File format: lon, lat, var 
            open(2,file=trim(file_now),status="old")
            read(2,*) tmpstr(1), tmpstr(2), tmpstr(3) 
            do i = 1, np 
                read(2,*) inp%lon(i), inp%lat(i), inp%var(i) 
            end do 
            close(2)

            ! ## MAP FIELD ##
            call map_field(map,vnm_now,inp%var,outvar,outmask,"nng",fill=.TRUE.,sigma=10.d0,missing_value=mv)

            write(*,*) "Range invar:  ",minval(inp%var,inp%var.ne.mv), maxval(inp%var,inp%var.ne.mv)
            write(*,*) "Range outvar: ",minval(outvar,outvar.ne.mv), maxval(outvar,outvar.ne.mv)
            
            ! Fill any missing values
            call fill_weighted(outvar,missing_value=missing_value)
        
            ! Write field to output file 
            call nc_write(filename,vnm_now,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))

            ! Write variable metadata
            call nc_write_attr(filename,vnm_now,"units",trim(units_now))
            call nc_write_attr(filename,vnm_now,"long_name",trim(long_name_now))
            call nc_write_attr(filename,vnm_now,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine ghfMartos_to_grid

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
        call map_field(map,"ghf",inp%var,outvar,outmask,"nng",fill=.TRUE.,sigma=40.d0,missing_value=mv)

        write(*,*) "Range outvar: ",minval(outvar,outvar.ne.mv), &
                                    maxval(outvar,outvar.ne.mv)
        
        ! Fill any missing values (Antarctica)
        call fill_weighted(outvar,missing_value=missing_value)
        
        ! Write field to output file 
        call nc_write(filename,"ghf",real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))

        ! Write variable metadata
        call nc_write_attr(filename,"ghf","units","mW m**-2")
        call nc_write_attr(filename,"ghf","long_name","Geothermal heat flux")
        call nc_write_attr(filename,"ghf","coordinates","lat2D lon2D")
            
        return 

    end subroutine ghfMaule_to_grid

    subroutine ghfDavies_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GEOTHERMAL HEAT FLUX DATA - Davies, 2013
        !       Davies, J. H.: Global map of solid Earth surface heat flow,
        !       Geochem. Geophys. Geosys., 14, 10, 
        !       doi:10.1002/ggge.20271, 2013
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
            double precision, allocatable :: lon(:), lat(:)
            double precision, allocatable :: mean(:), med(:), err(:)
        end type 

        type(inpts_type)     :: inp
        type(points_class)   :: pTOPO
        character(len=256)   :: file_in, pname, tmp 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: i, np 

        
        ! Define input points from global data
        
        ! Define the input filename
        file_in = "/data/sicopolis/data/GeothermalHeatFlux/Davies2013/Davies2013_ghf_sup-0001.txt"
        
        pname = "Davies2013"
        np = 10312  ! Number of data points

        desc    = "Geothermal heat flux"
        ref     = "Davies, J. H.: Global map of solid Earth surface heat flow, &
                  &Geochem. Geophys. Geosys., 14, 10, doi:10.1002/ggge.20271, 2013"

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                          "_GHF-D13.nc"

        allocate(inp%lon(np),inp%lat(np))
        allocate(inp%mean(np),inp%med(np),inp%err(np))

        ! File format: lon, lat, sed_thickness
        open(2,file=trim(file_in),status="old")
        read(2,*) tmp ! 1-line header
        do i = 1, np 
            read(2,*) inp%lon(i), inp%lat(i), inp%mean(i), inp%med(i), inp%err(i) 
        end do 
        close(2)

        write(*,*) "lon: ",minval(inp%lon),maxval(inp%lon)
        write(*,*) "lat: ",minval(inp%lat),maxval(inp%lat)
        write(*,*) "mean: ",minval(inp%mean),maxval(inp%mean)
        write(*,*) "med: ",minval(inp%med),maxval(inp%med)
        write(*,*) "err: ",minval(inp%err),maxval(inp%err)

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

        ! ## MAP FIELDS ##

        ! ## ghf_mean ##
        call map_field(map,"ghf_mean",inp%mean,outvar,outmask,"nng", &
                       fill=.TRUE.,sigma=80.d0,missing_value=missing_value)

        ! Fill any missing values (Antarctica)
        call fill_weighted(outvar,missing_value=missing_value)
        
        ! Write field to output file 
        call nc_write(filename,"ghf_mean",real(outvar),dim1="xc",dim2="yc", &
                      missing_value=real(missing_value))

        ! Write variable metadata
        call nc_write_attr(filename,"ghf_mean","units","mW m**-2")
        call nc_write_attr(filename,"ghf_mean","long_name","Geothermal heat flux (mean)")
        call nc_write_attr(filename,"ghf_mean","coordinates","lat2D lon2D")
        
        ! ## ghf_median ##
        call map_field(map,"ghf_med",inp%med,outvar,outmask,"nng", &
                       fill=.TRUE.,sigma=80.d0,missing_value=missing_value)

        ! Fill any missing values (Antarctica)
        call fill_weighted(outvar,missing_value=missing_value)
        
        ! Write field to output file 
        call nc_write(filename,"ghf_med",real(outvar),dim1="xc",dim2="yc", &
                      missing_value=real(missing_value))

        ! Write variable metadata
        call nc_write_attr(filename,"ghf_med","units","mW m**-2")
        call nc_write_attr(filename,"ghf_med","long_name","Geothermal heat flux (median)")
        call nc_write_attr(filename,"ghf_med","coordinates","lat2D lon2D")
        
        ! ## ghf_err ##
        call map_field(map,"ghf_err",inp%err,outvar,outmask,"nng", &
                       fill=.TRUE.,sigma=80.d0,missing_value=missing_value)

        ! Fill any missing values (Antarctica)
        call fill_weighted(outvar,missing_value=missing_value)
        
        ! Write field to output file 
        call nc_write(filename,"ghf_err",real(outvar),dim1="xc",dim2="yc", &
                      missing_value=real(missing_value))

        ! Write variable metadata
        call nc_write_attr(filename,"ghf_err","units","mW m**-2")
        call nc_write_attr(filename,"ghf_err","long_name","Geothermal heat flux (error)")
        call nc_write_attr(filename,"ghf_err","coordinates","lat2D lon2D")
            
        return 

    end subroutine ghfDavies_to_grid

    subroutine ghfShapiro_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       GEOTHERMAL HEAT FLUX DATA - Shapiro and Ritzwoller, 2004
        !       Shapiro, N. M. and Ritzwoller, M. H.: Inferring surface heat
        !       flux distributions guided by a global seismic model: 
        !       particular application to Antarctica, Earth Planet. Sci. Lett., 
        !       223(1-2), 213–224, doi:10.1016/j.epsl.2004.04.011, 2004.
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
            double precision, allocatable :: lon(:), lat(:)
            double precision, allocatable :: mean(:), med(:), err(:)
        end type 

        type(inpts_type)     :: inp
        type(points_class)   :: pTOPO
        character(len=256)   :: file_in, pname, tmp 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: i, np 

        
        ! Define input points from global data
        
        ! Define the input filename
        file_in = "/data/sicopolis/data/GeothermalHeatFlux/Shapiro2004/hfmap.asc"
        
        pname = "Shapiro2004"
        np = 65160  ! Number of data points

        desc    = "Geothermal heat flux"
        ref     = "Shapiro, N. M. and Ritzwoller, M. H.: Inferring surface heat &
                  &flux distributions guided by a global seismic model: &
                  &particular application to Antarctica, Earth Planet. Sci. Lett., &
                  &223(1-2), 213–224, doi:10.1016/j.epsl.2004.04.011, 2004."

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                          "_GHF-S04.nc"

        allocate(inp%lon(np),inp%lat(np))
        allocate(inp%mean(np),inp%med(np),inp%err(np))

        ! File format: lon, lat, ghf, ghf_sigma
        open(2,file=trim(file_in),status="old")
        do i = 1, np 
            read(2,*) inp%lon(i), inp%lat(i), inp%mean(i), inp%err(i) 
        end do 
        close(2)

        write(*,*) "lon: ",minval(inp%lon),maxval(inp%lon)
        write(*,*) "lat: ",minval(inp%lat),maxval(inp%lat)
        write(*,*) "mean: ",minval(inp%mean),maxval(inp%mean)
        write(*,*) "err: ",minval(inp%err),maxval(inp%err)

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

        ! ## MAP FIELDS ##

        ! ## ghf_mean ##
        call map_field(map,"ghf",inp%mean,outvar,outmask,"nng", &
                       fill=.TRUE.,sigma=40.d0,missing_value=missing_value)

        ! Fill any missing values (Antarctica)
        call fill_weighted(outvar,missing_value=missing_value)
        
        ! Write field to output file 
        call nc_write(filename,"ghf",real(outvar),dim1="xc",dim2="yc", &
                      missing_value=real(missing_value))

        ! Write variable metadata
        call nc_write_attr(filename,"ghf","units","mW m**-2")
        call nc_write_attr(filename,"ghf","long_name","Geothermal heat flux (mean)")
        call nc_write_attr(filename,"ghf","coordinates","lat2D lon2D")
        
        ! ## ghf_median ##
        call map_field(map,"ghf_sigma",inp%err,outvar,outmask,"nng", &
                       fill=.TRUE.,sigma=40.d0,missing_value=missing_value)

        ! Fill any missing values (Antarctica)
        call fill_weighted(outvar,missing_value=missing_value)
        
        ! Write field to output file 
        call nc_write(filename,"ghf_sigma",real(outvar),dim1="xc",dim2="yc", &
                      missing_value=real(missing_value))

        ! Write variable metadata
        call nc_write_attr(filename,"ghf_sigma","units","mW m**-2")
        call nc_write_attr(filename,"ghf_sigma","long_name","Geothermal heat flux (standard deviation)")
        call nc_write_attr(filename,"ghf_sigma","coordinates","lat2D lon2D")
            
        return 

    end subroutine ghfShapiro_to_grid

end module GeothermalHeatFlux 
