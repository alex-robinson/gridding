module topographies_grl

    use gridding_datasets
    use coordinates
    use interp2D 
    use polygons 
    use ncio 
    
    implicit none 

    private 
    public :: Bamber13_to_grid
    public :: Morlighem14_to_grid 

contains 

    subroutine Bamber13_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,grad_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       TOPO DATA
        !
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim, grad_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid0
        character(len=256) :: file_in
        type(var_defs), allocatable :: vars(:)
        double precision, allocatable :: invar(:,:) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: q, k, m, i, l, n_var 
        integer :: thin_by = 5 
        character(len=128) :: method, grad_lim_str  

        logical, allocatable :: mask_reg(:,:)
        real(4), allocatable :: xp(:), yp(:) 

        grad_lim_str = "" 
        if (grad_lim .gt. 0.09d0) then 
            write(grad_lim_str,"(a,f3.1)") "_gl", grad_lim 
        else if (grad_lim .gt. 0.d0) then 
            write(grad_lim_str,"(a,f4.2)") "_gl", grad_lim 
        end if 

        ! Define input grid
        if (trim(domain) .eq. "Greenland") then 
            
            ! Define topography (Bamber et al. 2013) grid and input variable field
            select case(thin_by)
                case(10)
                    call grid_init(grid0,name="TOPO-B13-10KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-1300.d0,dx=10.d0,nx=251,y0=-3500.d0,dy=10.d0,ny=301, &
                            lambda=-39.d0,phi=90.d0,alpha=7.5d0)

                case(5)
                    call grid_init(grid0,name="TOPO-B13-5KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-1300.d0,dx=5.d0,nx=501,y0=-3500.d0,dy=5.d0,ny=601, &
                            lambda=-39.d0,phi=90.d0,alpha=7.5d0)

                case DEFAULT
                    write(*,*) "Bamber13_to_grid:: error: thin_by can only be 5 or 10."
                    stop 

            end select 

            ! Define the input filenames
            file_in = "/data/sicopolis/data/Greenland/Greenland_bedrock_topography_V3.nc"
            desc    = "Greenland bedrock and surface topography (V3)"
            ref     = "Bamber, J. L., Griggs, J. A., Hurkmans, R. T. W. L., &
                      &Dowdeswell, J. A., Gogineni, S. P., Howat, I., Mouginot, J., &
                      &Paden, J., Palmer, S., Rignot, E., and Steinhage, D.: &
                      &A new bed elevation dataset for Greenland, &
                      &The Cryosphere, 7, 499-510, doi:10.5194/tc-7-499-2013, 2013."

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-B13"//trim(grad_lim_str)//".nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(vars(4))
        call def_var_info(vars(1),trim(file_in),"BedrockElevation","zb",units="m",long_name="Bedrock elevation",method="nng")
        call def_var_info(vars(2),trim(file_in),"SurfaceElevation","zs",units="m",long_name="Surface elevation",method="nng")
        call def_var_info(vars(3),trim(file_in),"IceThickness",     "H",units="m",long_name="Ice thickness",method="nng")
        call def_var_info(vars(4),trim(file_in),"LandMask",      "mask",units="(0 - 4)", &
                          long_name="Land mask",method="nn")

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp(2501,3001))

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ## FIELDS ##
        do i = 1, size(vars)
            var_now = vars(i) 
            call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=mv)
            call thin(invar,tmp,by=thin_by)

            if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
                where( invar .eq. mv ) invar = 0.d0 
            end if
            if (trim(var_now%nm_out) .eq. "zb") then 
                call fill_mean(invar,missing_value=mv,fill_value=-1500.d0)
            end if 

            method = "nng"
            if (trim(var_now%nm_out) .eq. "mask") method = "nn" 

            call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
                           sigma=grid%G%dx*0.5d0,fill=.TRUE.,missing_value=mv)
            
            if (var_now%method .eq. "nn") then 
                call fill_nearest(outvar,missing_value=mv)
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",missing_value=int(mv))
            else
                call fill_weighted(outvar,missing_value=mv)
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
            end if 
            
            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        ! Modify variables for consistency and gradient limit 

        ! Allocate helper arrays
        call grid_allocate(grid,zs)
        call grid_allocate(grid,zb)
        call grid_allocate(grid,H)

        ! Re-load data
        call nc_read(filename,"zs",zs)
        call nc_read(filename,"zb",zb)
        call nc_read(filename,"H",H)
        
        ! Eliminate problematic regions for this domain ========
        call grid_allocate(grid,mask_reg)    
        
        ! Baffin Bay
        allocate(xp(4),yp(4))
        xp = [-63.5,-57.7,-53.9,-57.7]
        yp = [ 69.6, 67.3, 63.3, 55.0]
        mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
        where (mask_reg .and. zb .gt. -600.d0) zb = mv 

        ! Iceland 
        xp = [-17.0,-23.8,-31.2,-22.1]
        yp = [ 69.1, 68.6, 64.1, 63.2]
        mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
        where (mask_reg .and. zb .gt. -200.d0) zb = mv 
        where (mask_reg .and. zb .gt. -200.d0) zs = mv 

        ! Land in the Northeast corner
        xp = [  5.1, -4.5, 22.0, 22.1]
        yp = [ 78.8, 83.5, 80.1, 80.2]
        mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
        where (mask_reg .and. zb .gt. -200.d0) zb = mv 
        where (mask_reg .and. zb .gt. -200.d0) zs = mv 

        ! Replaces problematic regions with regional mean values or zero for surface
        call fill_weighted(zb,missing_value=mv)
        call fill_weighted(zs,missing_value=mv,fill_value=0.d0)

        ! ======================================================


        ! Apply gradient limit as needed
        if (grad_lim .gt. 0.d0) then 
            ! Limit the gradient (m/m) to below threshold 
            call limit_gradient(zs,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            call limit_gradient(zb,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            
        end if 

        ! First clean up the field based on original H
        where (H .lt. 1.d0) 
            H  = 0.d0 
            zs = zb 
        end where 
        
        ! Also make sure zs is always higher than zb 
        where (zs .lt. zb) zs = zb 
        
        ! Now make sure zs is zero over the ocean
        where (zs .lt. 0.d0) zs = 0.d0 

        ! Adjust H again for consistency
        H = 0.d0 
        where (zs .ne. 0.d0) H = zs-zb 

        ! Re-write fields 
        call nc_write(filename,"zs",real(zs),dim1="xc",dim2="yc",missing_value=real(mv))
        call nc_write(filename,"zb",real(zb),dim1="xc",dim2="yc",missing_value=real(mv))
        call nc_write(filename,"H", real(H), dim1="xc",dim2="yc",missing_value=real(mv))

        ! Define new masks ==========

        ! ocean-land-ice-shelf (0,1,2,3) mask 
        outmask = 0     ! Ocean
        where (zs .gt. 0.d0) outmask = 1    ! Land
        where ( H .gt. 0.d0) outmask = 2    ! Grounded ice
        where (zs .gt. 0.d0 .and. zs-zb .gt. H) outmask = 3   ! Floating ice 

        call nc_write(filename,"mask", outmask, dim1="xc",dim2="yc",missing_value=int(mv), &
                      long_name="Mask (ocean=0,land=1,grounded-ice=2,floating-ice=3)")

        ! Region masks

        ! Greenland
        outmask = 1 

        ! Ellesmere Island
        if (allocated(xp)) deallocate(xp)
        if (allocated(yp)) deallocate(yp)
        allocate(xp(7),yp(7)) 
        xp = [-59.3,-60.0,-70.8,-78.0,-99.5,-96.7,-88.7 ]
        yp = [ 85.0, 82.4, 79.7, 76.0, 81.0, 81.6, 85.0 ]
        mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
        where (mask_reg) outmask = 2 

        call nc_write(filename,"mask_reg", outmask, dim1="xc",dim2="yc",missing_value=int(mv), &
                      long_name="Region mask (Greenland=1,Ellesmere Island=2)")

        ! ===========================


        return 

    end subroutine Bamber13_to_grid

    subroutine Morlighem14_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       TOPO DATA
        !       http://sites.uci.edu/morlighem/dataproducts/mass-conservation-dataset/
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid0
        character(len=512) :: file_in
        type(var_defs), allocatable :: vars(:)
        double precision, allocatable :: invar(:,:) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: q, k, m, i, l, n_var 
        integer :: thin_by = 10 
        character(len=128) :: method 

        integer :: status, ncid 

        ! Define input grid
        if (trim(domain) .eq. "Greenland") then 
            
            ! Define grid and input variable field
            select case(thin_by)
                case(10)  ! 150m => 1.5km 
                    call grid_init(grid0,name="ESPG-3413-1.5KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-637.925d0,dx=1.5d0,nx=1001,y0=-3349.425d0,dy=1.5d0,ny=1794, &
                            lambda=-45.d0,phi=70.d0,alpha=7.2d0)

                case DEFAULT
                    write(*,*) "Morlighem14_to_grid:: error: thin_by can only be 10."
                    stop 

            end select 

            ! Define the input filenames
            file_in = "/data/sicopolis/data/Greenland/Morlighem2014_topo/MCdataset-2014-11-19_NetCDF3.nc"
            desc    = "BedMachine: Greenland dataset based on mass conservation, 19-Nov-2014 (v1.7)"
            ref     = "Morlighem, M., Rignot, E., Mouginot, J., Seroussi, H. and Larour, E., &
                      &Deeply incised submarine glacial valleys beneath the Greenland Ice Sheet, &
                      &Nat. Geosci., 7, 418-422, doi:10.1038/ngeo2167, 2014. \n&
                      &http://sites.uci.edu/morlighem/dataproducts/mass-conservation-dataset"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-M14.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(vars(6))
        call def_var_info(vars(1),trim(file_in),"bed",      "zb",units="m",long_name="Bedrock elevation")
        call def_var_info(vars(2),trim(file_in),"surface",  "zs",units="m",long_name="Surface elevation")
        call def_var_info(vars(3),trim(file_in),"thickness","H",units="m",long_name="Ice thickness")
        call def_var_info(vars(4),trim(file_in),"errbed",   "zb_err",units="m",long_name="Bedrock / ice thickness error")
        call def_var_info(vars(5),trim(file_in),"mask",     "mask",units="(0 - 3)", &
                          long_name="(0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice)",method="nn")
        call def_var_info(vars(6),trim(file_in),"source",     "mask_source",units="(0 - 3)", &
                          long_name="data source (0 = none, 1 = gimpdem, 2 = Mass conservation, &
                                    &4 = interpolation, 5 = hydrostatic equilibrium, 6=kriging)",method="nn")

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp(10018,17946))

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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

        ! ## FIELDS ##
        do i = 1, size(vars)
            var_now = vars(i) 
            call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=mv)
            call thin(invar,tmp,by=thin_by)

            if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
                where( invar .eq. mv ) invar = 0.d0 
            end if
            if (trim(var_now%nm_out) .eq. "zb") then 
                call fill_mean(invar,missing_value=mv,fill_value=-1500.d0)
            end if 

            method = "radius"
            if (trim(var_now%nm_out) .eq. "mask")        method = "nn" 
            if (trim(var_now%nm_out) .eq. "mask_source") method = "nn" 

            call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
                           radius=grid%G%dx*grid%xy_conv*0.75d0,fill=.TRUE.,missing_value=mv)
            
            if (var_now%method .eq. "nn") then 
                call fill_nearest(outvar,missing_value=mv)
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",missing_value=int(mv))
            else
                call fill_weighted(outvar,missing_value=mv)
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
            end if 
            
            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine Morlighem14_to_grid

end module topographies_grl 
