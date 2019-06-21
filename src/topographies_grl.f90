module topographies_grl

    use gridding_datasets
    use coord
    use ncio 
    
    use gaussian_filter 

    use regions 

    implicit none 

    private 
    public :: Morlighem17_to_grid 
    public :: Morlighem14_to_grid 
    public :: Bamber13_to_grid !, Bamber13_to_grid_conserv

contains 

        subroutine Morlighem17_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,grad_lim,thin_by)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       TOPO DATA
        !       http://sites.uci.edu/morlighem/dataproducts/bedmachine-greenland/
        !       https://daacdata.apps.nsidc.org/pub/DATASETS/IDBMG4_BedMachineGr/
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim
        double precision :: grad_lim   
        integer, intent(IN) :: thin_by 

        ! Local variables 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid0
        character(len=512) :: file_in
        type(var_defs), allocatable :: vars(:)
        double precision, allocatable :: invar(:,:) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp(:,:), tmp_rev(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: q, k, m, i, l, n_var, j 
        character(len=128) :: method, grad_lim_str  
        integer :: status, ncid 

        double precision, allocatable :: var_fill(:,:)
        double precision, allocatable :: zs(:,:), zb(:,:), H(:,:)
        character(len=512) :: filename0 
        double precision :: sigma 

        grad_lim_str = "" 
        if (grad_lim .gt. 0.09d0) then 
            write(grad_lim_str,"(a,f3.1)") "_gl", grad_lim 
        else if (grad_lim .gt. 0.d0) then 
            write(grad_lim_str,"(a,f4.2)") "_gl", grad_lim 
        end if 
        
        ! Define input grid
        if (trim(domain) .eq. "Greenland") then 
            
            ! Define grid and input variable field
            select case(thin_by)
                case(15)  ! 150m => 2.25km 
                    call grid_init(grid0,name="ESPG-3413-2.25KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-652.925d0,dx=2.25d0,nx=682,y0=-3384.425d0,dy=2.25d0,ny=1224, &
                            lambda=-45.d0,phi=70.d0)
                case(10)  ! 150m => 1.5km 
                    call grid_init(grid0,name="ESPG-3413-1.5KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-652.925d0,dx=1.5d0,nx=1022,y0=-3384.425d0,dy=1.5d0,ny=1835, &
                            lambda=-45.d0,phi=70.d0)

                case DEFAULT
                    call grid_init(grid0,name="ESPG-3413-150M",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-652.925d0,dx=0.150d0,nx=10218,y0=-3384.425d0,dy=0.150d0,ny=18346, &
                            lambda=-45.d0,phi=70.d0)
                    
            end select 

            ! Define the input filenames
            file_in = "/data/sicopolis/data/Greenland/BedMachineGreenland-2017-09-20.nc"
            desc    = "BedMachine v3 (2017-09-20): Greenland dataset based on mass conservation"
            ref     = "Morlighem M. et al., (2017), BedMachine v3: Complete bed topography and &
                      &ocean bathymetry mapping of Greenland from multi-beam echo sounding combined &
                      &with mass conservation, Geophys. Res. Lett., 44, doi:10.1002/2017GL074954. &
                      &(http://onlinelibrary.wiley.com/doi/10.1002/2017GL074954/full)."

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-M17"//trim(grad_lim_str)//".nc"

            ! Define filename holding RTOPO-2.0.1 data
            write(filename0,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-RTOPO-2.0.1"//".nc"
        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(vars(6))
        call def_var_info(vars(1),trim(file_in),"bed",      "z_bed",units="m",long_name="Bedrock elevation")
        call def_var_info(vars(2),trim(file_in),"surface",  "z_srf",units="m",long_name="Surface elevation")
        call def_var_info(vars(3),trim(file_in),"thickness","H_ice",units="m",long_name="Ice thickness")
        call def_var_info(vars(4),trim(file_in),"errbed",   "z_bed_err",units="m",long_name="Bedrock / ice thickness error")
        call def_var_info(vars(5),trim(file_in),"mask",     "mask",units="(0 - 3)", &
                          long_name="(0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice)",method="nn")
        call def_var_info(vars(6),trim(file_in),"source",     "mask_source",units="(0 - 3)", &
                          long_name="data source (0 = none, 1 = gimpdem, 2 = Mass conservation, &
                                    &4 = interpolation, 5 = hydrostatic equilibrium, 6=kriging)",method="nn")

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)
        
        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp_rev(10218,18346))
        allocate(tmp(10218,18346))

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

            call nc_read(trim(var_now%filename),var_now%nm_in,tmp_rev,missing_value=mv)
            do j = 1, size(tmp_rev,2)
                tmp(:,j) = tmp_rev(:,size(tmp_rev,2)-j+1)
            end do 
            call thin(invar,tmp,by=thin_by,missing_value=mv)

            if (trim(var_now%nm_out) .eq. "H_ice" .or. trim(var_now%nm_out) .eq. "z_srf") then 
                where( invar .eq. mv ) invar = 0.d0 
            end if

if (.FALSE.) then 

            method = "nn"
            if (trim(var_now%nm_out) .eq. "mask")        method = "nn" 
            if (trim(var_now%nm_out) .eq. "mask_source") method = "nn" 

            ! Perform Gaussian smoothing at high resolution 
            if (trim(var_now%nm_out) .eq. "H_ice" .or. &
                trim(var_now%nm_out) .eq. "z_srf" .or. &
                trim(var_now%nm_out) .eq. "z_bed") then 
                
                sigma = grid%G%dx / 2.0 

                call filter_gaussian(var=invar,sigma=sigma,dx=grid0%G%dx)

            end if 

            outvar = mv 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
                           radius=grid%G%dx, &
                           sigma=grid%G%dx*0.5d0,fill=.FALSE.,missing_value=mv)

else 
            ! Perform conservative interpolation 

            method = "mean"
            if (trim(var_now%nm_out) .eq. "mask")        method = "count" 
            if (trim(var_now%nm_out) .eq. "mask_source") method = "count" 

            call map_field_conservative_map1(map%map,var_now%nm_in,invar,outvar, &
                                                            method=method,missing_value=mv)

end if 

            if (trim(var_now%nm_out) .eq. "z_srf") then
                write(*,"(a,3f10.2)") "maxval(z_srf): ", maxval(outvar), maxval(invar), maxval(tmp)
            end if 

            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",missing_value=int(mv))
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
            end if 
            
            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 


        ! === Update specific variables from other information ===

        call grid_allocate(grid,zs)
        call grid_allocate(grid,zb)
        call grid_allocate(grid,H)
        call grid_allocate(grid,var_fill)
        
        ! Fill missing data with rtopo2 

        ! rtopo2 bedrock
        call nc_read(filename0,"z_bed",var_fill)
        call nc_read(filename, "z_bed",zb,missing_value=mv)
        where(zb .eq. mv) zb = var_fill 
        var_now = vars(1)
        call nc_write(filename,var_now%nm_out,real(zb),dim1="xc",dim2="yc",missing_value=real(mv))
        
        ! rtopo2 surface
        call nc_read(filename0,"z_srf",var_fill)
        call nc_read(filename, "z_srf",zs,missing_value=mv)
        where(zs .eq. mv) zs = var_fill 
        var_now = vars(2)
        call nc_write(filename,var_now%nm_out,real(zs),dim1="xc",dim2="yc",missing_value=real(mv))
    
        ! rtopo2 ice thickness 
        call nc_read(filename0,"H_ice",var_fill)
        call nc_read(filename, "H_ice",H,missing_value=mv)
        where(H .eq. mv) H = var_fill 
        var_now = vars(3)
        call nc_write(filename,var_now%nm_out,real(H),dim1="xc",dim2="yc",missing_value=real(mv))
                

        ! Modify variables for consistency

        ! Re-load data
        call nc_read(filename,"z_srf",zs)
        call nc_read(filename,"z_bed",zb)
        call nc_read(filename,"H_ice",H)
        
        ! Apply gradient limit as needed
        if (grad_lim .gt. 0.d0) then 
            ! Limit the gradient (m/m) to below threshold 
            call limit_gradient(zs,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            call limit_gradient(zb,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            
        end if 

        ! Eliminate problematic regions for this domain ========
        call clean_greenland(zs,zb,grid)

        call clean_thickness(zs,zb,H)

        ! Re-write fields 
        call nc_write(filename,"z_srf",real(zs),dim1="xc",dim2="yc",missing_value=real(mv))
        call nc_write(filename,"z_bed",real(zb),dim1="xc",dim2="yc",missing_value=real(mv))
        call nc_write(filename,"H_ice",real(H), dim1="xc",dim2="yc",missing_value=real(mv))

        ! Define new masks ==========

        ! ocean-land-ice-shelf (0,1,2,3) mask 
        outmask = 0     ! Ocean
        where (zs .gt. 0.d0) outmask = 1    ! Land
        where ( H .gt. 0.d0) outmask = 2    ! Grounded ice
        where (zs .gt. 0.d0 .and. zs-zb .gt. H) outmask = 3   ! Floating ice 

        call nc_write(filename,"mask", outmask, dim1="xc",dim2="yc",missing_value=int(mv), &
                      long_name="Mask (ocean=0,land=1,grounded-ice=2,floating-ice=3)")

        return 

    end subroutine Morlighem17_to_grid

        subroutine Morlighem14_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,grad_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       TOPO DATA
        !       http://sites.uci.edu/morlighem/dataproducts/mass-conservation-dataset/
        !       ftp://sidads.colorado.edu/DATASETS/IDBMG4_BedMachineGr/
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim, grad_lim  
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid0
        character(len=512) :: file_in
        type(var_defs), allocatable :: vars(:)
        double precision, allocatable :: invar(:,:) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp(:,:), tmp_rev(:,:)
        integer, allocatable          :: outmask(:,:)
        integer :: q, k, m, i, l, n_var, j 
        integer :: thin_by = 10 
        character(len=128) :: method, grad_lim_str  
        integer :: status, ncid 

        double precision, allocatable :: var_fill(:,:)
        double precision, allocatable :: zs(:,:), zb(:,:), H(:,:)
        character(len=512) :: filename0 

        grad_lim_str = "" 
        if (grad_lim .gt. 0.09d0) then 
            write(grad_lim_str,"(a,f3.1)") "_gl", grad_lim 
        else if (grad_lim .gt. 0.d0) then 
            write(grad_lim_str,"(a,f4.2)") "_gl", grad_lim 
        end if 

        ! Define input grid
        if (trim(domain) .eq. "Greenland") then 
            
            ! Define grid and input variable field
            select case(thin_by)
                case(10)  ! 150m => 1.5km 
                    call grid_init(grid0,name="ESPG-3413-1.5KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-637.175d0,dx=1.5d0,nx=1001,y0=-3348.675d0,dy=1.5d0,ny=1794, &
                            lambda=-45.d0,phi=70.d0,alpha=20.0d0)

                case DEFAULT
                    write(*,*) "Morlighem14_to_grid:: error: thin_by can only be 10."
                    stop 

            end select 

!             ! Default grid at 150 m resolution
!             call grid_init(grid0,name="ESPG-3413-150M",mtype="polar_stereographic", &
!                             units="kilometers",lon180=.TRUE., &
!                             x0=-637.925d0,dx=0.15d0,nx=10018,y0=-3349.425d0,dy=0.15d0,ny=17946, &
!                             lambda=-45.d0,phi=70.d0,alpha=20.0d0)

            ! Define the input filenames
            file_in = "/data/sicopolis/data/Greenland/Morlighem2014_topo/MCdataset-2015-04-27.nc"
            desc    = "BedMachine: Greenland dataset based on mass conservation, 2015-04-27 (v2.0)"
            ref     = "Morlighem, M., Rignot, E., Mouginot, J., Seroussi, H. and Larour, E., &
                      &Deeply incised submarine glacial valleys beneath the Greenland Ice Sheet, &
                      &Nat. Geosci., 7, 418-422, doi:10.1038/ngeo2167, 2014. \n&
                      &http://sites.uci.edu/morlighem/dataproducts/mass-conservation-dataset"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-M14"//trim(grad_lim_str)//".nc"

            ! Define filename holding ETOPO1 data
            write(filename0,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-ETOPO1"//trim(grad_lim_str)//".nc"
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
        allocate(tmp_rev(10018,17946))
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

            method = "radius"
            if (trim(var_now%nm_out) .eq. "mask")        method = "nn" 
            if (trim(var_now%nm_out) .eq. "mask_source") method = "nn" 

            call nc_read(trim(var_now%filename),var_now%nm_in,tmp_rev,missing_value=mv)
            do j = 1, size(tmp_rev,2)
                tmp(:,j) = tmp_rev(:,size(tmp_rev,2)-j+1)
            end do 
            if (var_now%method .eq. "nn") then 
                call thin(invar,tmp,by=thin_by,missing_value=mv)
            else 
                call thin(invar,tmp,by=thin_by,missing_value=mv)
!                 call thin_ave(invar,tmp,by=thin_by,missing_value=mv)  ! Diffuses zs too much!!
            end if 
            if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
                where( invar .eq. mv ) invar = 0.d0 
            end if
            
            outvar = mv 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
                           radius=grid%G%dx, &
                           sigma=grid%G%dx*0.5d0,fill=.FALSE.,missing_value=mv)
            
            if (trim(var_now%nm_out) .eq. "zs") then
                write(*,"(a,3f10.2)") "maxval(zs): ", maxval(outvar), maxval(invar), maxval(tmp)
            end if 

            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",missing_value=int(mv))
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
            end if 
            
            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 


        ! === Update specific variables from other information ===

        call grid_allocate(grid,zs)
        call grid_allocate(grid,zb)
        call grid_allocate(grid,H)
        call grid_allocate(grid,var_fill)
        
        ! Bedrock from ETOPO-1 
        ! Also load etopo bedrock, use it to replace high latitude regions 
        call nc_read(filename0,"zb",var_fill)
        call nc_read(filename, "zb",zb,missing_value=mv)
        where(zb .eq. mv) zb = var_fill 
        var_now = vars(1)
        call nc_write(filename,var_now%nm_out,real(zb),dim1="xc",dim2="yc",missing_value=real(mv))
        

        ! Modify variables for consistency and gradient limit 

        ! Re-load data
        call nc_read(filename,"zs",zs)
        call nc_read(filename,"zb",zb)
        call nc_read(filename,"H",H)
        
        ! Eliminate problematic regions for this domain ========
        call clean_greenland(zs,zb,grid)

        ! Apply gradient limit as needed
        if (grad_lim .gt. 0.d0) then 
            ! Limit the gradient (m/m) to below threshold 
            call limit_gradient(zs,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            call limit_gradient(zb,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            
        end if 

        call clean_thickness(zs,zb,H)

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

        return 

    end subroutine Morlighem14_to_grid

!     subroutine Bamber13_to_grid_conserv(outfldr,grid)
!         ! Convert the variables to the desired grid format and write to file
!         ! =========================================================
!         !
!         !       TOPO DATA
!         !
!         ! =========================================================
        
!         implicit none 

!         character(len=*), intent(IN) :: outfldr 
!         type(grid_class), intent(IN) :: grid 

!         ! Local variables
!         character(len=512) :: filename 
!         character(len=1024) :: desc, ref 

!         type(grid_class)   :: grid00, grid0
!         character(len=256) :: file_in
!         type(var_defs), allocatable :: vars(:)
!         double precision, allocatable :: invar(:,:) 

!         type(map_class)  :: map 
!         type(var_defs) :: var_now 
!         double precision, allocatable :: outvar(:,:), tmp(:,:)
!         integer, allocatable          :: outmask(:,:)
!         double precision, allocatable :: var_fill(:,:)
!         double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
!         double precision, allocatable :: zb_neg(:,:), zs_sl(:,:)

!         integer :: q, k, m, i, l, n_var  
!         character(len=128) :: method 
!         character(len=512) :: filename0 

!         ! Original grid (independent projection)
!         call grid_init(grid00,name="TOPO-B13-1KM",mtype="polar_stereographic", &
!                     units="kilometers",lon180=.TRUE., &
!                     x0=-1300.d0,dx=1.d0,nx=2501,y0=-3500.d0,dy=1.d0,ny=3001, &
!                     lambda=-39.d0,phi=71.d0,alpha=19.0d0)

!         ! Intermediate (very) high resolution grid (target projection)
!         call grid_init(grid0,name="TOPO-B13-1KM",mtype="polar_stereographic", &
!                     units="kilometers",lon180=.TRUE., &
!                     x0=-1300.d0,dx=1.d0,nx=2501,y0=-3500.d0,dy=1.d0,ny=3001, &
!                     lambda=-39.d0,phi=71.d0,alpha=19.0d0)


!             ! Define the input filenames
!             file_in = "/data/sicopolis/data/Greenland/Greenland_bedrock_topography_V3.nc"
!             desc    = "Greenland bedrock and surface topography (V3)"
!             ref     = "Bamber, J. L., Griggs, J. A., Hurkmans, R. T. W. L., &
!                       &Dowdeswell, J. A., Gogineni, S. P., Howat, I., Mouginot, J., &
!                       &Paden, J., Palmer, S., Rignot, E., and Steinhage, D.: &
!                       &A new bed elevation dataset for Greenland, &
!                       &The Cryosphere, 7, 499-510, doi:10.5194/tc-7-499-2013, 2013."

!             ! Define the output filename 
!             write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
!                               "_TOPO-B13"//trim(grad_lim_str)//".nc"

!             ! Define filename holding ETOPO1 data
!             write(filename0,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
!                               "_TOPO-ETOPO1"//trim(grad_lim_str)//".nc"
!         else

!             write(*,*) "Domain not recognized: ",trim(domain)
!             stop 
!         end if 

!         ! Define the variables to be mapped 
!         allocate(vars(4))
!         call def_var_info(vars(1),trim(file_in),"BedrockElevation","zb",units="m",long_name="Bedrock elevation",method="nn")
!         call def_var_info(vars(2),trim(file_in),"SurfaceElevation","zs",units="m",long_name="Surface elevation",method="nn")
!         call def_var_info(vars(3),trim(file_in),"IceThickness",     "H",units="m",long_name="Ice thickness",method="nn")
!         call def_var_info(vars(4),trim(file_in),"LandMask",      "mask",units="(0 - 4)", &
!                           long_name="Land mask",method="nn")

!         ! Allocate the input grid variable
!         call grid_allocate(grid0,invar)

!         ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
!         allocate(tmp(2501,3001))

!         ! Initialize mapping
!         call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

!         ! Initialize output variable arrays
!         call grid_allocate(grid,outvar)
!         call grid_allocate(grid,outmask)    
        
!         ! Initialize the output file
!         call nc_create(filename)
!         call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
!         call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
!         call nc_write_dim(filename,"month",x=[1,2,3,4,5,6,7,8,9,10,11,12],units="month")
!         call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
!         ! Write meta data 
!         call nc_write_attr(filename,"Description",desc)
!         call nc_write_attr(filename,"Reference",ref)

!         ! ## FIELDS ##
!         do i = 1, size(vars)
!             var_now = vars(i) 

!             method = "radius"
!             if (trim(var_now%nm_out) .eq. "mask") method = "nn" 

!             call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=mv)

!             if (thin_by .gt. 1) then 
!                 if (var_now%method .eq. "nn") then 
!                     call thin(invar,tmp,by=thin_by,missing_value=mv)
!                 else 
!                     call thin_ave(invar,tmp,by=thin_by,missing_value=mv)
!                 end if 
!             else 
!                 invar = tmp 
!             end if 

!             if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
!                 where( invar .eq. mv ) invar = 0.d0 
!             end if

!             outvar = mv 
!             call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
!                            radius=grid%G%dx,sigma=grid%G%dx*0.5d0,fill=.FALSE.,missing_value=mv)
            
!             if (var_now%method .eq. "nn") then 
!                 call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",missing_value=int(mv))
!             else
!                 call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
!             end if 
            
!             ! Write variable metadata
!             call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
!             call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
!             call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
!         end do 

!         ! Interpolate only below sea-level points for fjord filling
!         call grid_allocate(grid,zb_neg)
!         call grid_allocate(grid,zs_sl)

!         if (.FALSE.) then 

!             var_now = vars(1)   ! BedrockElevation 
!             method = "radius"

!             call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=mv)

!             if (thin_by .gt. 1) then 
!                 call thin_ave(invar,tmp,by=thin_by,missing_value=mv) 
!             else 
!                 invar = tmp 
!             end if 
    
!             ! Set land points to missing (to avoid making fjords overly shallow)
!             where( invar .gt. 0.d0 ) invar = mv 

!             call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
!                            radius=grid%G%dx,sigma=grid%G%dx*0.5d0,fill=.TRUE.,missing_value=mv)
            
!             zb_neg = outvar 

!             var_now = vars(2)   ! SurfaceElevation 
!             method = "nn"

!             call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=mv)

!             if (thin_by .gt. 1) then 
!                 call thin_ave(invar,tmp,by=thin_by,missing_value=mv) 
!             else 
!                 invar = tmp 
!             end if 
    
!             ! Set land points to missing (to avoid making fjords overly shallow)
!             where( invar .gt. 0.d0 ) invar = mv 

!             call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
!                            radius=grid%G%dx,sigma=grid%G%dx*0.5d0,fill=.TRUE.,missing_value=mv)
            
!             zs_sl = outvar 

!         else 
!             zb_neg = mv
!             zs_sl  = mv 

!         end if 

!         ! Modify variables for consistency and gradient limit 

!         ! Allocate helper arrays
!         call grid_allocate(grid,zs)
!         call grid_allocate(grid,zb)
!         call grid_allocate(grid,H)
!         call grid_allocate(grid,var_fill)
        
!         ! Bedrock from ETOPO-1
!         ! Also load etopo bedrock, use it to replace high latitude regions
!         call nc_read(filename0,"zb",var_fill)
!         call nc_read(filename, "zb",zb,missing_value=mv)
!         call nc_read(filename, "zs",zs)  
!         where(zb .eq. mv) zb = var_fill 
!         var_now = vars(1)
!         call nc_write(filename,var_now%nm_out,real(zb),dim1="xc",dim2="yc",missing_value=real(mv))
        

!         ! Re-load data
!         call nc_read(filename,"zs",zs)
!         call nc_read(filename,"zb",zb)
!         call nc_read(filename,"H",H)
        
!         ! Fill in fjords from second zb and zs fields
!         where (zb_neg .lt. 0.d0 .and. zb_neg .ne. mv) zb = zb_neg 
!         where (zs_sl  .ne. mv) zs = zs_sl 
        
!         ! Eliminate problematic regions for this domain ========
!         call clean_greenland(zs,zb,grid)

!         ! Apply gradient limit as needed
!         if (grad_lim .gt. 0.d0) then 
!             ! Limit the gradient (m/m) to below threshold 
!             call limit_gradient(zs,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
!             call limit_gradient(zb,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            
!         end if 

!         call clean_thickness(zs,zb,H)

!         ! Re-write fields 
!         call nc_write(filename,"zs",real(zs),dim1="xc",dim2="yc",missing_value=real(mv))
!         call nc_write(filename,"zb",real(zb),dim1="xc",dim2="yc",missing_value=real(mv))
!         call nc_write(filename,"H", real(H), dim1="xc",dim2="yc",missing_value=real(mv))

!         ! Define new masks ==========

!         ! ocean-land-ice-shelf (0,1,2,3) mask 
!         outmask = 0     ! Ocean
!         where (zs .gt. 0.d0) outmask = 1    ! Land
!         where ( H .gt. 0.d0) outmask = 2    ! Grounded ice
!         where (zs .gt. 0.d0 .and. zs-zb .gt. H) outmask = 3   ! Floating ice 

!         call nc_write(filename,"mask", outmask, dim1="xc",dim2="yc",missing_value=int(mv), &
!                       long_name="Mask (ocean=0,land=1,grounded-ice=2,floating-ice=3)")

!         return 

!     end subroutine Bamber13_to_grid_conserv

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
        double precision, allocatable :: var_fill(:,:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        double precision, allocatable :: zb_neg(:,:), zs_sl(:,:)

        integer :: q, k, m, i, l, n_var 
        integer :: thin_by = 1 
        character(len=128) :: method, grad_lim_str  
        character(len=512) :: filename0 

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
                case(20)
                    call grid_init(grid0,name="TOPO-B13-20KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-1290.d0,dx=20.d0,nx=126,y0=-3490.d0,dy=20.d0,ny=151, &
                            lambda=-39.d0,phi=71.d0,alpha=19.0d0)

                case(10)
                    call grid_init(grid0,name="TOPO-B13-10KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-1295.d0,dx=10.d0,nx=251,y0=-3495.d0,dy=10.d0,ny=301, &
                            lambda=-39.d0,phi=71.d0,alpha=19.0d0)

                case(5)
                    call grid_init(grid0,name="TOPO-B13-5KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-1297.5d0,dx=5.d0,nx=501,y0=-3497.5d0,dy=5.d0,ny=601, &
                            lambda=-39.d0,phi=71.d0,alpha=19.0d0)

                case(1)
                    ! Original Bamber grid (1KM)
                    call grid_init(grid0,name="TOPO-B13-1KM",mtype="polar_stereographic", &
                                units="kilometers",lon180=.TRUE., &
                                x0=-1300.d0,dx=1.d0,nx=2501,y0=-3500.d0,dy=1.d0,ny=3001, &
                                lambda=-39.d0,phi=71.d0,alpha=19.0d0)

                case DEFAULT
                    write(*,*) "Bamber13_to_grid:: error: thin_by can only be 10, 5 or 1."
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

            ! Define filename holding ETOPO1 data
            write(filename0,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-ETOPO1"//trim(grad_lim_str)//".nc"
        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(vars(4))
        call def_var_info(vars(1),trim(file_in),"BedrockElevation","zb",units="m",long_name="Bedrock elevation",method="nn")
        call def_var_info(vars(2),trim(file_in),"SurfaceElevation","zs",units="m",long_name="Surface elevation",method="nn")
        call def_var_info(vars(3),trim(file_in),"IceThickness",     "H",units="m",long_name="Ice thickness",method="nn")
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

            method = "radius"
            if (trim(var_now%nm_out) .eq. "mask") method = "nn" 

            call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=mv)

            if (thin_by .gt. 1) then 
                if (var_now%method .eq. "nn") then 
                    call thin(invar,tmp,by=thin_by,missing_value=mv)
                else 
                    call thin_ave(invar,tmp,by=thin_by,missing_value=mv)
                end if 
            else 
                invar = tmp 
            end if 

            if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
                where( invar .eq. mv ) invar = 0.d0 
            end if

            outvar = mv 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
                           radius=grid%G%dx,sigma=grid%G%dx*0.5d0,fill=.FALSE.,missing_value=mv)
            
            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",missing_value=int(mv))
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
            end if 
            
            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        ! Interpolate only below sea-level points for fjord filling
        call grid_allocate(grid,zb_neg)
        call grid_allocate(grid,zs_sl)

        if (.FALSE.) then 

            var_now = vars(1)   ! BedrockElevation 
            method = "radius"

            call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=mv)

            if (thin_by .gt. 1) then 
                call thin_ave(invar,tmp,by=thin_by,missing_value=mv) 
            else 
                invar = tmp 
            end if 
    
            ! Set land points to missing (to avoid making fjords overly shallow)
            where( invar .gt. 0.d0 ) invar = mv 

            call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
                           radius=grid%G%dx,sigma=grid%G%dx*0.5d0,fill=.TRUE.,missing_value=mv)
            
            zb_neg = outvar 

            var_now = vars(2)   ! SurfaceElevation 
            method = "nn"

            call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=mv)

            if (thin_by .gt. 1) then 
                call thin_ave(invar,tmp,by=thin_by,missing_value=mv) 
            else 
                invar = tmp 
            end if 
    
            ! Set land points to missing (to avoid making fjords overly shallow)
            where( invar .gt. 0.d0 ) invar = mv 

            call map_field(map,var_now%nm_in,invar,outvar,outmask,method, &
                           radius=grid%G%dx,sigma=grid%G%dx*0.5d0,fill=.TRUE.,missing_value=mv)
            
            zs_sl = outvar 

        else 
            zb_neg = mv
            zs_sl  = mv 

        end if 

        ! Modify variables for consistency and gradient limit 

        ! Allocate helper arrays
        call grid_allocate(grid,zs)
        call grid_allocate(grid,zb)
        call grid_allocate(grid,H)
        call grid_allocate(grid,var_fill)
        
        ! Bedrock from ETOPO-1
        ! Also load etopo bedrock, use it to replace high latitude regions
        call nc_read(filename0,"zb",var_fill)
        call nc_read(filename, "zb",zb,missing_value=mv)
        call nc_read(filename, "zs",zs)  
        where(zb .eq. mv) zb = var_fill 
        var_now = vars(1)
        call nc_write(filename,var_now%nm_out,real(zb),dim1="xc",dim2="yc",missing_value=real(mv))
        

        ! Re-load data
        call nc_read(filename,"zs",zs)
        call nc_read(filename,"zb",zb)
        call nc_read(filename,"H",H)
        
        ! Fill in fjords from second zb and zs fields
        where (zb_neg .lt. 0.d0 .and. zb_neg .ne. mv) zb = zb_neg 
        where (zs_sl  .ne. mv) zs = zs_sl 
        
        ! Eliminate problematic regions for this domain ========
        call clean_greenland(zs,zb,grid)

        ! Apply gradient limit as needed
        if (grad_lim .gt. 0.d0) then 
            ! Limit the gradient (m/m) to below threshold 
            call limit_gradient(zs,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            call limit_gradient(zb,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            
        end if 

        call clean_thickness(zs,zb,H)

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

        return 

    end subroutine Bamber13_to_grid

    subroutine clean_thickness(zs,zb,H)

        implicit none 

        double precision, intent(INOUT) :: zs(:,:), zb(:,:), H(:,:)
        
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

        return 

    end subroutine clean_thickness

    subroutine clean_greenland(zs,zb,grid)

        implicit none 

        double precision, intent(INOUT) :: zs(:,:), zb(:,:)
        type(grid_class), intent(IN)    :: grid 

        real(4), allocatable :: xp(:), yp(:) 
        logical, allocatable :: in_reg(:,:)

        call grid_allocate(grid,in_reg)    

        ! Baffin Bay
        if (allocated(xp)) deallocate(xp)
        if (allocated(yp)) deallocate(yp)
        allocate(xp(4),yp(4))
        xp = [-63.5,-57.7,-53.9,-57.7]
        yp = [ 69.6, 67.3, 63.3, 55.0]
        in_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
        where (in_reg .and. zb .gt. -600.d0) zb = mv 

!         ! Iceland 
!         if (allocated(xp)) deallocate(xp)
!         if (allocated(yp)) deallocate(yp)
!         allocate(xp(4),yp(4))
!         xp = [-17.0,-23.8,-31.2,-22.1]
!         yp = [ 69.1, 68.6, 64.1, 63.2]
!         in_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
!         where (in_reg .and. zb .gt. -200.d0) zb = mv 
!         where (in_reg .and. zb .gt. -200.d0) zs = mv 

!         ! Svalbard
!         if (allocated(xp)) deallocate(xp)
!         if (allocated(yp)) deallocate(yp)
!         allocate(xp(5),yp(5))
!         xp = [ 40.0, 40.0, 20.0,  0.0, -10.0 ]
!         yp = [ 85.0, 80.0, 73.0, 75.0,  85.0 ]
!         in_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
!         where (in_reg .and. zb .gt. -200.d0) zb = mv 
!         where (in_reg .and. zb .gt. -200.d0) zs = mv 
    
        ! Note: Iceland and Svalbard can be masked out via the regions mask,
        !       so there is no need to modify the topography here. 

        ! Replaces problematic regions with regional mean values or zero for surface
        call fill_weighted(zb,missing_value=mv)
        call fill_weighted(zs,missing_value=mv,fill_value=0.d0)

        return 

    end subroutine clean_greenland

end module topographies_grl 
