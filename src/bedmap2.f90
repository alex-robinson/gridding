module bedmap2 

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: bedmap2_to_grid
    public :: bedmap2vel_to_grid
    public :: bedmap2acc_to_grid

contains 

    subroutine bedmap2_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
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
        character(len=512) :: filename, infldr, prefix  
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid0
        character(len=256) :: file_invariant
        type(var_defs), allocatable :: invariant(:) 
        double precision, allocatable :: invar(:,:) 
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp1(:,:), tmp2(:,:), tmp3(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 
        character(len=128) :: grad_lim_str  

        double precision, allocatable :: zb0(:,:)
        character(len=512) :: filename0

        logical, allocatable :: mask_reg(:,:)
        real(4), allocatable :: xp(:), yp(:) 

        grad_lim_str = "" 
        if (grad_lim .gt. 0.09d0) then 
            write(grad_lim_str,"(a,f3.1)") "_gl", grad_lim 
        else if (grad_lim .gt. 0.d0) then 
            write(grad_lim_str,"(a,f4.2)") "_gl", grad_lim 
        end if 


        ! Define input grid
        if (trim(domain) .eq. "Antarctica") then 
            
            ! Define topography (BEDMAP2) grid and input variable field
            call grid_init(grid0,name="BEDMAP2-10KM",mtype="polar_stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-3328.d0,dx=10.d0,nx=666,y0=-3328.d0,dy=10.d0,ny=666,lambda=0.d0,phi=-71.d0)

            ! Original x0,y0 values
!             x0=-3333.d0, y0=-3333.d0 

            ! Define the input filenames
            infldr         = "output/Antarctica/"
            file_invariant = trim(infldr)//"ANT-10KM_TOPO-BEDMAP2.nc"
            desc    = "Antarctica bedrock and surface topography (BEDMAP2)"
            ref     = "Fretwell et al.: Bedmap2: improved ice bed, surface and &
                      &thickness datasets for Antarctica, The Cryosphere, 7, 375-393, &
                      &doi:10.5194/tc-7-375-2013, 2013."

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-BEDMAP2"//trim(grad_lim_str)//".nc"

            ! Define filename holding ETOPO1 data
            write(filename0,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-ETOPO1"//trim(grad_lim_str)//".nc"
        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(4))
        call def_var_info(invariant(1),file_invariant,"zs","zs",units="m",long_name="Surface elevation",method="con")
        call def_var_info(invariant(2),file_invariant,"zb","zb",units="m",long_name="Bedrock elevation",method="con")
        call def_var_info(invariant(3),file_invariant,"H","H",units="m",long_name="Ice thickness",method="con")
        call def_var_info(invariant(4),file_invariant,"mask_ice","mask",units="(0 - 1)",long_name="Ice mask",method="nn")

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)
        
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

        ! ## INVARIANT FIELDS ##
        do i = 1, size(invariant)
            var_now = invariant(i) 
            call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=mv)
            call thin(invar,tmp1,by=10,missing_value=mv)        ! Diffuses too much - check?
!             call thin_ave(invar,tmp1,by=10,missing_value=mv)
            if (trim(var_now%nm_out) .eq. "H" .or. &
                trim(var_now%nm_out) .eq. "zs") then 
                where( invar .eq. mv ) invar = 0.d0 
            end if

            if (trim(var_now%nm_out) .eq. "zb") then 
!                 call fill_mean(invar,missing_value=mv,fill_value=-1001.d0)
!                 call fill_mean(invar,missing_value=mv)
                call fill_nearest(invar,missing_value=mv)
            end if 

            if (trim(var_now%nm_out) .eq. "mask_ice") then 
                where ( invar .eq. 1.d0 ) invar = 3.d0 
                where ( invar .eq. 0.d0 ) invar = 2.d0 
                where ( invar .eq. missing_value ) invar = 0.d0 
            end if 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,radius=50.d0, &
                          fill=.TRUE.,sigma=grid%G%dx*0.5,missing_value=mv)
            call fill_mean(outvar,missing_value=mv)
            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",missing_value=nint(mv))
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
            end if 

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
!             call nc_write_attr(filename,var_now%nm_out,"grid_mapping",trim(grid%mtype))
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        ! Modify variables for consistency and gradient limit 

        ! Allocate helper arrays
        call grid_allocate(grid,zs)
        call grid_allocate(grid,zb)
        call grid_allocate(grid,H)
        call grid_allocate(grid,zb0)

        ! Re-load data
        call nc_read(filename,"zs",zs)
        call nc_read(filename,"zb",zb)
        call nc_read(filename,"H",H)
        call nc_read(filename,"mask_ice",outvar)
        
        ! Also load etopo bedrock, use it to replace high latitude regions 
        call nc_read(filename0,"zb",zb0)
        where (grid%lat > -60.d0) zb = zb0 

        ! Eliminate problematic regions for this domain ========
        call grid_allocate(grid,mask_reg)    
        
        ! Bad island with ice 
        allocate(xp(4),yp(4))
        xp = [-44.51, -43.51, -45.49, -44.49]
        yp = [-59.87, -60.35, -60.36, -60.85]
        mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
        where (mask_reg .and. H .gt. 0) zb = mv 
        where (mask_reg .and. H .gt. 0) zs = mv 

        ! Bad island with ice 
        xp = [167.0, 159.0, 160.0, 169.5]
        yp = [-67.7, -68.0, -64.1, -65.0]
        mask_reg = point_in_polygon(real(grid%lon),real(grid%lat),xp,yp) 
        where (mask_reg .and. H .gt. 0) zb = mv 
        where (mask_reg .and. H .gt. 0) zs = mv 

        ! Replaces problematic regions with regional mean values or zero for surface
        where (zb .eq. mv) H = 0.d0 
        where (zb .eq. mv) outvar = 0.d0 

        call fill_weighted(zb,missing_value=mv)
        call fill_weighted(zs,missing_value=mv,fill_value=0.d0)

!         ! Apply gradient limit as needed
!         if (grad_lim .gt. 0.d0) then 
!             ! Limit the gradient (m/m) to below threshold 
!             call limit_gradient(zs,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
!             call limit_gradient(zb,grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv,grad_lim=grad_lim,iter_max=50)
            
!         end if 

        ! Update mask and H 
        where (H .lt. 1.d0) H  = 0.d0 
        where (H .eq. 0.d0) zs = 0.d0 
        where (outvar .eq. 2.d0 .and. zs .gt. 0.d0) H = zs-zb 
        where (zs     .eq. 0.d0) outvar = 0.d0 

        ! Re-write fields 
        call nc_write(filename,"zs",real(zs),dim1="xc",dim2="yc",missing_value=real(mv))
        call nc_write(filename,"zb",real(zb),dim1="xc",dim2="yc",missing_value=real(mv))
        call nc_write(filename,"H", real(H), dim1="xc",dim2="yc",missing_value=real(mv))
        call nc_write(filename,"mask",nint(outvar),dim1="xc",dim2="yc",missing_value=nint(mv))

        return 

    end subroutine bedmap2_to_grid

    subroutine bedmap2vel_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       VELOCITY DATA on bedmap2 grid
        !
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename, infldr, prefix  
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid0
        character(len=256) :: file_invariant
        type(var_defs), allocatable :: invariant(:)
        double precision, allocatable :: invar(:,:), invarb(:,:)
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp1(:,:), tmp2(:,:), tmp3(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 

        ! Define input grid
        if (trim(domain) .eq. "Antarctica") then 
            
            ! Define topography (BEDMAP2) grid and input variable field
            call grid_init(grid0,name="BEDMAP2-10KM",mtype="polar_stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-3328.d0,dx=10.d0,nx=666,y0=-3328.d0,dy=10.d0,ny=666, &
                   lambda=0.d0,phi=-90.d0,alpha=24.7d0)

            ! Define the input filenames
            infldr         = "output/Antarctica/BEDMAP2-netcdf/"
            file_invariant = trim(infldr)//"ANT-1KM_BEDMAP2_vel.nc"
            desc    = "Antarctica surface velocity"
            ref     = "Rignot, E., Mouginot, J. and Scheuchl, B.: &
                      &Ice Flow of the Antarctic Ice Sheet, Science, &
                      &doi 10.1126/science.1208336, 2011."

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_VEL-R11.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(3))
        call def_var_info(invariant(1),file_invariant,  "u","ux_srf",units="m*a-1", &
                          long_name="Surface velocity, u-component",method="radius")
        call def_var_info(invariant(2),file_invariant,  "v","uy_srf",units="m*a-1", &
                          long_name="Surface velocity, v-component",method="radius")
        call def_var_info(invariant(3),file_invariant,"uv","uxy_srf",units="m*a-1", &
                          long_name="Surface velocity, magnitude",method="radius")

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)
        call grid_allocate(grid0,invarb)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp1(6667,6667))  ! bedmap2 array

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

        ! ## INVARIANT FIELDS ##
        do i = 1, size(invariant)
            var_now = invariant(i) 
            if (trim(var_now%nm_out) .eq. "uxy_srf") then 
                call nc_read(var_now%filename,"u",tmp1,missing_value=missing_value)
                call thin(invar,tmp1,by=10,missing_value=mv)
!                 call thin_ave(invar,tmp1,by=10,missing_value=mv)
                where( invar .eq. missing_value ) invar = 0.d0 
                call nc_read(var_now%filename,"v",tmp1,missing_value=missing_value)
                call thin(invarb,tmp1,by=10,missing_value=mv)
!                 call thin_ave(invarb,tmp1,by=10,missing_value=mv)
                where( invarb .eq. missing_value ) invarb = 0.d0 
                invar = dsqrt(invar**2 + invarb**2)
            else
                call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=missing_value)
                call thin(invarb,tmp1,by=10,missing_value=mv)
!                 call thin_ave(invar,tmp1,by=10,missing_value=mv)
                where( invar .eq. missing_value ) invar = 0.d0 
            end if 

            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,radius=grid%G%dx, &
                          fill=.FALSE.,missing_value=missing_value)
            call fill_mean(outvar,missing_value=missing_value)
            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc")
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc")
            end if 
        
            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
!             call nc_write_attr(filename,var_now%nm_out,"grid_mapping",trim(grid%mtype))
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine bedmap2vel_to_grid

    subroutine bedmap2acc_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       ACCUMULATION DATA on bedmap2 grid
        !
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename, infldr, prefix  
        character(len=1024) :: desc, ref 

        type(grid_class)   :: grid0
        character(len=256) :: file_invariant
        type(var_defs), allocatable :: invariant(:)
        double precision, allocatable :: invar(:,:), invarb(:,:)
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp1(:,:), tmp2(:,:), tmp3(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 

        ! Define input grid
        if (trim(domain) .eq. "Antarctica") then 
            
            ! Define topography (BEDMAP2) grid and input variable field
            call grid_init(grid0,name="BEDMAP2-10KM",mtype="polar_stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-3328.d0,dx=10.d0,nx=666,y0=-3328.d0,dy=10.d0,ny=666, &
                   lambda=0.d0,phi=-90.d0,alpha=24.7d0)

            ! Define the input filenames
            infldr         = "output/Antarctica/BEDMAP2-netcdf/"
            file_invariant = trim(infldr)//"ANT-1KM_BEDMAP2_acc.nc"
            desc    = "Antarctica climatological accumulation"
            ref     = "Arthern, R. J., Winebrenner, D. P. and Vaughan, D. G.: &
                      &Antarctic snow accumulation mapped using polarization of &
                      &4.3-cm wavelength microwave emission, J. Geophys. Res., 111, D06107, &
                      &doi:10.1029/2004JD005667, 2006."

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_ACC-A06.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(1))
        call def_var_info(invariant(1),file_invariant,  "accum","accum",units="mm*a-1", &
                          long_name="Annual accumulation")

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)
        call grid_allocate(grid0,invarb)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp1(6667,6667))  ! bedmap2 array

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
        var_now = invariant(1) 
        call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=missing_value)
        call thin(invar,tmp1,by=10,missing_value=mv)
!         call thin_ave(invar,tmp1,by=10,missing_value=mv)
!         where( invar .eq. missing_value ) invar = 0.d0 
    
        call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,radius=20.d0, &
                      fill=.TRUE.,missing_value=missing_value)
        call fill_weighted(outvar,missing_value=missing_value)
        call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc")
    
        ! Write variable metadata
        call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
        call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
!         call nc_write_attr(filename,var_now%nm_out,"grid_mapping",trim(grid%mtype))
        call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")

        return 

    end subroutine bedmap2acc_to_grid





end module bedmap2
