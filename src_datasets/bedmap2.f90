module bedmap2 

    use gridding_datasets
    use coordinates 
    use interp2D
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
        double precision :: lat_lim 
        character(len=512) :: filename, infldr, prefix  
        character(len=1024) :: desc, ref 

        type(grid_class)   :: gTOPO
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

        ! Define input grid
        if (trim(domain) .eq. "Antarctica") then 
            
            ! Define topography (BEDMAP2) grid and input variable field
            call grid_init(gTOPO,name="BEDMAP2-10KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-3400.d0,dx=10.d0,nx=681,y0=-3400.d0,dy=10.d0,ny=681, &
                   lambda=0.d0,phi=-90.d0,alpha=19.0d0)

            ! Define the input filenames
            infldr         = "output/Antarctica/"
            file_invariant = trim(infldr)//"ANT-1KM_BEDMAP2_topo.nc"
            desc    = "Antarctica bedrock and surface topography (BEDMAP2)"
            ref     = "Fretwell et al.: Bedmap2: improved ice bed, surface and &
                      &thickness datasets for Antarctica, The Cryosphere, 7, 375-393, &
                      &doi:10.5194/tc-7-375-2013, 2013."

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO-BEDMAP2.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(4))
        call def_var_info(invariant(1),file_invariant,"zs","zs",units="m",long_name="Surface elevation")
        call def_var_info(invariant(2),file_invariant,"zb","zb",units="m",long_name="Bedrock elevation")
        call def_var_info(invariant(3),file_invariant,"H","H",units="m",long_name="Ice thickness")
        call def_var_info(invariant(4),file_invariant,"mask_ice","mask_ice",units="(0 - 1)", &
                          long_name="Ice mask",method="nn")

        ! Allocate the input grid variable
        call grid_allocate(gTOPO,invar)
        
        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp1(6667,6667))  ! bedmap2 array

        ! Initialize mapping
        call map_init(map,gTOPO,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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
            call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=missing_value)
            call thin(invar,tmp1,by=10)
            if (trim(var_now%nm_out) .eq. "H" .or. &
                trim(var_now%nm_out) .eq. "zs") then 
                where( invar .eq. missing_value ) invar = 0.d0 
            end if
            if (trim(var_now%nm_out) .eq. "zb") then 
                call fill_mean(invar,missing_value=missing_value,fill_value=-1001.d0)
!                 call fill_mean(invar,missing_value=missing_value)
            end if 
            if (trim(var_now%nm_out) .eq. "mask_ice") then 
                where ( invar .eq. 1.d0 ) invar = 3.d0 
                where ( invar .eq. 0.d0 ) invar = 2.d0 
                where ( invar .eq. missing_value ) invar = 0.d0 
            end if 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,20.d3, &
                          fill=.TRUE.,missing_value=missing_value)
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

        type(grid_class)   :: gTOPO
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
            call grid_init(gTOPO,name="BEDMAP2-10KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-3400.d0,dx=10.d0,nx=681,y0=-3400.d0,dy=10.d0,ny=681, &
                   lambda=0.d0,phi=-90.d0,alpha=19.0d0)

            ! Define the input filenames
            infldr         = "output/Antarctica/"
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
        call def_var_info(invariant(1),file_invariant,  "u","u",units="m*a-1", &
                          long_name="Surface velocity, u-component")
        call def_var_info(invariant(2),file_invariant,  "v","v",units="m*a-1", &
                          long_name="Surface velocity, v-component")
        call def_var_info(invariant(3),file_invariant,"uv","uv",units="m*a-1", &
                          long_name="Surface velocity, magnitude")

        ! Allocate the input grid variable
        call grid_allocate(gTOPO,invar)
        call grid_allocate(gTOPO,invarb)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp1(6667,6667))  ! bedmap2 array

        ! Initialize mapping
        call map_init(map,gTOPO,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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
            if (trim(var_now%nm_out) .eq. "uv") then 
                call nc_read(var_now%filename,"u",tmp1,missing_value=missing_value)
                call thin(invar,tmp1,by=10)
                where( invar .eq. missing_value ) invar = 0.d0 
                call nc_read(var_now%filename,"v",tmp1,missing_value=missing_value)
                call thin(invarb,tmp1,by=10)
                where( invarb .eq. missing_value ) invarb = 0.d0 
                invar = dsqrt(invar**2 + invarb**2)
            else
                call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=missing_value)
                call thin(invar,tmp1,by=10)
                where( invar .eq. missing_value ) invar = 0.d0 
            end if 

            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,20.d3, &
                          fill=.TRUE.,missing_value=missing_value)
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

        type(grid_class)   :: gTOPO
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
            call grid_init(gTOPO,name="BEDMAP2-10KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-3400.d0,dx=10.d0,nx=681,y0=-3400.d0,dy=10.d0,ny=681, &
                   lambda=0.d0,phi=-90.d0,alpha=19.0d0)

            ! Define the input filenames
            infldr         = "output/Antarctica/"
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
        call grid_allocate(gTOPO,invar)
        call grid_allocate(gTOPO,invarb)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp1(6667,6667))  ! bedmap2 array

        ! Initialize mapping
        call map_init(map,gTOPO,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

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
        call thin(invar,tmp1,by=10)
!         where( invar .eq. missing_value ) invar = 0.d0 
    
        call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,20.d3, &
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
