module Bamber13 

    use gridding_datasets
    use coordinates 
    use ncio 
    
    implicit none 

    private 
    public :: Bamber13_to_grid
    
contains 

    subroutine Bamber13_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
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
        character(len=512) :: filename 

        type(grid_class)   :: gTOPO
        character(len=256) :: file_invariant, file_surface, file_prefix(2)
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:) 
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 

        ! Define ECMWF input grid
        if (trim(domain) .eq. "Greenland") then 
            
            ! Define topography (Bamber et al. 2013) grid and input variable field
            call grid_init(gTOPO,name="TOPO-10KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-1300.d0,dx=10.d0,nx=251,y0=-3500.d0,dy=10.d0,ny=301, &
                   lambda=-39.d0,phi=90.d0,alpha=7.5d0)

            ! Define the input filenames
            file_invariant = "data/Greenland/Greenland_bedrock_topography_V3.nc"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_TOPO.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(4))
        call def_var_info(invariant(1),trim(file_invariant),"BedrockElevation","zb",units="m")
        call def_var_info(invariant(2),trim(file_invariant),"SurfaceElevation","zs",units="m")
        call def_var_info(invariant(3),trim(file_invariant),"IceThickness",    "H", units="m")
        call def_var_info(invariant(4),trim(file_invariant),"LandMask",      "mask",units="(0 - 4",method="nn")

        ! Allocate the input grid variable
        call grid_allocate(gTOPO,invar)
        
        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp(2501,3001))

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
    
        ! ## INVARIANT FIELDS ##
        do i = 1, size(invariant)
            var_now = invariant(i) 
            call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=missing_value)
            call thin(invar,tmp,by=10)
            if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
                where( invar .eq. missing_value ) invar = 0.d0 
            end if
            if (trim(var_now%nm_out) .eq. "zb") then 
                call fill_mean(invar,missing_value=missing_value,fill_value=-1000.d0)
            end if 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,20.d3, &
                          fill=.TRUE.,missing_value=missing_value)
            call fill_mean(outvar,missing_value=missing_value)
            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            end if 
        end do 

!         ! Fix the mask to be consistent with interpolated fields 
!         ! Initialize variable arrays
!         call grid_allocate(grid,zb)
!         call grid_allocate(grid,zs)
!         call grid_allocate(grid,H)
    
!         call nc_read(trim(filename),"zb",zb,missing_value=missing_value)
!         call nc_read(trim(filename),"zs",zs,missing_value=missing_value)
!         call nc_read(trim(filename),"H",H,missing_value=missing_value)
        
!         where(zs .lt. 0.d0) zs = 0.d0 
!         where(zs .lt. zb)   zs = zb 
!         H = zs - zb 
!         where(H  .lt. 1.d0) H  = 0.d0 

!         outvar = 0.d0 
!         where (zs .gt. 0.d0) outvar = 1.d0 
        
        ! Also perform interpolations to get drainage basins
        call basins_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)

        return 

    end subroutine Bamber13_to_grid

end module Bamber13 
