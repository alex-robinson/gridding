module Rignot13_BasalMelt

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: Rignot13_BasalMelt_to_grid
    
contains 

    subroutine Rignot13_BasalMelt_to_grid(outfldr,grid,domain,max_neighbors,lat_lim,fill,sigma)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       BASAL MELTING DATA on bedmap2 grid
        !       Reference: 
        !           Rignot, E., Jacobs, S., Mouginot, J. and 
        !           Scheuchl, B., Ice-Shelf Melting Around Antarctica, 
        !           Science, doi: 10.1126/science.1235798, 2013
        ! =========================================================
        
        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim 
        character(len=512) :: filename, infldr, prefix  
        character(len=1024) :: desc, ref 
        double precision :: sigma 
        logical :: fill

        type(grid_class)   :: grid0
        character(len=256) :: file_in
        type(var_defs), allocatable :: vars(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:)
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp1(:,:), tmp2(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: xax(:), yax(:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 
        character(len=256) :: method, method_str, fill_str 

        character(len=56)  :: nm_out 
        character(len=512) :: file_basins 
        integer, allocatable :: basins(:,:)
        double precision :: basin_ave 
        logical, allocatable :: mask_basin(:,:) 

        ! Determine the method to use here 
        method = "nn" 
        if (sigma .gt. 0.d0) method = "nng" 

        method_str = "nn" 
        if (sigma .gt. 0.d0) then 
            if (sigma .gt. 99.d0) then 
                write(method_str,"(a,i3)") trim(method), int(sigma) 
            else
                write(method_str,"(a,i2)") trim(method), int(sigma) 
            end if 
        end if 

        fill_str = "_nofill"
        if (fill) fill_str = "_fill" 

            ! Define the input filenames
            infldr  = "/data/sicopolis/data/Antarctica/"
            file_in = trim(infldr)//"Ant_MeltingRate.nc"

            desc    = "Ice shelf basal melting dataset"
            ref     = "Rignot, E., Jacobs, S., Mouginot, J. and Scheuchl, B.: &
                      &Ice-Shelf Melting Around Antarctica, Science, &
                      &doi: 10.1126/science.1235798, 2013."

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_BMELT-R13.nc"

            ! Define topography (BEDMAP2/rignot) grid and input variable field
            call grid_init(grid0,name="rignot-10KM",mtype="polar_stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-2795.d0,dx=10.d0,nx=561,y0=2795.d0,dy=-10.d0,ny=561, &
                   lambda=0.d0,phi=-90.d0,alpha=19.0d0)

        ! Define the variables to be mapped 
        allocate(vars(2))
        call def_var_info(vars(1),file_in,"melt_actual","bm_actual",units="m*a**-1", &
                          long_name="Basal melt rate, actual present day",method=trim(method))
        call def_var_info(vars(2),file_in,"melt_steadystate","bm_equil",units="m*a**-1", &
                          long_name="Basal melt rate, shelf equilibrium",method=trim(method))

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp1(5601,5601))  ! bedmap2-rignot array

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")

        ! Write the standard grid information (x2D,y2D,lon2D,lat2D,area,etc)
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write dataset meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ## INVARIANT FIELDS ##
        do i = 1, size(vars)
            var_now = vars(i) 
            call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=mv)
            call thin(invar,tmp1,by=10,missing_value=mv)
!             call thin_ave(invar,tmp1,by=10,missing_value=mv)   ! Diffuses too much - check?
            where( invar .eq. 0.d0 ) invar = mv

            ! Make sure outvar is initialized with missing values 
            outvar = mv 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,"nn",radius=40.d0, &
                          fill=.FALSE.,missing_value=mv,sigma=sigma)
            where(outvar .eq. mv) outvar = 0.d0 

            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))

            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        ! Now output some variants to the original data
        call grid_allocate(grid,basins)    
        call grid_allocate(grid,mask_basin)    
        
        ! First load basin mask 
        write(file_basins,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_BASINS-nasa.nc"
        write(*,*) "Reading basin mask: ", trim(file_basins)
        call nc_read(file_basins,"basin",basins)

        ! ===== Smoothed with extrapolation of mean basin value ===== 
        do i = 1, size(vars)

            var_now = vars(i) 
            call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=mv)
            call thin(invar,tmp1,by=10,missing_value=mv)
!             call thin_ave(invar,tmp1,by=10,missing_value=mv)   ! Diffuses too much - check?
            where( invar .eq. 0.d0 ) invar = mv

            ! Make sure outvar is initialized with missing values 
            outvar = mv 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,"nng",radius=40.d0, &
                          fill=.FALSE.,missing_value=mv,sigma=sigma)
            
            do q = 1, maxval(basins)
                mask_basin = (outvar .ne. mv) .and. (basins .eq. q)
                basin_ave = mv
                if (count(mask_basin) .gt. 0) then 
                    basin_ave = sum(outvar,mask=mask_basin) / count(mask_basin)
                end if 
                where ( (outvar .eq. mv) .and. (basins .eq. q) ) outvar = basin_ave 
            end do 

            ! Fill in missing values 
            call fill_mean(outvar,missing_value=mv)

            nm_out = trim(var_now%nm_out)//"_sm"
            call nc_write(filename,nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name", &
                                trim(var_now%long_name)//" (nng + basin average)")
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
                
            ! ===== Mean basin value everywhere, then smoothed ==========
            var_now = vars(i) 
            call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=mv)
            call thin(invar,tmp1,by=10)
            where( invar .eq. 0.d0 ) invar = mv

            ! Make sure outvar is initialized with missing values 
            outvar = mv 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,"nn",radius=40.d0, &
                          fill=.FALSE.,missing_value=mv,sigma=sigma)

            do q = 1, maxval(basins)
                mask_basin = (outvar .ne. mv) .and. (basins .eq. q)
                basin_ave = mv
                if (count(mask_basin) .gt. 0) then 
                    basin_ave = sum(outvar,mask=mask_basin) / count(mask_basin)
                end if 
                where ( (basins .eq. q) ) outvar = basin_ave 
            end do 

            ! Fill in missing values 
            call fill_mean(outvar,missing_value=mv)

!             call grid_allocate(grid,tmp2)
!             call filter_gaussian(var=outvar,sigma=sigma,dx=grid%G%dx,mask=tmp2.ne.mv)
        
            nm_out = trim(var_now%nm_out)//"_ave"
            call nc_write(filename,nm_out,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name", &
                                trim(var_now%long_name)//" (basin average)")
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        return 

    end subroutine Rignot13_BasalMelt_to_grid

end module Rignot13_BasalMelt
