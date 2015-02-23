module Rignot13_BasalMelt

    use gridding_datasets
    use coordinates 
    use interp2D
    use ncio 
    
    implicit none 

    private 
    public :: Rignot13_BasalMelt_to_grid
    
contains 

    subroutine Rignot13_BasalMelt_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
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

        type(grid_class)   :: gTOPO
        character(len=256) :: file_invariant, file_surface, file_prefix(2)
        type(var_defs), allocatable :: invariant(:), surf(:), pres(:) 
        double precision, allocatable :: invar(:,:)
        integer :: plev(9) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp1(:,:), tmp2(:,:), tmp3(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: xax(:), yax(:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: nyr, nm, q, k, year, m, i, l, year0, year_switch, n_prefix, n_var 

        ! Define input grid
        if (trim(domain) .eq. "Antarctica") then 
            
            ! Define the input filenames
            infldr         = "sicodata/Antarctica_data/"
            file_invariant = trim(infldr)//"Ant_MeltingRate.nc"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_BMELT.nc"

!             ! Define input grid vectors 
!             allocate(xax(561),yax(561))
!             call nc_read(file_invariant,"xaxis",xax)
!             call nc_read(file_invariant,"yaxis",yax)
            
            ! Define topography (BEDMAP2/rignot) grid and input variable field
            call grid_init(gTOPO,name="rignot-10KM",mtype="polar stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-2800.d0,dx=10.d0,nx=561,y0=2800.d0,dy=-10.d0,ny=561, &
                   lambda=0.d0,phi=-90.d0,alpha=19.0d0)

            
        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(2))
        call def_var_info(invariant(1),file_invariant,  "melt_actual","bm_actual",units="m*a-1")
        call def_var_info(invariant(2),file_invariant,  "melt_steadystate","bm_equil",units="m*a-1")

        ! Allocate the input grid variable
        call grid_allocate(gTOPO,invar)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp1(5601,5601))  ! bedmap2-rignot array

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
            call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=missing_value)
            call thin(invar,tmp1,by=10)
            where( invar .eq. missing_value ) invar = 0.d0 

            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method,20.d3, &
                          fill=.TRUE.,missing_value=missing_value)
            call fill_mean(outvar,missing_value=missing_value)
            if (var_now%method .eq. "nn") then 
                call nc_write(filename,var_now%nm_out,nint(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            else
                call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc",units=var_now%units_out)
            end if 
        end do 

        return 

    end subroutine Rignot13_BasalMelt_to_grid



end module Rignot13_BasalMelt
