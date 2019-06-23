module AntarcticaVelocity 

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: antevelr11_to_grid

contains 

    subroutine antevelr11_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       MEaSURES Antarctica Ice Velocity Map 450m spacing
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
        double precision, allocatable :: ux(:,:), uy(:,:) 
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: nyr, nm, q, k, year, m, i, j, l, year0, year_switch, n_prefix, n_var 

        character(len=512) :: filename_H_ice 
        double precision, allocatable :: H_ice(:,:)
        
        ! Define input grid
        if (trim(domain) .eq. "Antarctica") then 
            
            ! Define input variable grid
            call grid_init(grid0,name="ANTVEL-4.5KM",mtype="polar_stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-2800.d0,dx=10.d0,nx=1246,y0=-2800.d0,dy=10.d0,ny=1246, &
                   lambda=0.d0,phi=-71.d0)

            ! Define the input filenames
            infldr         = "data/Antarctica/"
            file_invariant = trim(infldr)//"antarctica_ice_velocity_450m_v2.nc"
            desc    = "MEaSURES Antarctica Ice Velocity, version 2, https://nsidc.org/data/nsidc-0484"
            ref     = "Rignot, E., Mouginot, J. and Scheuchl, B.: &
                      &Ice Flow of the Antarctic Ice Sheet, Science, &
                      &doi 10.1126/science.1208336, 2011."

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_VEL-R11_2.nc"

            write(filename_H_ice,"(a)") trim(outfldr)//"/"// &
                              trim(grid%name)//"_TOPO-RTOPO-2.0.1.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(invariant(3))
        call def_var_info(invariant(1),file_invariant,  "U","ux_srf",units="m*a-1", &
                          long_name="Surface velocity, u-component",method="radius")
        call def_var_info(invariant(2),file_invariant,  "V","uy_srf",units="m*a-1", &
                          long_name="Surface velocity, v-component",method="radius")
        call def_var_info(invariant(3),file_invariant,"uv","uxy_srf",units="m*a-1", &
                          long_name="Surface velocity, magnitude",method="radius")

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)
        call grid_allocate(grid0,invarb)

        ! Allocate tmp array to hold full data (that will be trimmed to smaller size)
        allocate(tmp1(12445,12445))  ! 450m-array
        allocate(tmp2(12445,12445))  ! 450m-array

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        call grid_allocate(grid,H_ice)
        

        ! First load H_ice field for filtering where velocity should be zero 
        ! (assumes topography is already available on output domain)
        call nc_read(filename_H_ice,"H_ice",H_ice)


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
            if (trim(var_now%nm_out) .ne. "uxy_srf") then 
                call nc_read(var_now%filename,var_now%nm_in,tmp1,missing_value=missing_value)
                do j = 1, size(tmp2,2)
                    tmp1(:,j) = tmp2(:,size(tmp2,2)-j+1)
                end do 
                call thin(invar,tmp1,by=10,missing_value=mv)
                where( invar .eq. missing_value ) invar = 0.d0 
            end if 

            outvar = mv 

            call map_field_conservative_map1(map%map,var_now%nm_in,invar,outvar, &
                                                            method="mean",missing_value=mv)
            where(H_ice .eq. 0.0) outvar = mv 
!             call fill_mean(outvar,missing_value=missing_value)
            call nc_write(filename,var_now%nm_out,real(outvar),dim1="xc",dim2="yc")
            
            ! Write variable metadata
            call nc_write_attr(filename,var_now%nm_out,"units",var_now%units_out)
            call nc_write_attr(filename,var_now%nm_out,"long_name",var_now%long_name)
!             call nc_write_attr(filename,var_now%nm_out,"grid_mapping",trim(grid%mtype))
            call nc_write_attr(filename,var_now%nm_out,"coordinates","lat2D lon2D")
            
        end do 

        ! === Now velocity magnitude ===

        call grid_allocate(grid,ux)
        call grid_allocate(grid,uy)

        call nc_read(filename,"ux_srf",ux)
        call nc_read(filename,"uy_srf",uy)

        outvar = mv
        where(ux .ne. mv .and. uy .ne. mv) outvar = sqrt(ux**2 + uy**2)

        call nc_write(filename,"uxy_srf",outvar,dim1="xc",dim2="yc",missing_value=mv)

        ! Write variable metadata
        call nc_write_attr(filename,"uxy_srf","units","m/a")
        call nc_write_attr(filename,"uxy_srf","long_name", &
                    "Surface velocity, magnitude")
        call nc_write_attr(filename,"uxy_srf","coordinates","lat2D lon2D")

        return 

    end subroutine antevelr11_to_grid


end module AntarcticaVelocity 
