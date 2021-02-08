module Morlighem2013

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: morlighem2013_taub_to_grid

contains

    subroutine morlighem2013_taub_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       Morlighem et al (2013) Antarctica inversion data 
        !         - Basal stress 
        !         - 3D velocity field 
        !         - f_shr 
        !
        ! =========================================================
        
        implicit none 

        character(len=*), intent(IN) :: outfldr 
        type(grid_class), intent(IN) :: grid 
        character(len=*), intent(IN) :: domain 
        integer,          intent(IN) :: max_neighbors 
        double precision, intent(IN) :: lat_lim


        ! Local variables 
        character(len=512) :: filename, file_input  
        character(len=1024) :: desc, ref 

        type(points_class) :: pts0 

        type points_vector_type 
            double precision, allocatable :: lon(:), lat(:), var(:)
            double precision, allocatable :: taub(:), taud(:), f_vbvs(:)
        end type 

        type(points_vector_type) :: inp
        integer :: i, np, k, nz 
        character(len=256) :: vnm_now 
        character(len=256) :: units_now 
        character(len=256) :: long_name_now 
        
        type(map_class)  :: map 

        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)

        if (trim(domain) .ne. "Antarctica") then 
            write(*,*) "Error: Morlighem 2013 dataset is only valid for domain=Antarctica."
            stop 
        end if 

        ! ============================================================================
        !
        ! Input data 
        !
        ! ============================================================================
        
        ! Define the input data information
        file_input = "/data/sicopolis/data/Antarctica/Morlighem2013/Morlighem2013_data.txt"
        desc = "Antarctic basal stress inversion."
        ref  = "Morlighem et al.: Inversion of basal friction in Antarctica using exact and &
               &incomplete adjoints of a higher-order model, JGR: Earth Surface, 118, 1746–1753, &
               &doi:10.1002/jgrf.20125, 2013"

        ! Number of data points 
        np = 32933

        allocate(inp%lon(np),inp%lat(np),inp%var(np))
        allocate(inp%taub(np),inp%taud(np),inp%f_vbvs(np))

        ! File format: lat, lon, basin 
        open(2,file=trim(file_input),status="old")
        do i = 1, np 
            read(2,*) inp%lat(i), inp%lon(i), inp%taub(i), inp%taud(i), inp%f_vbvs(i) 
        end do 
        close(2)

        ! Convert f_vbvs to fraction 
        inp%f_vbvs = inp%f_vbvs/100.d0 

        write(*,*) "lon:    ", minval(inp%lon),    maxval(inp%lon)
        write(*,*) "lat:    ", minval(inp%lat),    maxval(inp%lat)
        write(*,*) "taub:   ", minval(inp%taub),   maxval(inp%taub)
        write(*,*) "taud:   ", minval(inp%taud),   maxval(inp%taud)
        write(*,*) "f_vbvs: ", minval(inp%f_vbvs), maxval(inp%f_vbvs)

        ! Define input points for mapping
        call points_init(pts0,name="M13-taub",mtype="latlon",units="degrees",x=inp%lon,y=inp%lat,lon180=.TRUE.)
        

        ! ============================================================================
        !
        ! Output data 
        !
        ! ============================================================================
        
        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)//"_TAUB-M13.nc"

        ! Initialize mapping
        call map_init(map,pts0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Define number of vertical layers 
        nz = 14 

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="km")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="km")
        call nc_write_dim(filename,"zeta", x=0.d0,dx=1.d0/dble(nz-1),nx=nz,units="1")
        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! === tau_b ===

        inp%var       = inp%taub 
        vnm_now       = "tau_b"
        units_now     = "kPa"
        long_name_now = "Basal shear stress"
        
        ! ## MAP FIELD ##
        outvar = mv 
        call map_field(map,vnm_now,inp%var,outvar,outmask,"nn",fill=.FALSE.,sigma=grid%G%dx,missing_value=mv)

        write(*,*) "Range invar:  ",minval(inp%var,inp%var.ne.mv), maxval(inp%var,inp%var.ne.mv)
        write(*,*) "Range outvar: ",minval(outvar,outvar.ne.mv),   maxval(outvar,outvar.ne.mv)
        
        ! Fill any missing values
        !call fill_weighted(outvar,missing_value=missing_value)
    
        ! Write field to output file 
        call nc_write(filename,vnm_now,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))

        ! Write variable metadata
        call nc_write_attr(filename,vnm_now,"units",trim(units_now))
        call nc_write_attr(filename,vnm_now,"long_name",trim(long_name_now))
        call nc_write_attr(filename,vnm_now,"coordinates","lat2D lon2D")
        
        ! === tau_d ===

        inp%var       = inp%taud 
        vnm_now       = "tau_d"
        units_now     = "kPa"
        long_name_now = "Driving stress"
        
        ! ## MAP FIELD ##
        outvar = mv 
        call map_field(map,vnm_now,inp%var,outvar,outmask,"nn",fill=.FALSE.,sigma=grid%G%dx,missing_value=mv)

        write(*,*) "Range invar:  ",minval(inp%var,inp%var.ne.mv), maxval(inp%var,inp%var.ne.mv)
        write(*,*) "Range outvar: ",minval(outvar,outvar.ne.mv),   maxval(outvar,outvar.ne.mv)
        
        ! Fill any missing values
        !call fill_weighted(outvar,missing_value=missing_value)
    
        ! Write field to output file 
        call nc_write(filename,vnm_now,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))

        ! Write variable metadata
        call nc_write_attr(filename,vnm_now,"units",trim(units_now))
        call nc_write_attr(filename,vnm_now,"long_name",trim(long_name_now))
        call nc_write_attr(filename,vnm_now,"coordinates","lat2D lon2D")
        
        ! === f_vbvs ===

        inp%var       = inp%f_vbvs 
        vnm_now       = "f_vbvs"
        units_now     = "1"
        long_name_now = "Basal-to-surface velocity ratio"
        
        ! ## MAP FIELD ##
        outvar = mv 
        call map_field(map,vnm_now,inp%var,outvar,outmask,"nn",fill=.FALSE.,sigma=grid%G%dx,missing_value=mv)

        write(*,*) "Range invar:  ",minval(inp%var,inp%var.ne.mv), maxval(inp%var,inp%var.ne.mv)
        write(*,*) "Range outvar: ",minval(outvar,outvar.ne.mv),   maxval(outvar,outvar.ne.mv)
        
        ! Fill any missing values
        !call fill_weighted(outvar,missing_value=missing_value)
    
        ! Write field to output file 
        call nc_write(filename,vnm_now,real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))

        ! Write variable metadata
        call nc_write_attr(filename,vnm_now,"units",trim(units_now))
        call nc_write_attr(filename,vnm_now,"long_name",trim(long_name_now))
        call nc_write_attr(filename,vnm_now,"coordinates","lat2D lon2D")
        
        ! === 3D velocity, ux (14 layers) ===

        vnm_now       = "ux"
        units_now     = "m/yr"
        long_name_now = "Velocity, x-direction"
        
        file_input = "/data/sicopolis/data/Antarctica/Morlighem2013/Morlighem2013_vx.txt"
        
        do k = 1, nz 

            ! Load data for this layer (from given column in ascii file)
            inp%var = read_vector(file_input,np,col=nz-k+1,skip=0)

            ! ## MAP FIELD ##
            outvar = mv 
            call map_field(map,vnm_now,inp%var,outvar,outmask,"nn",fill=.FALSE.,sigma=grid%G%dx,missing_value=mv)

            write(*,*) "Range invar:  ",minval(inp%var,inp%var.ne.mv), maxval(inp%var,inp%var.ne.mv)
            write(*,*) "Range outvar: ",minval(outvar,outvar.ne.mv),   maxval(outvar,outvar.ne.mv)
            
            ! Fill any missing values
            !call fill_weighted(outvar,missing_value=missing_value)
        
            ! Write field to output file 
            call nc_write(filename,vnm_now,real(outvar),dim1="xc",dim2="yc",dim3="zeta",missing_value=real(mv), &
                            start=[1,1,k],count=[grid%G%nx,grid%G%nx,1])

        end do 

        ! Write variable metadata
        call nc_write_attr(filename,vnm_now,"units",trim(units_now))
        call nc_write_attr(filename,vnm_now,"long_name",trim(long_name_now))
        call nc_write_attr(filename,vnm_now,"coordinates","lat2D lon2D")
        
        ! === 3D velocity, ux (14 layers) ===

        vnm_now       = "uy"
        units_now     = "m/yr"
        long_name_now = "Velocity, y-direction"
        
        file_input = "/data/sicopolis/data/Antarctica/Morlighem2013/Morlighem2013_vy.txt"
        
        do k = 1, nz 

            ! Load data for this layer (from given column in ascii file)
            inp%var = read_vector(file_input,np,col=nz-k+1,skip=0)

            ! ## MAP FIELD ##
            outvar = mv 
            call map_field(map,vnm_now,inp%var,outvar,outmask,"nn",fill=.FALSE.,sigma=grid%G%dx,missing_value=mv)

            write(*,*) "Range invar:  ",minval(inp%var,inp%var.ne.mv), maxval(inp%var,inp%var.ne.mv)
            write(*,*) "Range outvar: ",minval(outvar,outvar.ne.mv),   maxval(outvar,outvar.ne.mv)
            
            ! Fill any missing values
            !call fill_weighted(outvar,missing_value=missing_value)
        
            ! Write field to output file 
            call nc_write(filename,vnm_now,real(outvar),dim1="xc",dim2="yc",dim3="zeta",missing_value=real(mv), &
                            start=[1,1,k],count=[grid%G%nx,grid%G%nx,1])

        end do 

        ! Write variable metadata
        call nc_write_attr(filename,vnm_now,"units",trim(units_now))
        call nc_write_attr(filename,vnm_now,"long_name",trim(long_name_now))
        call nc_write_attr(filename,vnm_now,"coordinates","lat2D lon2D")
        
        return 

    end subroutine morlighem2013_taub_to_grid

end module Morlighem2013
