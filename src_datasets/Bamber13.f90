module Bamber13 

    use gridding_datasets
    use coordinates
    use interp2D 
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
        character(len=1024) :: desc, ref 

        type(grid_class)   :: gTOPO
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

        ! Define input grid
        if (trim(domain) .eq. "Greenland") then 
            
            ! Define topography (Bamber et al. 2013) grid and input variable field
            select case(thin_by)
                case(10)
                    call grid_init(gTOPO,name="TOPO-B13-10KM",mtype="polar_stereographic", &
                            units="kilometers",lon180=.TRUE., &
                            x0=-1300.d0,dx=10.d0,nx=251,y0=-3500.d0,dy=10.d0,ny=301, &
                            lambda=-39.d0,phi=90.d0,alpha=7.5d0)

                case(5)
                    call grid_init(gTOPO,name="TOPO-B13-5KM",mtype="polar_stereographic", &
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
                              "_TOPO-B13.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Define the variables to be mapped 
        allocate(vars(4))
        call def_var_info(vars(1),trim(file_in),"BedrockElevation","zb",units="m",long_name="Bedrock elevation")
        call def_var_info(vars(2),trim(file_in),"SurfaceElevation","zs",units="m",long_name="Surface elevation")
        call def_var_info(vars(3),trim(file_in),"IceThickness",     "H",units="m",long_name="Ice thickness")
        call def_var_info(vars(4),trim(file_in),"LandMask",      "mask",units="(0 - 4)", &
                          long_name="Land mask",method="nn")

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
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! ## FIELDS ##
        do i = 1, size(vars)
            var_now = vars(i) 
            call nc_read(trim(var_now%filename),var_now%nm_in,tmp,missing_value=mv)
            call thin(invar,tmp,by=thin_by)

            write(*,*) "1:   ", trim(var_now%nm_in), minval(invar), maxval(invar)

            if (trim(var_now%nm_out) .eq. "H" .or. trim(var_now%nm_out) .eq. "zs") then 
                where( invar .eq. mv ) invar = 0.d0 
            end if
            if (trim(var_now%nm_out) .eq. "zb") then 
                call fill_mean(invar,missing_value=mv,fill_value=-1500.d0)
            end if 
            call map_field(map,var_now%nm_in,invar,outvar,outmask,var_now%method, &
                           radius=grid%G%dx*5d0,fill=.TRUE.,missing_value=mv)
            
            write(*,*) "2:   ", trim(var_now%nm_in), minval(invar),  maxval(invar)
            write(*,*) "2:   ", trim(var_now%nm_in), minval(outvar), maxval(outvar)

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

    end subroutine Bamber13_to_grid

end module Bamber13 
