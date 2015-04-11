module AN1CRUST

    use gridding_datasets
    use coordinates
    use interp2D 
    use ncio 
    
    implicit none 

    private 
    public :: An15litho_to_grid 


contains 

    subroutine An15litho_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
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
        character(len=256) :: file_in
        type(var_defs), allocatable :: vars(:)
        double precision, allocatable :: invar(:,:) 

        type(map_class)  :: map 
        type(var_defs) :: var_now 
        double precision, allocatable :: outvar(:,:), tmp(:,:)
        integer, allocatable          :: outmask(:,:)
        double precision, allocatable :: zb(:,:), zs(:,:), H(:,:)
        integer :: q, k, m, i, l, n_var 
        integer :: thin_by = 10 
        character(len=128) :: method 

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
            file_in = "/data/sicopolis/data/Greenland/Morlighem2014_topo/MCdataset-2014-11-19.nc"
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
        allocate(vars(5))
        call def_var_info(vars(1),trim(file_in),"bed",      "zb",units="m",long_name="Bedrock elevation")
        call def_var_info(vars(2),trim(file_in),"surface",  "zs",units="m",long_name="Surface elevation")
        call def_var_info(vars(3),trim(file_in),"thickness","H",units="m",long_name="Ice thickness")
        call def_var_info(vars(4),trim(file_in),"errbed",   "zb_err",units="m",long_name="Bedrock / ice thickness error")
        call def_var_info(vars(5),trim(file_in),"mask",     "mask",units="(0 - 3)", &
                          long_name="(0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice)",method="nn")

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

        ! Make sure the file can be opened 
        write(*,*) "Reading: ",trim(file_in)
        write(*,*) "nx = ", nc_size(file_in,"x")
        write(*,*) "ny = ", nc_size(file_in,"y")
        
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
            if (trim(var_now%nm_out) .eq. "mask") method = "nn" 

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

    end subroutine An15litho_to_grid

end module AN1CRUST 
