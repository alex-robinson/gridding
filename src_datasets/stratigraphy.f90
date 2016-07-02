module stratigraphy

    use gridding_datasets
    use coordinates
    use interp2D 
    use ncio 
    
    implicit none 

    private 
    public :: MacGregor15_to_grid


contains 

    subroutine MacGregor15_to_grid(outfldr,grid,domain,max_neighbors,lat_lim)
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
        integer :: q, k, m, i, l, n_var 
        integer :: thin_by = 10 
        character(len=128) :: method 

        real(4) :: age_iso(4), depth_norm(25) 

        ! Define input grid
        if (trim(domain) .eq. "Greenland") then 
            
            ! Define grid and input variable field
            call grid_init(grid0,name="ESPG-3413-1.0KM",mtype="polar_stereographic", &
                    units="kilometers",lon180=.TRUE., &
                    x0=-632.000d0,dx=1.0d0,nx=1479,y0=-3344.000d0,dy=1.0d0,ny=2675, &
                    lambda=-45.d0,phi=70.d0,alpha=20.d0)

            ! Define the input filenames
            file_in = "/data/sicopolis/data/Greenland/MacGregor2015_stratigraphy/Greenland_age_grid_zyx.nc"
            desc    = "Radiostratigraphy and age structure of the Greenland ice sheet, 24-Jun-2015 (v1.2)"
            ref     = "MacGregor, J. A., Fahnestock, M. A., Catania, G. A., Paden, J. D., &
                      &Prasad Gogineni, S., Young, S. K., Rybarski, S. C., Mabrey, A. N., &
                      &Wagman, B. M. and Morlighem, M.: Radiostratigraphy and age structure &
                      &of the Greenland Ice Sheet., J. Geophys. Res. Earth Surf., 120(2), 212–241, &
                      &doi:10.1002/2014JF003215, 2015. \n&
                      &http://www-udc.ig.utexas.edu/external/joemac/Greenland_radiostratigraphy.html"

            ! Define the output filename 
            write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                              "_STRAT-M15.nc"

        else

            write(*,*) "Domain not recognized: ",trim(domain)
            stop 
        end if 

        ! Allocate the input grid variable
        call grid_allocate(grid0,invar)

        ! Allocate tmp array to hold full data (that will potentially be trimmed to smaller size)
        allocate(tmp(1479,2675))

        ! Initialize mapping
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)    
        
        ! Get some additional dimension info from the input file 
        call nc_read(file_in,"age_iso",   age_iso)
        call nc_read(file_in,"depth_norm",depth_norm)


        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename,"xc",   x=grid%G%x,units="kilometers")
        call nc_write_dim(filename,"yc",   x=grid%G%y,units="kilometers")
        call nc_write_dim(filename,"age_iso",   x=age_iso*1e-3,units="kiloyears")
        call nc_write_dim(filename,"depth_norm",x=depth_norm,  units="1")

        call grid_write(grid,filename,xnm="xc",ynm="yc",create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! First process normalized ages
        do q = 1, size(depth_norm)

            call nc_read(file_in,"age_norm",invar,missing_value=mv, &
                         start=[1,1,q],count=[grid0%G%nx,grid0%G%ny,1])
            where (abs(invar) .gt. 1d8) invar = mv 
            where (isnan(invar)) invar = mv 

            outvar = mv 
            call map_field(map,"age_norm",invar,outvar,outmask,method="radius", &
                           radius=grid%G%dx*grid%xy_conv,fill=.FALSE.,missing_value=mv)
            where(outvar .ne. mv) outvar = outvar*1d-3
            
            call nc_write(filename,"ice_age",real(outvar),dim1="xc",dim2="yc",dim3="depth_norm", &
                          missing_value=real(mv),start=[1,1,q],count=[grid%G%nx,grid%G%ny,1])
            
            call nc_read(file_in,"age_norm_uncert",invar,missing_value=mv, &
                         start=[1,1,q],count=[grid0%G%nx,grid0%G%ny,1])
            where (abs(invar) .gt. 1d8) invar = mv 
            where (isnan(invar)) invar = mv 

            outvar = mv 
            call map_field(map,"ice_age_err",invar,outvar,outmask,method="radius", &
                           radius=grid%G%dx*grid%xy_conv,fill=.FALSE.,missing_value=mv)
            where(outvar .ne. mv) outvar = outvar*1d-3

            call nc_write(filename,"ice_age_err",real(outvar),dim1="xc",dim2="yc",dim3="depth_norm", &
                          missing_value=real(mv),start=[1,1,q],count=[grid%G%nx,grid%G%ny,1])
        
        end do 

        ! Write variable metadata
        call nc_write_attr(filename,"ice_age","units","kiloyears")
        call nc_write_attr(filename,"ice_age","long_name","Age of ice at normalized depth")
        call nc_write_attr(filename,"ice_age","coordinates","lat2D lon2D")
        
        ! Write variable metadata
        call nc_write_attr(filename,"ice_age_err","units","kiloyears")
        call nc_write_attr(filename,"ice_age_err","long_name","Error in age of ice at normalized depth")
        call nc_write_attr(filename,"ice_age_err","coordinates","lat2D lon2D")
            
  

        ! Process isochronal depths and errors
        do q = 1, size(age_iso)

            call nc_read(file_in,"depth_iso",invar,missing_value=mv, &
                         start=[1,1,q],count=[grid0%G%nx,grid0%G%ny,1])
            where (abs(invar) .gt. 1d8) invar = mv 
            where (isnan(invar)) invar = mv 
            
            outvar = mv 
            call map_field(map,"depth_iso",invar,outvar,outmask,method="radius", &
                           radius=grid%G%dx*grid%xy_conv,fill=.FALSE.,missing_value=mv)
            
            call nc_write(filename,"depth_iso",real(outvar),dim1="xc",dim2="yc",dim3="age_iso", &
                          missing_value=real(mv),start=[1,1,q],count=[grid%G%nx,grid%G%ny,1])
            
            call nc_read(file_in,"depth_iso_uncert",invar,missing_value=mv, &
                         start=[1,1,q],count=[grid0%G%nx,grid0%G%ny,1])
            where (abs(invar) .gt. 1d8) invar = mv 
            where (isnan(invar)) invar = mv 
            
            outvar = mv 
            call map_field(map,"depth_iso_err",invar,outvar,outmask,method="radius", &
                           radius=grid%G%dx*grid%xy_conv,fill=.FALSE.,missing_value=mv)
            
            call nc_write(filename,"depth_iso_err",real(outvar),dim1="xc",dim2="yc",dim3="age_iso", &
                          missing_value=real(mv),start=[1,1,q],count=[grid%G%nx,grid%G%ny,1])
        
        end do 

        ! Write variable metadata
        call nc_write_attr(filename,"depth_iso","units","normalized depth")
        call nc_write_attr(filename,"depth_iso","long_name","Depth of isochronal layer")
        call nc_write_attr(filename,"depth_iso","coordinates","lat2D lon2D")
        
        ! Write variable metadata
        call nc_write_attr(filename,"depth_iso_err","units","normalized depth")
        call nc_write_attr(filename,"depth_iso_err","long_name","Error in depth of isochronal layer")
        call nc_write_attr(filename,"depth_iso_err","coordinates","lat2D lon2D")
        
        ! Process reference thickness too
        call nc_read(file_in,"thick",invar,missing_value=mv,start=[1,1],count=[grid0%G%nx,grid0%G%ny])
        where (abs(invar) .gt. 1d8) invar = mv 
        where (isnan(invar)) invar = mv 
        
        outvar = mv 
        call map_field(map,"thick",invar,outvar,outmask,method="radius", &
                       radius=grid%G%dx*grid%xy_conv,fill=.FALSE.,missing_value=mv)
        where(outvar .eq. mv) outvar = 0.d0 

        call nc_write(filename,"H",real(outvar),dim1="xc",dim2="yc",missing_value=real(mv))
        
        ! Write variable metadata
        call nc_write_attr(filename,"depth_iso","units","m")
        call nc_write_attr(filename,"depth_iso","long_name","Ice thickness")
        call nc_write_attr(filename,"depth_iso","coordinates","lat2D lon2D")
        
        return 

    end subroutine MacGregor15_to_grid

end module stratigraphy 
