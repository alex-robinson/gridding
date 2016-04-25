module nasaBasins 

    use gridding_datasets
    use coordinates 
    use polygons 
    use index 
    use interp_time
    use interp2D 
    use ncio 
    
    implicit none 

    private 
    public :: nasaBasins_to_grid
    
contains 


    subroutine nasaBasins_to_grid(outfldr,grid,domain)

        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid  
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 

        type(points_class) :: pTOPO
        character(len=256) :: file_invariant, tmp 

        type basin_type 
            double precision, allocatable :: lon(:), lat(:), basin(:)
        end type 

        type(basin_type) :: inb
        double precision, allocatable :: outvar(:,:)
        integer, allocatable :: outmask(:,:)

        integer :: nh, nl, np, i, j, k  

        double precision, allocatable :: basins(:) 
        integer, allocatable :: inds(:)
        integer :: nb, q 
        logical :: in_basin 
        double precision, parameter :: tol = 1d-5 

        ! Get input data (from polygon files)
        if (trim(domain) .eq. "Antarctica") then 

            ! Define the input filenames
            file_invariant = "/data/sicopolis/data/Antarctica/nasa_basins/Ant_Full_DrainageSystem_Polygons.txt"
            desc = "Antarctic drainage basins mapped by NASA."
            ref  = "Zwally, H. Jay, Mario B. Giovinetto, Matthew A. Beckley, &
                   &and Jack L. Saba, 2012, Antarctic and Greenland Drainage &
                   & Systems, GSFC Cryospheric Sciences Laboratory, at &
                   http://icesat4.gsfc.nasa.gov/cryo_data/ant_grn_drainage_systems.php."

            nh = 7        ! Header length
            nl = 901329   ! File length 
            np = nl - nh  ! Number of data points 

            allocate(inb%lon(np),inb%lat(np),inb%basin(np))

            ! File format: lat, lon, basin 
            open(2,file=trim(file_invariant),status="old")
            do i = 1, nh
                read(2,"(a)") tmp 
            end do 
            do i = 1, np 
                read(2,*) inb%lat(i), inb%lon(i), inb%basin(i) 
            end do 
            close(2)

            write(*,*) "lon: ",minval(inb%lon),maxval(inb%lon)
            write(*,*) "lat: ",minval(inb%lat),maxval(inb%lat)
            write(*,*) "var: ",minval(inb%basin),maxval(inb%basin)

            ! Define input points for mapping
            call points_init(pTOPO,grid,name="NASA-ANT",x=inb%lon,y=inb%lat,latlon=.TRUE.)

        else if (trim(domain) .eq. "Greenland") then 

            ! Define the input filenames
            file_invariant = "/data/sicopolis/data/Greenland/nasa_basins/GrnDrainageSystems_Ekholm.txt"
            desc = "Greenland drainage basins mapped by NASA."
            ref  = "Zwally, H. Jay, Mario B. Giovinetto, Matthew A. Beckley, &
                   &and Jack L. Saba, 2012, Antarctic and Greenland Drainage &
                   & Systems, GSFC Cryospheric Sciences Laboratory, at &
                   http://icesat4.gsfc.nasa.gov/cryo_data/ant_grn_drainage_systems.php."

            nh = 7        ! Header length
            nl = 272972   ! File length 
            np = nl - nh  ! Number of data points 

            allocate(inb%lon(np),inb%lat(np),inb%basin(np))

            ! File format: basin, lat, lon 
            open(2,file=trim(file_invariant),status="old")
            do i = 1, nh
                read(2,"(a)") tmp 
            end do 
            do i = 1, np 
                read(2,*) inb%basin(i) , inb%lat(i), inb%lon(i)
            end do 
            close(2)

            write(*,*) "lon: ",minval(inb%lon),maxval(inb%lon)
            write(*,*) "lat: ",minval(inb%lat),maxval(inb%lat)
            write(*,*) "var: ",minval(inb%basin),maxval(inb%basin)

            ! Define input points for mapping
            call points_init(pTOPO,grid,name="NASA-GRL",x=inb%lon,y=inb%lat,latlon=.TRUE.)

        else 

            write(*,*) "nasaBasins_to_grid:: error: "
            write(*,*) "Domain not recognized: ",trim(domain)
            stop 

        end if 

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                          "_BASINS-nasa.nc"

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

        ! ## MAP FIELDS ##
        ! Map polygons onto new grid points 

        ! Initially set all output values to missing
        outvar = missing_value 

        ! Determine the input basins available 
        call unique(basins,inb%basin)
        nb = size(basins)

        ! Go through each output point and determine if it fits inside a polygon 
        write(*,*) "Mapping polygons..." 
        k = 0 
        do i = 1, grid%G%nx 
            do j = 1, grid%G%ny  

                ! Initialize basin info for current point
                in_basin      = .FALSE.

                ! Loop over basins and check point in polygon
                do q = 1, nb 
                    call which(abs(inb%basin-basins(q)).lt.tol,inds)
                    in_basin = point_in_polygon(real(grid%x(i,j)), real(grid%y(i,j)), &
                                                real(pTOPO%x(inds)), real(pTOPO%y(inds))) 
                    if (in_basin) exit 
                end do 

                ! If basin was found, save it
                if (in_basin) outvar(i,j) = basins(q)

                k = k+1 
                if (mod(k,1000) .eq. 0) write(*,"(i10,a3,i10)") k, " / ", grid%npts 

            end do 
        end do 

        ! Write a mask of the original basin extent (no ocean points)
        outmask = 0 
        where (outvar .ne. missing_value) outmask = 1 
        call nc_write(filename,"basin_mask",outmask,dim1="xc",dim2="yc", &
                      units="1",missing_value=int(missing_value))
        call nc_write_attr(filename,"basin_mask","long_name","Mask of original basin extent")
!         call nc_write_attr(filename,"basin_mask","grid_mapping",trim(grid%mtype))
        call nc_write_attr(filename,"basin_mask","coordinates","lat2D lon2D")
            
        ! Fill in basins over ocean too
        call fill_nearest(outvar,missing_value)

        ! First write basins including sub basins if available
        if (trim(domain) .eq. "Greenland") then 
            call nc_write(filename,"basin_sub",real(outvar),dim1="xc",dim2="yc", &
                          units="1",missing_value=real(missing_value))
            call nc_write_attr(filename,"basin_sub","long_name","Basins and sub-basins")
!               call nc_write_attr(filename,"basin_sub","grid_mapping",trim(grid%mtype))
            call nc_write_attr(filename,"basin_sub","coordinates","lat2D lon2D")
            
        end if 

        ! Now whole number basins (aggregate basins)
        outvar = floor(outvar)
        call nc_write(filename,"basin",nint(outvar),dim1="xc",dim2="yc", &
                      units="1",missing_value=int(missing_value))
        call nc_write_attr(filename,"basin","long_name","Basins")
!         call nc_write_attr(filename,"basin","grid_mapping",trim(grid%mtype))
        call nc_write_attr(filename,"basin","coordinates","lat2D lon2D")
                    
        return 

    end subroutine nasaBasins_to_grid 

end module nasaBasins
