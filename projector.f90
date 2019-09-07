program projector

    use coord 
    use control 
    
    implicit none 


    type(points_class) :: pts
    character(len=56)  :: mtype  
    real(8) :: lambda, phi 

    integer :: narg 
    character(len=56)  :: arg_lon, arg_lat 
    real(8) :: x, y, lon, lat 

    ! GRL domain 
    mtype  = "polar_stereographic"
    lambda = -45.d0 
    phi    =  70.d0 

    ! Initial points object with projection information 
    call points_init(pts,name="projector",mtype=trim(mtype),units="kilometers", &
                        lon180=.TRUE.,x=[0.0d0],y=[0.0d0],lambda=lambda,phi=phi,latlon=.FALSE.,verbose=.FALSE.)
    
    ! Test lon/lat points
    lon = -48.37d0
    lat =  58.21d0

    ! Get lon/lat from command line 
    narg = command_argument_count() 
    if (narg .gt. 0) then 
        call get_command_argument(1,arg_lon)
        call get_command_argument(2,arg_lat)

        read(arg_lon,*)  lon
        read(arg_lat,*)  lat
    end if 

    ! Project from lon/lat => x/y
    call oblimap_projection(lon,lat,x,y,pts%proj)
    x = x/pts%xy_conv
    y = y/pts%xy_conv

    write(*,"(a,2f10.2,a4,2f10.2)") "lonlat to xy: ", lon, lat, " >> ", x, y  


                


end program

