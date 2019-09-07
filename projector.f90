program projector

    use coord 
    use control 
    
    implicit none 


    type(points_class) :: pts
    character(len=56)  :: mtype  
    real(8) :: lambda, phi 

    real(8) :: x, y, lon, lat 

    write(*,*) 


    ! GRL domain 
    mtype  = "polar_stereographic"
    lambda = -45.d0 
    phi    =  70.d0 

    ! Initial points object with projection information 
    call points_init(pts,name="projector",mtype=trim(mtype),units="kilometers", &
                        lon180=.TRUE.,x=[0.0d0],y=[0.0d0],lambda=lambda,phi=phi,latlon=.FALSE.)
    
    lon = -48.37d0
    lat =  58.21d0

    ! Project from lon/lat => x/y
    call oblimap_projection(lon,lat,x,y,pts%proj)
    x = x/pts%xy_conv
    y = y/pts%xy_conv

    write(*,"(a,2f10.2,a4,2f10.2)") "lon/lat => x/y: ", lon, lat, " == ", x, y  


                


end program

