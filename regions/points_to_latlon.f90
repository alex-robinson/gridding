
program points_to_latlon

    use coordinates 

    implicit none 

    type(points_class) :: pts 


    call points_init(pts,name="pts-B13",mtype="polar_stereographic",units="kilometers", &
                         filename="regions/Bamber2013/polygon_grl_Bamber2013.txt",skip=1, &
                         lon180=.TRUE.,lambda=-39.d0,phi=71.d0,alpha=19.0d0)
        

end program points_to_latlon 

