
program gentopo

    use ncio 
    use coordinates
    use gridding_datasets

    implicit none

    type(grid_class)   :: grid
    character(len=256) :: outfldr 
    
    write(*,*) 

    ! =========================================================
    !
    ! GRL-50KM Dataset
    !
    ! =========================================================

    if ( .TRUE. ) then 

        ! ## Define clim grid and output variable field ##
        call grid_init(grid,name="GRL-50KM",mtype="stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=50.d0,nx=37,dy=50.d0,ny=61, &
                       lambda=-40.d0,phi=72.d0,alpha=7.5d0)

        outfldr = "output/Greenland/"//trim(grid%name)

!         call Bamber13_to_grid(outfldr, grid, "Greenland",max_neighbors=20,lat_lim=2.d0)
!         call ecmwf_to_grid(   outfldr, grid, "GRL075",max_neighbors=8,lat_lim=2.d0)
!         call CERES_to_grid(   outfldr, grid, "Global",max_neighbors=8,lat_lim=2.d0)

!         call MARv33_to_grid(  outfldr, grid, "Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
!         call MARv33_to_grid(  outfldr, grid, "Greenland-MIROC5-RCP85",max_neighbors=20,lat_lim=2.d0)
!         call MARv32_to_grid(  outfldr, grid, "Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
        
        ! Climatologlies
!         call MARv33_to_grid(outfldr,grid,"Greenland-ERA",clim_range=[1981,2010])
        call MARv33_to_grid(outfldr,grid,"Greenland-MIROC5-RCP85",clim_range=[1981,2010])
        call MARv33_to_grid(outfldr,grid,"Greenland-MIROC5-RCP85",clim_range=[2071,2100])
        
    end if 

    ! =========================================================
    !
    ! GRL-20KM Dataset
    !
    ! =========================================================

    if ( .FALSE. ) then 
        
        ! ## Define ice grid and output variable field ##
        call grid_init(grid,name="GRL-20KM",mtype="stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=20.d0,nx=90,dy=20.d0,ny=150, &
                       lambda=-40.d0,phi=72.d0,alpha=7.5d0)

        outfldr = "output/Greenland/"//trim(grid%name)

        call Bamber13_to_grid(outfldr, grid, "Greenland",max_neighbors=20,lat_lim=2.d0)
        call ecmwf_to_grid(   outfldr, grid, "GRL075",max_neighbors=8,lat_lim=2.d0)
        call CERES_to_grid(   outfldr, grid, "Global",max_neighbors=8,lat_lim=2.d0)

        call MARv33_to_grid(  outfldr, grid, "Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
        call MARv33_to_grid(  outfldr, grid, "Greenland-MIROC5-RCP85",max_neighbors=20,lat_lim=2.d0)
        call MARv32_to_grid(  outfldr, grid, "Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
    
    end if 

    ! =========================================================
    !
    ! GRL-20KMb Dataset (Bamber et al., 2001 grid)
    !
    ! =========================================================

    if ( .FALSE. ) then 
        
        ! Define Bamber et al. 2001 20KM grid and input variable field
        call grid_init(grid,name="Bamber01-20KM",mtype="stereographic",units="kilometers", &
                       lon180=.TRUE.,x0=-800.d0,dx=20.d0,nx=76,y0=-3400.d0,dy=20.d0,ny=141, &
                       lambda=-39.d0,phi=90.d0,alpha=7.5d0)

        outfldr = "output/Greenland/"//trim(grid%name)

        call Bamber13_to_grid(outfldr, grid, "Greenland",max_neighbors=20,lat_lim=2.d0)
        call ecmwf_to_grid(   outfldr, grid, "GRL075",max_neighbors=8,lat_lim=2.d0)
        call CERES_to_grid(   outfldr, grid, "Global",max_neighbors=8,lat_lim=2.d0)

        call MARv33_to_grid(  outfldr, grid, "Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
        call MARv33_to_grid(  outfldr, grid, "Greenland-MIROC5-RCP85",max_neighbors=20,lat_lim=2.d0)
        call MARv32_to_grid(  outfldr, grid, "Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
    
    end if 

    ! =========================================================
    !
    ! GRL-10KM Dataset
    !
    ! =========================================================

    if ( .FALSE. ) then 

        ! ## Define ice grid and output variable field ##
        call grid_init(grid,name="GRL-10KM",mtype="stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=10.d0,nx=180,dy=10.d0,ny=300, &
                       lambda=-40.d0,phi=72.d0,alpha=7.5d0)

        outfldr = "output/Greenland/"//trim(grid%name)
        
        call Bamber13_to_grid(outfldr, grid, "Greenland",max_neighbors=20,lat_lim=2.d0)
        call ecmwf_to_grid(   outfldr, grid, "GRL075",   max_neighbors=8, lat_lim=2.d0)
        call CERES_to_grid(   outfldr, grid, "Global",   max_neighbors=8, lat_lim=2.d0)

        call MARv33_to_grid(  outfldr, grid, "Greenland-ERA",         max_neighbors=20,lat_lim=2.d0)
        call MARv33_to_grid(  outfldr, grid, "Greenland-MIROC5-RCP85",max_neighbors=20,lat_lim=2.d0)
        call MARv32_to_grid(  outfldr, grid, "Greenland-ERA",         max_neighbors=20,lat_lim=2.d0)
    
    end if 


    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gentopo

