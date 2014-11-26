
program gentopo

    use ncio 
    use coordinates
    use gridding_datasets

    implicit none

    type(grid_class)   :: grid
    character(len=256) :: outfldr 
    
    double precision, dimension(6801,6801) :: var2D 

    write(*,*) 

    ! =========================================================
    !
    ! ANT-20KM Dataset
    !
    ! =========================================================

    if ( .TRUE. ) then 

        ! ## Define output variable field ##
        call grid_init(grid,name="ANT-20KM",mtype="stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=20.d0,nx=281,dy=20.d0,ny=281, &
                       lambda=0.d0,phi=-90.d0,alpha=19.0d0)

        outfldr = "output/Antarctica/"//trim(grid%name)

!         call bedmap2_to_grid(outfldr,grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0)

!         call CERES_to_grid(outfldr,grid,"Global",max_neighbors=8, lat_lim=2.d0)
        
!         call ecmwf_to_grid(outfldr,grid,"ANT075",max_neighbors=8, lat_lim=2.d0)
!         call ecmwf_to_grid(outfldr,grid,"ANT075",clim_range=[1981,2010])
        
        call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",max_neighbors=20,lat_lim=0.5d0)
        call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2000,2010])
        call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2001,2030])
        call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2071,2100])

        ! Note: Antartica-c20 doesn't work because some files only contain 239 months of
        !       data while they should all have 240 months (1980-1999)
        !       This dataset is not used for now...
!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-c20",max_neighbors=20,lat_lim=0.5d0)
!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-c20",clim_range=[1980,1999])

    end if 

    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gentopo

