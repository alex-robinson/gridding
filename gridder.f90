
program gridder

    use ncio 
    use coordinates
    use control 

    ! Datasets that can be gridded by this program (alphabetical)
    use AN1CRUST
    use bedmap2
    use CERES
    use climber2
    use climber3a 
    use davini2015
    use ECMWF 
    use ETOPO 
    use GeothermalHeatFlux 
    use GreenlandVelocity
    use grisli_g40
    use MAR 
    use NasaBasins 
    use RACMO2 
    use regions
    use Rignot13_BasalMelt  
    use sediments 
    use stratigraphy
    use topo_reconstructions 
    use topographies_grl 

    implicit none

    type(grid_class)   :: grid
    character(len=256) :: domain, grid_name, outfldr 
    character(len=256) :: path, subfldr  

    write(*,*) 

    ! =========================================================
    !
    ! USER DEFINITIONS
    !
    ! =========================================================
    
    domain    = "Greenland"
    grid_name = "GRL-5KM"
    outfldr   = "output/"//trim(domain)//"/"//trim(grid_name)

    ! =========================================================
    !
    ! OUTPUT GRID DEFINITION
    !
    ! =========================================================

    call domain_definition(grid,grid_name) 

    ! Write a regional mask 
!     call write_regions(outfldr,grid,domain)

    ! =========================================================
    !
    ! DATASET TO GRID CALCULATIONS
    !
    ! =========================================================

    ! == Global datasets - applicable to all domains ==

!     call CERES_to_grid(outfldr,     grid,"Global",max_neighbors=4,lat_lim=2.d0)
    call etopo1_to_grid(outfldr,    grid,"Global",max_neighbors=1,lat_lim=1.d0,grad_lim=0d0) !0.05d0)
!     call ICE6GC_to_grid(outfldr,    grid,"Global",max_neighbors=4,lat_lim=2.d0)
!     call ICE5G_to_grid(outfldr,     grid,"Global",max_neighbors=4,lat_lim=2.d0)
!     call sedLaske_to_grid(outfldr,  grid,"Global",max_neighbors=4,lat_lim=2.d0)
!     call ghfDavies_to_grid(outfldr, grid,"Global",max_neighbors=4,lat_lim=2.d0)
!     call ghfShapiro_to_grid(outfldr,grid,"Global",max_neighbors=4,lat_lim=2.d0)
    
!     call ecmwf_to_grid(outfldr,grid,sigma=30.d0,max_neighbors=1,lat_lim=2.d0)
!     call ecmwf_to_grid(outfldr,grid,clim_range=[1981,2010])
!     call ecmwf_to_grid(outfldr,grid,clim_range=[1979,1998])

!     call ecmwf_ocn_to_grid(outfldr,grid,sigma=30.d0,max_neighbors=4,lat_lim=2.d0)
!     call ecmwf_ocn_to_grid(outfldr,grid,clim_range=[1981,2010])
!     call ecmwf_ocn_to_grid(outfldr,grid,clim_range=[1979,1998])

!     call ecmwf40_to_grid(outfldr,grid,sigma=100.d0,max_neighbors=1,lat_lim=2.d0)
!     call ecmwf40_to_grid(outfldr,grid,clim_range=[1958,2001])
!     call ecmwf40_to_grid(outfldr,grid,clim_range=[1961,1990])

    if (trim(domain) .eq. "Antarctica") then 
        ! == Antarctica only datasets ==
        write(*,*) "Processing Antarctica..."
        
!         call An15litho_to_grid(outfldr, grid,"Antarctica", max_neighbors=5, lat_lim=1.0d0)
!         call bedmap2_to_grid(outfldr,   grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0,grad_lim=0.05d0)
!         call bedmap2vel_to_grid(outfldr,grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0)
!         call bedmap2acc_to_grid(outfldr,grid,"Antarctica",max_neighbors=20,lat_lim=0.5d0)
!         call ghfMaule_to_grid(outfldr,  grid,"Antarctica",max_neighbors=4,lat_lim=2.d0)
!         call nasaBasins_to_grid(outfldr,grid,"Antarctica")
!         call Rignot13_BasalMelt_to_grid(outfldr,grid,"Antarctica",max_neighbors=10,lat_lim=1.d0, &
!                                         fill=.TRUE.,sigma=100.d0)

!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",max_neighbors=20,lat_lim=0.5d0)
!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2000,2010])
!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2001,2030])
!         call RACMO2rot_to_grid( outfldr, grid, "Antarctica-A1B",clim_range=[2071,2100])

!         call RACMO23_to_grid( outfldr, grid, "ANT27",max_neighbors=20,lat_lim=0.5d0)
!         call RACMO23_to_grid( outfldr, grid, "ANT27",clim_range=[1981,2010])

        ! Note: Antartica-c20 doesn't work because some files only contain 239 months of
        !      data while they should all have 240 months (1980-1999)
        !      This dataset is not used for now...
        !call RACMO2rot_to_grid( outfldr, grid, "Antarctica-c20",max_neighbors=20,lat_lim=0.5d0)
        !call RACMO2rot_to_grid( outfldr, grid, "Antarctica-c20",clim_range=[1980,1999])
        
        ! Testing old grisli fields 
!         call g40_topo_to_grid(outfldr,subfldr,grid,domain,path_in,sigma,max_neighbors,lat_lim)

    end if 

    if (trim(domain) .eq. "Greenland") then 
        ! == Greenland only datasets ==
        write(*,*) "Processing Greenland..."

        call Bamber13_to_grid(outfldr,grid,"Greenland",max_neighbors=20,lat_lim=0.5d0,grad_lim=0d0) !0.05d0)   
!         call ghfMaule_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=2.d0)
!         call MacGregor15_to_grid(outfldr,grid,"Greenland",max_neighbors=20,lat_lim=0.5d0)
!         call Morlighem14_to_grid(outfldr,grid,"Greenland",max_neighbors=20,lat_lim=0.5d0,grad_lim=0d0) !0.05d0)
!         call MARv35_to_grid(outfldr,grid,"Greenland-ERA",max_neighbors=20,lat_lim=2.d0)
!         call MARv35_to_grid(outfldr,grid,"Greenland-ERA",clim_range=[1981,2010])
!         call nasaBasins_to_grid(outfldr,grid,"Greenland")
!         call LGMsimpson_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=1.d0)
            
!         call huy3_to_grid(outfldr,grid,"Greenland",max_neighbors=4,lat_lim=1.d0)
!         call grlvel_to_grid(outfldr,grid,"Greenland",max_neighbors=10,lat_lim=1.d0)

!         path = "data/Davini_GreenlandAMOC/"
!         call davini2015_to_grid(outfldr,"Davini2015",grid,domain="control", &
!                                 path_in=path,sigma=40.d0,max_neighbors=4,lat_lim=5.d0)
!         call davini2015_to_grid(outfldr,"Davini2015",grid,domain="bedrock", &
!                                 path_in=path,sigma=40.d0,max_neighbors=4,lat_lim=5.d0)
    
    end if 


    write(*,*)
    write(*,*) "Regridding program finished."
    write(*,*)

end program gridder

