program gridder_help
    ! This program will help get high resolution datasets onto the same projection
    ! as the basic ice_data grids. That way conservative interpolation is possible and easy.

    use coord 
    use generic 

    implicit none 

    type(grid_class)   :: grid0, grid1 
    character(len=512) :: outfldr, dataset, path_in
    character(len=56), allocatable :: vname(:), vname_int(:)
    integer :: thin_fac

    character(len=512), allocatable :: datasets(:)
    integer :: q 

    ! ====================================================
    !
    ! Bamber et al. 2013 - GREENLAND 
    !
    ! ====================================================
    
    if (.FALSE.) then 

!         ! Original Bamber grid (1KM)
!         call grid_init(grid0,name="B13-1KM",mtype="polar_stereographic", &
!                     units="kilometers",lon180=.TRUE., &
!                     x0=-1300.d0,dx=1.d0,nx=2501,y0=-3500.d0,dy=1.d0,ny=3001, &
!                     lambda=-39.d0,phi=71.d0,alpha=19.0d0)

!         call grid_init(grid1,name="GRL-1KM",mtype="stereographic",units="kilometers", &
!                        lon180=.TRUE.,dx=1.d0,nx=1801,dy=1.d0,ny=3001, &
!                        lambda=-40.d0,phi=72.d0,alpha=8.4d0)

!         thin_fac = 1 

        ! Original Bamber grid - thinned (2KM)
        call grid_init(grid0,name="B13-2KM",mtype="polar_stereographic", &
                    units="kilometers",lon180=.TRUE., &
                    x0=-1299.0d0,dx=2.d0,nx=1251,y0=-3499.0d0,dy=2.d0,ny=1501, &
                    lambda=-39.d0,phi=71.d0,alpha=19.0d0)

        call grid_init(grid1,name="GRL-2KM",mtype="stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=2.d0,nx=901,dy=2.d0,ny=1501, &
                       lambda=-40.d0,phi=72.d0,alpha=8.4d0)

        thin_fac = 2 

        outfldr = "output/Greenland"
        dataset = "TOPO-B13"
        path_in = "/data/sicopolis/data/Greenland/Greenland_bedrock_topography_V3.nc"

        allocate(vname(7),vname_int(1))
        vname(1)  = "BedrockElevation"
        vname(2)  = "SurfaceElevation"
        vname(3)  = "IceThickness"
        vname(4)  = "SurfaceRMSE"
        vname(5)  = "BedrockError"
        vname(6)  = "LandMask"
        vname(7)  = "Geoid"

        vname_int(1) = "LandMask"

        call generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)


    end if 

    ! ====================================================
    !
    ! BEDMAP2 - ANTARCTICA 
    !
    ! ====================================================

    if (.FALSE.) then 

        ! Original BEDMAP2 grid - thinned (10KM)
        call grid_init(grid0,name="BEDMAP2-10KM",mtype="polar_stereographic",units="kilometers",lon180=.TRUE., &
               x0=-3328.d0,dx=10.d0,nx=666,y0=-3328.d0,dy=10.d0,ny=666,lambda=0.d0,phi=-71.d0)

        call grid_init(grid1,name="ANT-10KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=10.d0,nx=625,dy=10.d0,ny=585,lambda=0.d0,phi=-71.d0)

        ! ==== TOPOGRAPHY (BEDMAP2) =========
        outfldr = "output/Antarctica"
        dataset = "TOPO-BEDMAP2"
        path_in = "/data/sicopolis/data/gridding_output/Antarctica/BEDMAP2-netcdf/BEDMAP2-1KM_BEDMAP2_topo.nc"

        allocate(vname(4),vname_int(1))
        vname(1)     = "zs"
        vname(2)     = "zb"
        vname(3)     = "H"
        vname(4)     = "mask_ice"

        vname_int(1) = "mask_ice"

        call generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)

        ! ==== VELOCITY (Rignot et al., 2011) =========
        outfldr = "output/Antarctica"
        dataset = "VEL-R11"
        path_in = "/data/sicopolis/data/gridding_output/Antarctica/BEDMAP2-netcdf/BEDMAP2-1KM_BEDMAP2_vel.nc"

        if (allocated(vname)) deallocate(vname)
        if (allocated(vname_int)) deallocate(vname_int)
        allocate(vname(3),vname_int(1))
        vname(1)     = "u"
        vname(2)     = "v"
        vname(3)     = "err"
        vname_int(1) = "None" 

        call generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)

        ! ==== ACCUMULATION (Rignot et al., 2011) =========
        outfldr = "output/Antarctica"
        dataset = "ACC-A06"
        path_in = "/data/sicopolis/data/gridding_output/Antarctica/BEDMAP2-netcdf/BEDMAP2-1KM_BEDMAP2_acc.nc"

        if (allocated(vname)) deallocate(vname)
        if (allocated(vname_int)) deallocate(vname_int)
        allocate(vname(1),vname_int(1))
        vname(1)     = "accum"
        vname_int(1) = "None"

        call generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)

    end if 



    ! ====================================================
    !
    ! Climber3 - NORTH
    ! Banderas2015 datafiles using old NH-40KM domain
    ! translated to new NH-40KM domain 
    !
    ! ====================================================
    
    if (.FALSE.) then 

        ! Original NH-40KM grid
        call grid_init(grid0,name="NH-40KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=40.d0,nx=224,dy=40.d0,ny=208, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

        call grid_init(grid1,name="NH-40KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=40.d0,nx=225,dy=40.d0,ny=211, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

        thin_fac = 1 

        outfldr = "output/North"

        ! == ATMOSPHERIC DATASETS == 

        allocate(datasets(6))
        datasets(1) = "holo_c3a_atm_ice5g_corr"
        datasets(2) = "holo_c3a_atm"
        datasets(3) = "interstadial_c3a_atm_ice5g_corr"
        datasets(4) = "interstadial_c3a_atm"
        datasets(5) = "stadial_c3a_atm_ice5g_corr"
        datasets(6) = "stadial_c3a_atm"
        
        allocate(vname(5),vname_int(1))
        vname(1)  = "zs"
        vname(2)  = "mask_land"
        vname(3)  = "t2m_ann"
        vname(4)  = "t2m_sum"
        vname(5)  = "pr_ann"

        vname_int(1) = "mask_land"

        do q = 1, size(datasets)
            dataset = datasets(q)
            path_in = "output/North/NH-40KM_old/Banderas2015/NH-40KM_"//trim(dataset)//".nc"
            call generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)

        end do 

        ! == OCEAN DATASETS == 

        deallocate(datasets)
        allocate(datasets(3))
        datasets(1) = "holo_c3a_ocn"
        datasets(2) = "interstadial_c3a_ocn"
        datasets(3) = "stadial_c3a_ocn"

        deallocate(vname,vname_int)
        allocate(vname(2),vname_int(1))
        vname(1)  = "to"
        vname(2)  = "mask_ocn"

        vname_int(1) = "mask_ocn"

        do q = 1, size(datasets)
            dataset = datasets(q)
            path_in = "output/North/NH-40KM_old/Banderas2015/NH-40KM_"//trim(dataset)//".nc"
            call generic_to_grid_3D_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)

        end do 

    end if 


    ! ====================================================
    !
    ! BAMBER 2013 - GREENLAND
    ! *To get old ice_data into new domain quickly
    !
    ! ====================================================

    if (.FALSE.) then 

        ! Original grid
        call grid_init(grid0,name="GRL-20KM_old",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=20.d0,nx=90,dy=20.d0,ny=150, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

        ! New grid
        call grid_init(grid1,name="GRL-20KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=20.d0,nx=91,dy=20.d0,ny=151, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

        thin_fac = 1 

        ! ==== TOPOGRAPHY (BEDMAP2) =========
        outfldr = "output/Greenland/GRL-20KM"
        dataset = "TOPO-M14"
        path_in = "/data/sicopolis/data/gridding_output/Greenland_old/GRL-20KM/GRL-20KM_TOPO-M14.nc"

        allocate(vname(4),vname_int(1))
        vname(1)     = "zs"
        vname(2)     = "zb"
        vname(3)     = "H"
        vname(4)     = "mask"

        vname_int(1) = "mask"

        call generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)

    end if 

    ! ====================================================
    !
    ! BEDMAP2 - ANTARCTICA
    ! *To get old ice_data into new domain quickly
    !
    ! ====================================================

    if (.FALSE.) then 

        ! Original grid
        call grid_init(grid0,name="ANT-40KM_old",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=40.d0,nx=156,dy=40.d0,ny=146,lambda=0.d0,phi=-71.d0)

        ! New grid
        call grid_init(grid1,name="ANT-40KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=40.d0,nx=157,dy=40.d0,ny=147,lambda=0.d0,phi=-71.d0)

        thin_fac = 1 

        ! ==== TOPOGRAPHY (BEDMAP2) =========
        outfldr = "output/Antarctica/ANT-40KM"
        dataset = "BEDMAP2"
        path_in = "/data/sicopolis/data/gridding_output/Antarctica_old/ANT-40KM/ANT-40KM_TOPO-BEDMAP2.nc"

        allocate(vname(5),vname_int(2))
        vname(1)     = "zs"
        vname(2)     = "zb"
        vname(3)     = "H"
        vname(4)     = "mask_ice"
        vname(5)     = "mask"

        vname_int(1) = "mask_ice"
        vname_int(2) = "mask"

        call generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)

    end if 

    if (.FALSE.) then 

        ! Original grid
        call grid_init(grid0,name="ANT-40KM_old",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=40.d0,nx=156,dy=40.d0,ny=146,lambda=0.d0,phi=-71.d0)

        ! New grid
        call grid_init(grid1,name="ANT-40KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=40.d0,nx=157,dy=40.d0,ny=147,lambda=0.d0,phi=-71.d0)

        thin_fac = 1 

        ! ==== TOPOGRAPHY (BEDMAP2) =========
        outfldr = "output/Antarctica/ANT-40KM"
        dataset = "VEL-R11"
        path_in = "/data/sicopolis/data/gridding_output/Antarctica_old/ANT-40KM/ANT-40KM_VEL-R11.nc"

        allocate(vname(3),vname_int(1))
        vname(1)     = "u"
        vname(2)     = "v"
        vname(3)     = "uv"

        vname_int(1) = "None"

        call generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)

    end if 

    if (.TRUE.) then 

        ! Original grid
        call grid_init(grid0,name="ANT-40KM_old",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=40.d0,nx=156,dy=40.d0,ny=146,lambda=0.d0,phi=-71.d0)

        ! New grid
        call grid_init(grid1,name="ANT-40KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=40.d0,nx=157,dy=40.d0,ny=147,lambda=0.d0,phi=-71.d0)

        thin_fac = 1 

        ! ==== TOPOGRAPHY (BEDMAP2) =========
        outfldr = "output/Antarctica/ANT-40KM"
        dataset = "ACC-A06"
        path_in = "/data/sicopolis/data/gridding_output/Antarctica_old/ANT-40KM/ANT-40KM_ACC-A06.nc"

        allocate(vname(1),vname_int(1))
        vname(1)     = "acc"

        vname_int(1) = "None"

        call generic_to_grid_nn(grid0,grid1,outfldr,dataset,path_in,vname,vname_int,thin_fac=thin_fac)

    end if 

end program gridder_help 