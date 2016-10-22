program gridder_help
    
    use coord 
    use generic 

    implicit none 

    type(grid_class)   :: grid0, grid1 
    character(len=512) :: outfldr, dataset, path_in
    character(len=56), allocatable :: vname(:), vname_int(:)
    integer :: thin_fac

!     ! Original Bamber grid (1KM)
!     call grid_init(grid0,name="B13-1KM",mtype="polar_stereographic", &
!                 units="kilometers",lon180=.TRUE., &
!                 x0=-1300.d0,dx=1.d0,nx=2501,y0=-3500.d0,dy=1.d0,ny=3001, &
!                 lambda=-39.d0,phi=71.d0,alpha=19.0d0)

!     call grid_init(grid1,name="GRL-1KM",mtype="stereographic",units="kilometers", &
!                    lon180=.TRUE.,dx=1.d0,nx=1801,dy=1.d0,ny=3001, &
!                    lambda=-40.d0,phi=72.d0,alpha=8.4d0)

!     thin_fac = 1 

    ! Original Bamber grid - thinned (3KM)
    call grid_init(grid0,name="B13-3KM",mtype="polar_stereographic", &
                units="kilometers",lon180=.TRUE., &
                x0=-1298.5d0,dx=3.d0,nx=1251,y0=-3498.5d0,dy=3.d0,ny=1501, &
                lambda=-39.d0,phi=71.d0,alpha=19.0d0)

    call grid_init(grid1,name="GRL-3KM",mtype="stereographic",units="kilometers", &
                   lon180=.TRUE.,dx=3.d0,nx=901,dy=3.d0,ny=1501, &
                   lambda=-40.d0,phi=72.d0,alpha=8.4d0)

    thin_fac = 3 

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



end program gridder_help 