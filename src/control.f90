
module control

    use nml 
    use coord
    use ncio 
    
    implicit none

    type control_param_class 
        character(len=256) :: name
        character(len=256) :: grid_name 
        character(len=512) :: path_ref 
        character(len=512) :: topo_ref 
        character(len=512) :: clim_ref 
        character(len=512) :: forcing 
        integer :: year_range(2), ref_range(2)
    end type 
    
    type timer_class
        real (8) :: dtime1, dtime2 
    end type 

contains

    subroutine domain_definition(grid,grid_name)

        implicit none 

        type(grid_class) :: grid 
        character(len=*) :: grid_name 

        select case(trim(grid_name))

            ! Note - all North projections now use the ESPG-3413
            ! polar stereographic projection with (lambda=-45.d0,phi=70.d0)
            ! Smaller Northern domains like Eurasia and Greenland use
            ! the same projection for consistency. 
            ! ESPG-3413 (lambda=-45.d0,phi=70.d0) is used for Greenland in 
            ! model intercomparison exercises, eg ISMIP6. 

            ! NORTH DOMAINS ======================= 

!             case("NH-40KM")
!                 call grid_init(grid,name="NH-40KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=-4900.d0,dx=40.0d0,nx=221,y0=-5400.d0,dy=40.0d0,ny=221, &
!                         lambda=-45.d0,phi=70.d0)
! !                 call grid_init(grid,name="NH-40KM",mtype="stereographic",units="kilometers", &
! !                                lon180=.TRUE.,dx=40.d0,nx=225,dy=40.d0,ny=211, &
! !                                lambda=-53.d0,phi=78.d0,alpha=32.7d0)
! !                 call grid_init(grid,name="NH-40KM",mtype="stereographic",units="kilometers", &
! !                                lon180=.TRUE.,dx=40.d0,nx=224,dy=40.d0,ny=208, &
! !                                lambda=-53.d0,phi=78.d0,alpha=32.7d0)

!             case("NH-20KM")
!                 call grid_init(grid,name="NH-20KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=-4900.d0,dx=20.0d0,nx=441,y0=-5400.d0,dy=20.0d0,ny=441, &
!                         lambda=-45.d0,phi=70.d0)
            
!             case("NH-10KM")
!                 call grid_init(grid,name="NH-10KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=-4900.d0,dx=10.0d0,nx=881,y0=-5400.d0,dy=10.0d0,ny=881, &
!                         lambda=-45.d0,phi=70.d0)
            
!             case("NH-5KM")
!                 call grid_init(grid,name="NH-5KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=-4900.d0,dx=5.0d0,nx=1761,y0=-5400.d0,dy=5.0d0,ny=1761, &
!                         lambda=-45.d0,phi=70.d0)
            

            case("NH-32KM")
                call grid_init(grid,name="NH-32KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-4900.d0,dx=32.0d0,nx=276,y0=-5400.d0,dy=32.0d0,ny=276, &
                        lambda=-45.d0,phi=70.d0)

            case("NH-16KM")
                call grid_init(grid,name="NH-16KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-4900.d0,dx=16.0d0,nx=551,y0=-5400.d0,dy=16.0d0,ny=551, &
                        lambda=-45.d0,phi=70.d0)

            case("NH-8KM")
                call grid_init(grid,name="NH-8KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-4900.d0,dx=8.0d0,nx=1101,y0=-5400.d0,dy=8.0d0,ny=1101, &
                        lambda=-45.d0,phi=70.d0)

            case("NH-4KM")
                call grid_init(grid,name="NH-4KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-4900.d0,dx=4.0d0,nx=2201,y0=-5400.d0,dy=4.0d0,ny=2201, &
                        lambda=-45.d0,phi=70.d0)

            ! EURASIA DOMAINS ======================= 

!             case("EIS-40KM")
!                 call grid_init(grid,name="EIS-40KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=380.d0,dx=40.0d0,nx=89,y0=-5000.d0,dy=40.0d0,ny=161, &
!                         lambda=-45.d0,phi=70.d0)
            
!             case("EIS-20KM")
!                 call grid_init(grid,name="EIS-20KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=380.d0,dx=20.0d0,nx=177,y0=-5000.d0,dy=20.0d0,ny=321, &
!                         lambda=-45.d0,phi=70.d0)
            
!             case("EIS-10KM")
!                 call grid_init(grid,name="EIS-10KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=380.d0,dx=10.0d0,nx=353,y0=-5000.d0,dy=10.0d0,ny=641, &
!                         lambda=-45.d0,phi=70.d0)
            
!             case("EIS-5KM")
!                 call grid_init(grid,name="EIS-5KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=380.d0,dx=5.0d0,nx=705,y0=-5000.d0,dy=5.0d0,ny=1281, &
!                         lambda=-45.d0,phi=70.d0)
            
            case("EIS-32KM")
                call grid_init(grid,name="EIS-32KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=380.d0,dx=32.0d0,nx=111,y0=-5000.d0,dy=32.0d0,ny=201, &
                        lambda=-45.d0,phi=70.d0)
            
            case("EIS-16KM")
                call grid_init(grid,name="EIS-16KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=380.d0,dx=16.0d0,nx=221,y0=-5000.d0,dy=16.0d0,ny=401, &
                        lambda=-45.d0,phi=70.d0)
            
            case("EIS-8KM")
                call grid_init(grid,name="EIS-8KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=380.d0,dx=8.0d0,nx=441,y0=-5000.d0,dy=8.0d0,ny=801, &
                        lambda=-45.d0,phi=70.d0)
            
            case("EIS-4KM")
                call grid_init(grid,name="EIS-4KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=380.d0,dx=4.0d0,nx=881,y0=-5000.d0,dy=4.0d0,ny=1601, &
                        lambda=-45.d0,phi=70.d0)
            
            ! LAURENTIDE DOMAINS ======================

            case("LIS-32KM")
                call grid_init(grid,name="LIS-32KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-5400.d0,dx=32.0d0,nx=165,y0=-5300.d0,dy=32.0d0,ny=245, &
                        lambda=-45.d0,phi=70.d0)
            
            case("LIS-20KM")
                call grid_init(grid,name="LIS-20KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-5400.d0,dx=20.0d0,nx=263,y0=-5300.d0,dy=20.0d0,ny=389, &
                        lambda=-45.d0,phi=70.d0)
            
            case("LIS-16KM")
                call grid_init(grid,name="LIS-16KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-5400.d0,dx=16.0d0,nx=329,y0=-5300.d0,dy=16.0d0,ny=489, &
                        lambda=-45.d0,phi=70.d0)
            
            case("LIS-8KM")
                call grid_init(grid,name="LIS-8KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-5400.d0,dx=8.0d0,nx=657,y0=-5300.d0,dy=8.0d0,ny=977, &
                        lambda=-45.d0,phi=70.d0)
            
            case("LIS-4KM")
                call grid_init(grid,name="LIS-4KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-5400.d0,dx=4.0d0,nx=1313,y0=-5300.d0,dy=4.0d0,ny=1953, &
                        lambda=-45.d0,phi=70.d0)
            
            ! GREENLAND DOMAINS =======================

!             case("GRL-80KM")
!                 call grid_init(grid,name="GRL-80KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=-720.d0,dx=80.0d0,nx=22,y0=-3450.d0,dy=80.0d0,ny=37, &
!                         lambda=-45.d0,phi=70.d0)
            
!             case("GRL-40KM")
!                 call grid_init(grid,name="GRL-40KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=-720.d0,dx=40.0d0,nx=43,y0=-3450.d0,dy=40.0d0,ny=73, &
!                         lambda=-45.d0,phi=70.d0)
            
!             case("GRL-20KM")
!                 call grid_init(grid,name="GRL-20KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=-720.d0,dx=20.0d0,nx=85,y0=-3450.d0,dy=20.0d0,ny=145, &
!                         lambda=-45.d0,phi=70.d0)
            
!             case("GRL-10KM")
!                 call grid_init(grid,name="GRL-10KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=-720.d0,dx=10.0d0,nx=169,y0=-3450.d0,dy=10.0d0,ny=289, &
!                         lambda=-45.d0,phi=70.d0)
            
!             case("GRL-5KM")
!                 call grid_init(grid,name="GRL-5KM",mtype="polar_stereographic",units="kilometers", &
!                         lon180=.TRUE.,x0=-720.d0,dx=5.0d0,nx=337,y0=-3450.d0,dy=5.0d0,ny=577, &
!                         lambda=-45.d0,phi=70.d0)
            
            case("GRL-32KM")
                call grid_init(grid,name="GRL-32KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-720.d0,dx=32.0d0,nx=54,y0=-3450.d0,dy=32.0d0,ny=91, &
                        lambda=-45.d0,phi=70.d0)
            
            case("GRL-16KM")
                call grid_init(grid,name="GRL-16KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-720.d0,dx=16.0d0,nx=106,y0=-3450.d0,dy=16.0d0,ny=181, &
                        lambda=-45.d0,phi=70.d0)
            
            case("GRL-8KM")
                call grid_init(grid,name="GRL-8KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-720.d0,dx=8.0d0,nx=211,y0=-3450.d0,dy=8.0d0,ny=361, &
                        lambda=-45.d0,phi=70.d0)
            
            case("GRL-2KM")
                call grid_init(grid,name="GRL-2KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-720.d0,dx=2.0d0,nx=841,y0=-3450.d0,dy=2.0d0,ny=1441, &
                        lambda=-45.d0,phi=70.d0)
            
            case("GRL-1KM")
                call grid_init(grid,name="GRL-1KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-720.d0,dx=1.0d0,nx=1681,y0=-3450.d0,dy=1.0d0,ny=2881, &
                        lambda=-45.d0,phi=70.d0)

            case("Bamber01-20KM")
                call grid_init(grid,name="Bamber01-20KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-800.d0,dx=20.d0,nx=76,y0=-3400.d0,dy=20.d0,ny=141, &
                        lambda=-39.d0,phi=90.d0)

            case("Bamber01-10KM")
                call grid_init(grid,name="Bamber01-10KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-800.d0,dx=10.d0,nx=151,y0=-3400.d0,dy=10.d0,ny=281, &
                        lambda=-39.d0,phi=90.d0)


            case("GRL-2KM-EXT")
                call grid_init(grid,name="GRL-2KM-EXT",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-820.d0,dx=2.0d0,nx=941,y0=-4000.d0,dy=2.0d0,ny=1716, &
                        lambda=-45.d0,phi=70.d0)
            
            ! ANTARCTICA DOMAINS ======================= 

!             case("ANT-80KM")
!                 call grid_init(grid,name="ANT-80KM",mtype="polar_stereographic",units="kilometers", &
!                        lon180=.TRUE.,dx=80.d0,nx=79,dy=80.d0,ny=74,lambda=0.d0,phi=-71.d0)

            case("ANT-40KM")
                call grid_init(grid,name="ANT-40KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=40.d0,nx=157,dy=40.d0,ny=147,lambda=0.d0,phi=-71.d0)

!             case("ANT-20KM")
!                 call grid_init(grid,name="ANT-20KM",mtype="polar_stereographic",units="kilometers", &
!                        lon180=.TRUE.,dx=20.d0,nx=313,dy=20.d0,ny=293,lambda=0.d0,phi=-71.d0)

!             case("ANT-10KM")
!                 call grid_init(grid,name="ANT-10KM",mtype="polar_stereographic",units="kilometers", &
!                        lon180=.TRUE.,dx=10.d0,nx=625,dy=10.d0,ny=585,lambda=0.d0,phi=-71.d0)

!             case("ANT-5KM")
!                 call grid_init(grid,name="ANT-5KM",mtype="polar_stereographic",units="kilometers", &
!                        lon180=.TRUE.,dx=5.d0,nx=1249,dy=5.d0,ny=1169,lambda=0.d0,phi=-71.d0)

!             case("ANT-2KM")
!                 call grid_init(grid,name="ANT-2KM",mtype="polar_stereographic",units="kilometers", &
!                        lon180=.TRUE.,dx=2.d0,nx=3121,dy=2.d0,ny=2921,lambda=0.d0,phi=-71.d0)
            
!             case("ANT-1KM")
!                 call grid_init(grid,name="ANT-1KM",mtype="polar_stereographic",units="kilometers", &
!                        lon180=.TRUE.,dx=1.d0,nx=6241,dy=1.d0,ny=5841,lambda=0.d0,phi=-71.d0)

            ! ISIMIP6 grid definitions 

            case("ANT-32KM")
                call grid_init(grid,name="ANT-32KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,x0=-3040.d0,dx=32.d0,nx=191,y0=-3040.d0,dy=32.d0,ny=191, &
                       lambda=0.d0,phi=-71.d0)

            case("ANT-16KM")
                call grid_init(grid,name="ANT-16KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,x0=-3040.d0,dx=16.d0,nx=381,y0=-3040.d0,dy=16.d0,ny=381, &
                       lambda=0.d0,phi=-71.d0)

            case("ANT-8KM")
                call grid_init(grid,name="ANT-8KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,x0=-3040.d0,dx=8.d0,nx=761,y0=-3040.d0,dy=8.d0,ny=761, &
                       lambda=0.d0,phi=-71.d0)
            
            case("ANT-4KM")
                call grid_init(grid,name="ANT-4KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,x0=-3040.d0,dx=4.d0,nx=1521,y0=-3040.d0,dy=4.d0,ny=1521, &
                       lambda=0.d0,phi=-71.d0)
            
            case("ANT-2KM")
                call grid_init(grid,name="ANT-2KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,x0=-3040.d0,dx=2.d0,nx=3041,y0=-3040.d0,dy=2.d0,ny=3041, &
                       lambda=0.d0,phi=-71.d0)
                        

            case DEFAULT
                write(*,*) "domain_definition:: error: grid name not recognized: "//trim(grid_name)
                stop 

        end select

        return 

    end subroutine domain_definition

    subroutine control_par_load(par,filename,init)

        type(control_param_class)   :: par 
        character(len=*)            :: filename 
        logical, optional           :: init 
        logical                     :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        call nml_read(filename,"rembo_domain","name",       par%name,       init=init_pars)
        call nml_read(filename,"rembo_domain","grid_name",  par%grid_name,  init=init_pars)
        call nml_read(filename,"rembo_domain","path_ref",   par%path_ref,   init=init_pars)
        call nml_read(filename,"rembo_domain","topo_ref",   par%topo_ref,   init=init_pars)
        call nml_read(filename,"rembo_domain","clim_ref",   par%clim_ref,   init=init_pars)
        call nml_read(filename,"rembo_domain","forcing",    par%forcing,    init=init_pars)

        call nml_read(filename,"rembo_domain","ref_range",  par%ref_range,  init=init_pars)
        call nml_read(filename,"rembo_domain","year_range", par%year_range, init=init_pars)

        return

    end subroutine control_par_load

    subroutine timer_step(timer,step)

        implicit none 

        type(timer_class) :: timer 
        integer :: step 

        if (step .eq. 1) then 
            call cpu_time(timer%dtime1)           ! get current time in seconds   
        else if (step .eq. 2) then 
            call cpu_time(timer%dtime2)           ! get current time in seconds   
        else
            write(*,*) "timer_step:: unknown step: ", step 
            stop 
        end if 

        return 

    end subroutine timer_step 

    subroutine timer_print(timer,units,label)

        implicit none 

        type(timer_class) :: timer
        character(len=1)  :: units 
        character(len=*), optional :: label 
        character(len=512) :: time_label 

        double precision :: dtime 

        dtime = timer%dtime2 - timer%dtime1 

        select case(units)

            case("s")
                dtime = dtime 

            case("m")
                dtime = dtime / 60.d0 

            case("h") 
                dtime = dtime / (60.d0*60.d0)

            case DEFAULT
                write(*,*) "timer_print:: error: units not recognized: "//units

        end select

        time_label = "Calculation time: "
        if (present(label)) time_label = trim(label)//": "

        write(*,"(a,f18.3,1x,a1)") trim(time_label), dtime, units

        return 

    end subroutine timer_print

end module control 



! ================ OLD PROJECTIONS ================ 

!             case("GRL-120KM")
!                 call grid_init(grid,name="GRL-120KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=120.d0,nx=16,dy=120.d0,ny=26, &
!                                lambda=-40.d0,phi=72.d0,alpha=8.4d0)

!             case("GRL-40KM")
!                 call grid_init(grid,name="GRL-40KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=40.d0,nx=46,dy=40.d0,ny=76, &
!                                lambda=-40.d0,phi=72.d0,alpha=8.4d0)

!             case("GRL-20KM")
!                 call grid_init(grid,name="GRL-20KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=20.d0,nx=91,dy=20.d0,ny=151, &
!                                lambda=-40.d0,phi=72.d0,alpha=8.4d0)

!             case("GRL-10KM")
!                 call grid_init(grid,name="GRL-10KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=10.d0,nx=181,dy=10.d0,ny=301, &
!                                lambda=-40.d0,phi=72.d0,alpha=8.4d0)

!             case("GRL-5KM")
!                 call grid_init(grid,name="GRL-5KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=5.d0,nx=361,dy=5.d0,ny=601, &
!                                lambda=-40.d0,phi=72.d0,alpha=8.4d0)
            
!             case("GRL-2KM")
!                 call grid_init(grid,name="GRL-2KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=2.d0,nx=901,dy=2.d0,ny=1501, &
!                                lambda=-40.d0,phi=72.d0,alpha=8.4d0)

!             case("GRL-1KM")
!                 call grid_init(grid,name="GRL-1KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=1.d0,nx=1801,dy=1.d0,ny=3001, &
!                                lambda=-40.d0,phi=72.d0,alpha=8.4d0)

!             case("GRL-500M")
!                 call grid_init(grid,name="GRL-500M",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=0.5d0,nx=3601,dy=0.5d0,ny=6001, &
!                                lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            ! Old projection, before going back to polar_stereographic for grid homogeneity among NH domains...
!             case("NH-40KM")
!                 call grid_init(grid,name="NH-40KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=40.d0,nx=225,dy=40.d0,ny=211, &
!                                lambda=-53.d0,phi=78.d0,alpha=32.7d0)

!             case("NH-20KM")
!                 call grid_init(grid,name="NH-20KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=20.d0,nx=449,dy=20.d0,ny=421, &
!                                lambda=-53.d0,phi=78.d0,alpha=32.7d0)

!             case("NH-10KM")
!                 call grid_init(grid,name="NH-10KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=10.d0,nx=897,dy=10.d0,ny=841, &
!                                lambda=-53.d0,phi=78.d0,alpha=32.7d0)

!             case("NH-5KM")
!                 call grid_init(grid,name="NH-5KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=5.d0,nx=1793,dy=5.d0,ny=1681, &
!                                lambda=-53.d0,phi=78.d0,alpha=32.7d0)

!             case("NH-2KM")
!                 call grid_init(grid,name="NH-2KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=2.d0,nx=4481,dy=2.d0,ny=4201, &
!                                lambda=-53.d0,phi=78.d0,alpha=32.7d0)

!             case("NH-1KM")
!                 call grid_init(grid,name="NH-1KM",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=1.d0,nx=8961,dy=1.d0,ny=8401, &
!                                lambda=-53.d0,phi=78.d0,alpha=32.7d0)

!             case("NH-40KM-old")
!                 call grid_init(grid,name="NH-40KM-old",mtype="stereographic",units="kilometers", &
!                                lon180=.TRUE.,dx=40.d0,nx=224,dy=40.d0,ny=208, &
!                                lambda=-53.d0,phi=78.d0,alpha=32.7d0)
