
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

            ! GREENLAND DOMAINS =======================

            case("GRL-120KM")
                call grid_init(grid,name="GRL-120KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=120.d0,nx=16,dy=120.d0,ny=26, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-40KM")
                call grid_init(grid,name="GRL-40KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=40.d0,nx=46,dy=40.d0,ny=76, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-20KM")
                call grid_init(grid,name="GRL-20KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=20.d0,nx=91,dy=20.d0,ny=151, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-10KM")
                call grid_init(grid,name="GRL-10KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=10.d0,nx=181,dy=10.d0,ny=301, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-5KM")
                call grid_init(grid,name="GRL-5KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=5.d0,nx=361,dy=5.d0,ny=601, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)
            
            case("GRL-2KM")
                call grid_init(grid,name="GRL-2KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=2.d0,nx=901,dy=2.d0,ny=1501, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-1KM")
                call grid_init(grid,name="GRL-1KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=1.d0,nx=1801,dy=1.d0,ny=3001, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-500M")
                call grid_init(grid,name="GRL-500M",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=0.5d0,nx=3601,dy=0.5d0,ny=6001, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("Bamber01-20KM")
                call grid_init(grid,name="Bamber01-20KM",mtype="polar_stereographic",units="kilometers", &
                               lon180=.TRUE.,x0=-800.d0,dx=20.d0,nx=76,y0=-3400.d0,dy=20.d0,ny=141, &
                               lambda=-39.d0,phi=90.d0,alpha=19.d0)

            case("Bamber01-10KM")
                call grid_init(grid,name="Bamber01-10KM",mtype="polar_stereographic",units="kilometers", &
                               lon180=.TRUE.,x0=-800.d0,dx=10.d0,nx=151,y0=-3400.d0,dy=10.d0,ny=281, &
                               lambda=-39.d0,phi=90.d0,alpha=19.d0)

            ! ANTARCTICA DOMAINS ======================= 

            case("ANT-80KM")
                call grid_init(grid,name="ANT-80KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=80.d0,nx=79,dy=80.d0,ny=74,lambda=0.d0,phi=-71.d0)

            case("ANT-40KM")
                call grid_init(grid,name="ANT-40KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=40.d0,nx=157,dy=40.d0,ny=147,lambda=0.d0,phi=-71.d0)

            case("ANT-20KM")
                call grid_init(grid,name="ANT-20KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=20.d0,nx=313,dy=20.d0,ny=293,lambda=0.d0,phi=-71.d0)

            case("ANT-10KM")
                call grid_init(grid,name="ANT-10KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=10.d0,nx=625,dy=10.d0,ny=585,lambda=0.d0,phi=-71.d0)

            case("ANT-5KM")
                call grid_init(grid,name="ANT-5KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=5.d0,nx=1249,dy=5.d0,ny=1169,lambda=0.d0,phi=-71.d0)

            case("ANT-2KM")
                call grid_init(grid,name="ANT-2KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=2.d0,nx=3121,dy=2.d0,ny=2921,lambda=0.d0,phi=-71.d0)
            
            case("ANT-1KM")
                call grid_init(grid,name="ANT-1KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=1.d0,nx=6241,dy=1.d0,ny=5841,lambda=0.d0,phi=-71.d0)

            ! NORTH DOMAINS ======================= 

            case("NH-40KM")
                call grid_init(grid,name="NH-40KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=40.d0,nx=225,dy=40.d0,ny=211, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-20KM")
                call grid_init(grid,name="NH-20KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=20.d0,nx=449,dy=20.d0,ny=421, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-10KM")
                call grid_init(grid,name="NH-10KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=10.d0,nx=897,dy=10.d0,ny=841, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-5KM")
                call grid_init(grid,name="NH-5KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=5.d0,nx=1793,dy=5.d0,ny=1681, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-2KM")
                call grid_init(grid,name="NH-2KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=2.d0,nx=4481,dy=2.d0,ny=4201, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-1KM")
                call grid_init(grid,name="NH-1KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=1.d0,nx=8961,dy=1.d0,ny=8401, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-40KM-old")
                call grid_init(grid,name="NH-40KM-old",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=40.d0,nx=224,dy=40.d0,ny=208, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

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