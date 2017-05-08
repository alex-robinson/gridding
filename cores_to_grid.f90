
program cores_to_grid

    use coordinates 

    implicit none 

    type(points_class) :: pts 
    character(len=512) :: filename_in, filename_out 
    integer :: q 

    real(4), allocatable :: lat(:), lon(:) 
    character(len=56) :: grid_name 

    ! ===============================
    !
    ! Greenland ice cores 
    !
    ! ===============================
    
    if (.TRUE.) then 

    allocate(lat(7),lon(7))
    lon = [-37.6,-38.5,-42.3,-61.1,-43.8,-51.06,-45.0]
    lat = [ 72.6, 72.6, 75.1, 77.2, 65.2, 77.45, 76.0] 

    call points_init(pts,name="pts-GRL",mtype="stereographic",units="km", &
                         x=dble(lon),y=dble(lat),latlon=.TRUE.,lon180=.TRUE., &
                         lambda=-40.d0,phi=72.d0,alpha=8.4d0)
    
    filename_out = "regions/cores_greenland_"//trim(pts%name)//".txt"
    call write_ascii(filename_out,pts%x,pts%y,xnm="x",ynm="y")

    end if 
    
contains 

    subroutine read_ascii_column(filename,x,y,xnm,ynm)

        implicit none 

        character(len=*), intent(IN) :: filename 
        double precision, intent(IN) :: x(:), y(:) 
        character(len=*), intent(IN) :: xnm, ynm 
        integer :: i 
        open(7,file=trim(filename))

        write(7,"(2a14)") trim(xnm), trim(ynm)

        do i = 1, size(x)
            write(7,"(2f14.3)") x(i), y(i)
        end do 

        close(7)

        return 

    end subroutine read_ascii_column

    subroutine write_ascii(filename,x,y,xnm,ynm)

        implicit none 

        character(len=*), intent(IN) :: filename 
        double precision, intent(IN) :: x(:), y(:) 
        character(len=*), intent(IN) :: xnm, ynm 
        integer :: i 
        open(7,file=trim(filename))

        write(7,"(2a14)") trim(xnm), trim(ynm)

        do i = 1, size(x)
            write(7,"(2f14.3)") x(i), y(i)
        end do 

        close(7)

        return 

    end subroutine write_ascii 

end program cores_to_grid 

