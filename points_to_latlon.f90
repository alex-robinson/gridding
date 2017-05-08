
program points_to_latlon

    use coordinates 

    implicit none 

    type(points_class) :: pts 
    character(len=128), allocatable :: regions(:)
    character(len=512) :: filename_in, filename_out 
    integer :: q 

    ! ===============================
    !
    ! Greenland regions 
    !
    ! ===============================

    if (.FALSE.) then 

    if (allocated(regions)) deallocate(regions)
    allocate(regions(5))

    regions(1) = "ellesmere"
    regions(2) = "grl-inner"
    regions(3) = "grl"
    regions(4) = "iceland"
    regions(5) = "svalbard" 

    do q = 1, size(regions)
        filename_in  = "regions/Bamber2013/polygon_"//trim(regions(q))//"_Bamber2013.txt"
        filename_out = "regions/polygon_"//trim(regions(q))//".txt"
        
        call points_init(pts,name="pts-B13",mtype="polar_stereographic",units="meters", &
                             filename=filename_in,skip=1,lon180=.TRUE.,lambda=-39.d0,phi=71.d0,alpha=19.0d0)
        
        call write_ascii(filename_out,pts%lon,pts%lat,xnm="lon",ynm="lat")
    end do 
    
    end if 

    ! ===============================
    !
    ! Antarctica regions 
    !
    ! ===============================
    
    if (.FALSE.) then 

    if (allocated(regions)) deallocate(regions)
    allocate(regions(2))
    
    regions(1) = "ant"
    regions(2) = "ant-inner"

    do q = 1, 2
        filename_in  = "regions/BEDMAP2/polygon_"//trim(regions(q))//"_BEDMAP2.txt"
        filename_out = "regions/polygon_"//trim(regions(q))//".txt"
        
        call points_init(pts,name="pts-BEDMAP2",mtype="polar_stereographic",units="km", &
                             filename=filename_in,skip=1,lon180=.TRUE.,lambda=0.d0,phi=-71.d0,alpha=19.0d0)
        
        call write_ascii(filename_out,pts%lon,pts%lat,xnm="lon",ynm="lat")
    end do 

    end if 

    ! ===============================
    !
    ! Greenland ice cores 
    !
    ! ===============================
    
    if (.FALSE.) then 



    end if 
    
contains 

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

end program points_to_latlon 

