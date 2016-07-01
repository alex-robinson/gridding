

module gridding_datasets

    implicit none 

    double precision, parameter :: missing_value = -9999.d0
    double precision, parameter :: mv            = missing_value
    
    type var_defs
        character(len=512) :: filename, filenames(20)
        character(len=256) :: nm_in, nm_out  
        character(len=256) :: units_in, units_out 
        character(len=512) :: long_name
        character(len=256) :: method
        logical :: mask, dimextra
        character(len=256) :: plev
        double precision   :: conv 
        logical            :: fill 
    end type 

contains

    !##############################################
    !
    ! General subroutines related to the module
    !
    !##############################################

    ! Define some variable info for later manipulation
    subroutine def_var_info(var,filename,nm_in,nm_out,units,long_name, &
                            method,mask,dimextra,conv,plev,fill,filenames)
        implicit none 

        type(var_defs) :: var 
        character(len=*) :: filename, nm_in, nm_out, units, long_name
        character(len=*), optional :: method 
        logical, optional :: mask, dimextra
        character(len=*), optional :: plev 
        character(len=*), optional :: filenames(:)
        double precision, optional :: conv 
        logical, optional :: fill 

        var%filename  = trim(filename)
        var%nm_in     = trim(nm_in)
        var%nm_out    = trim(nm_out)
        var%long_name = trim(long_name)
        var%units_in  = trim(units)
        var%units_out = trim(units)

        var%method = "shepard"
        if (present(method)) var%method = trim(method)

        var%mask = .FALSE. 
        if (present(mask)) var%mask = mask 

        var%dimextra = .FALSE.
        if (present(dimextra)) var%dimextra = dimextra 

        var%plev = "None"
        if (present(plev)) var%plev = trim(plev)

        var%filenames(:) = "None"
        if (present(filenames)) var%filenames = filenames

        var%conv = 1.d0 
        if (present(conv)) var%conv = conv 

        var%fill = .FALSE. 
        if (present(fill)) var%fill = fill 

        return 

    end subroutine def_var_info

    function read_vector(filename,n,col,skip) result(var)
        ! Read a column of data from an ascii file

        implicit none 

        character(len=*) :: filename 
        integer :: n, col, skip 
        real(8) :: var(n), tmp(50)
        character(len=10) :: tmpc
        integer :: i 

        open(16,file=trim(filename),status="old")
        do i = 1, skip
            read(16,*) tmpc 
        end do 

        do i = 1, n 
            read(16,*) tmp(1:col-1), var(i)
        end do 

        close(16)

        return
    end function read_vector

    ! Extract a thinner version of an input array
    ! (new array should be a multiple of input array)
    subroutine thin(var1,var,by,missing_value)
        implicit none

        double precision, dimension(:,:) :: var, var1 
        integer :: by 
        double precision, optional :: missing_value 
        integer :: i,j, nx, ny 
        integer :: i1, j1
        double precision :: missing_val 

        missing_val = -9999.0 
        if (present(missing_value)) missing_val = missing_value

        nx = size(var,1)
        ny = size(var,2) 

        var1 = missing_val

        i1 = 0
        do i = 1, nx, by 
            i1 = i1+1 
            j1 = 0 
            do j = 1, ny, by  
                j1 = j1 + 1 
                if (i1 .le. size(var1,1) .and. j1 .le. size(var1,2)) &
                    var1(i1,j1) = var(i,j)
            end do 
        end do 

        return
    end subroutine thin 

    ! Extract a lower resolution average version of an input array
    ! (new array should be a multiple of input array)
    subroutine thin_ave(var1,var,by,missing_value)
        implicit none

        double precision, dimension(:,:) :: var, var1 
        integer :: by 
        double precision, optional :: missing_value 
        integer :: i,j, nx, ny 
        integer :: i1, j1

        double precision, allocatable :: wts(:,:), wts_now(:,:)
        double precision :: missing_val 
        integer :: nxn, cnt 

        missing_val = -9999.0 
        if (present(missing_value)) missing_val = missing_value

        nx = size(var,1)
        ny = size(var,2) 

        ! Define weights for neighbor averaging 
        nxn = (by-1)/2
        allocate(wts(by,by))
        allocate(wts_now(by,by))

        wts = 1.0 
        do i = 1, by 
            do j = 1, by 
                wts(i,j) = sqrt((i-1-real(by-1)/2.d0)**2+(j-1-real(by-1)/2.d0)**2)
            end do 
        end do 
        wts = 1.0 / (wts**2.0)    ! Shephard's distance weighting 
        wts = wts / sum(wts) 

        var1 = missing_val
        cnt  = 0 

        i1 = 0
        do i = nxn+1, nx-nxn, by 
            i1 = i1+1 
            j1 = 0 
            do j = nxn+1, ny-nxn, by  
                j1 = j1 + 1 
                if (i1 .le. size(var1,1) .and. j1 .le. size(var1,2)) then 
                    wts_now = wts 
                    where (var(i-nxn:i+nxn,j-nxn:j+nxn) .eq. missing_val) wts_now = 0.d0
                    if (sum(wts_now) .gt. 0.d0) then 
                        wts_now = wts_now / sum(wts_now) 
                        var1(i1,j1) = sum(var(i-nxn:i+nxn,j-nxn:j+nxn)*wts_now)
                    else 
                        cnt = cnt+1
                    end if 
                end if 
            end do 
        end do 

        if (cnt .gt. 0) then 
            write(*,*) "cnt = ", cnt 
            stop 
        end if 

        return
    end subroutine thin_ave 

    subroutine replace(s,text,rep,outs)
        ! Adapted from FUNCTION Replace_Text:
        ! http://fortranwiki.org/fortran/show/String_Functions
        CHARACTER(len=*)           :: s,text,rep
        CHARACTER(len=*), optional :: outs
        INTEGER :: i, nt, nr

        character(len=LEN(s)+100) :: tmps  ! Temp string to hold output 

        tmps = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
        DO
           i = INDEX(tmps,text(:nt)) ; IF (i == 0) EXIT
           tmps = tmps(:i-1) // rep(:nr) // tmps(i+nt:)
        END DO
        
        if (present(outs)) then 
            outs = trim(tmps)
        else 
            s = trim(tmps)
        end if 

        return 
        
    end subroutine replace

end module gridding_datasets
