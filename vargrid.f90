
module vargrid


    implicit none


    type var_defs
        character(len=256) :: filename
        character(len=256) :: nm_in, nm_out  
        character(len=256) :: units_in, units_out 
        character(len=256) :: method
        logical :: mask, dimextra
    end type 

    double precision, parameter :: missing_value = -9999.d0
    
contains

    ! Define some variable info for later manipulation
    subroutine def_var_info(var,filename,nm_in,nm_out,units,method,mask,dimextra)
        implicit none 

        type(var_defs) :: var 
        character(len=*) :: filename,nm_in,nm_out,units
        character(len=*), optional :: method 
        logical, optional :: mask, dimextra

        var%filename  = trim(filename)
        var%nm_in     = trim(nm_in)
        var%nm_out    = trim(nm_out)
        var%units_in  = trim(units)
        var%units_out = trim(units)

        var%method = "shepard"
        if (present(method)) var%method = trim(method)

        var%mask = .FALSE. 
        if (present(mask)) var%mask = mask 

        var%dimextra = .FALSE.
        if (present(dimextra)) var%dimextra = dimextra 

        return 

    end subroutine def_var_info

    ! Extract a thinner version of an input array
    ! (new array should be a multiple of input array)
    subroutine thin(var1,var,by)
        implicit none

        double precision, dimension(:,:) :: var, var1 
        integer :: by 
        integer :: i,j, nx, ny 
        integer :: i1, j1

        nx = size(var,1)
        ny = size(var,2) 

        var1 = missing_value 

        i1 = 0
        do i = 1, nx, by 
            i1 = i1+1 
            j1 = 0 
            do j = 1, ny, by  
                j1 = j1 + 1 
                var1(i1,j1) = var(i,j)
            end do 
        end do 

        return
    end subroutine thin 

    ! Fill in missing values of an array with neighbor averages
    ! or with a specified fill_value
    subroutine fill(var,missing_value,fill_value)
        implicit none 
        double precision, dimension(:,:) :: var 
        double precision :: missing_value 
        double precision, optional :: fill_value

        integer :: q, nx, ny, i, j 
        integer, parameter :: qmax = 50 ! Iterations 

        double precision, dimension (3,3) :: neighb, weight
        double precision :: wtot, mval 
        double precision, dimension(:,:), allocatable :: filled
        nx = size(var,1)
        ny = size(var,2) 

        allocate(filled(nx,ny))

        if (present(fill_value)) then
            where(var .eq. missing_value) var = fill_value 
        end if 

        do q = 1, qmax 

            filled = missing_value 

            do i = 2, nx-1 
                do j = 2, ny-1 
                    neighb = var(i-1:i+1,j-1:j+1)

                    weight = 0.d0 
                    where (neighb .ne. missing_value) weight = 1.d0
                    wtot = sum(weight)

                    if (wtot .gt. 0.d0) then 
                        mval = sum(neighb*weight)/wtot
                        where (neighb .eq. missing_value) neighb = mval 
                    end if 

                    filled(i-1:i+1,j-1:j+1) = neighb 

                end do 
            end do 

            where(filled .ne. missing_value) var = filled 

            write(*,*) q," : Missing values: ", count(var .eq. missing_value)
            if ( count(var .eq. missing_value) .eq. 0 ) exit 
        end do 

        return
    end subroutine fill 

end module vargrid

