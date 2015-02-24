

module gridding_datasets

    implicit none 

    double precision, parameter :: missing_value = -9999.d0
    
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
                if (i1 .le. size(var1,1) .and. j1 .le. size(var1,2)) &
                    var1(i1,j1) = var(i,j)
            end do 
        end do 

        return
    end subroutine thin 

end module gridding_datasets