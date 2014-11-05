
! This program will take raw MAR output and 
! clean it to be useable as input to REMBOv2 and SEMIC
! *Eliminating unneeded dimensions
! *Cleaning days to range from 1 to 365 (eliminates 29Feb)

! To compile:
! gfortran -I/opt/local/include -o mar_rawclean.x ncio.f90 mar_rawclean.f90 -L/opt/local/lib -lnetcdff -lnetcdf
!

program rawclean

    use ncio 

    implicit none 

    integer :: nx, ny, nlev, nt
    double precision, allocatable, dimension(:) :: x, y, time
    character(len=512) :: filename0, filename
    character(len=256) :: xnm0, ynm0, xnm, ynm, tnm, gm

    real, allocatable :: var2D(:,:), var3D(:,:,:), var3D365(:,:,:)
    real, parameter :: mv0 = -1.d34 
    real, parameter :: mv  = -9999.0 

    integer :: argcount
    character(len=512) :: args(2)

    argcount = command_argument_count()
    if (argcount .ne. 2) then 
        write(*,*) "Two arguments should be provided to program: INPUTFOLDER INPUTFILE"
        stop 
    end if 
    CALL get_command_argument(1, args(1))
    CALL get_command_argument(2, args(2))

    filename0 = trim(args(1))//"/"//trim(args(2))
    filename  = trim(args(1))//"clean/"//trim(args(2))

    if (trim(filename0) == trim(filename)) then 
        write(*,*) "Error: new filename is the same as the input filename."
        write(*,*) "Input filename:  "//trim(filename0)
        write(*,*) "Output filename: "//trim(filename)
        stop
    end if 

    xnm0      = "X10_60"
    ynm0      = "Y17_109"

    xnm       = "x"
    ynm       = "y"
    tnm       = "time"
    gm        = "stereographic"

    nx = nc_size(filename0,xnm0)
    ny = nc_size(filename0,ynm0)
    nt = nc_size(filename0,"TIME")
    write(*,*) "nx,ny,nt = ",nx, ny, nt 

    allocate(x(nx),y(ny),time(366))
    allocate(var2D(nx,ny),var3D(nx,ny,nt),var3D365(nx,ny,365))

    call nc_read(filename0,xnm0,x)
    call nc_read(filename0,ynm0,y)

    call nc_create(filename)
    call nc_write_map(filename,gm,lambda=-39.d0,phi=71.d0,x_e=0.d0,y_n=0.d0)
    call nc_write_dim(filename,xnm,x=x)
    call nc_write_dim(filename,ynm,x=y)
    call nc_write_dim(filename,"time",x=1,dx=1,nx=365,units="day",calendar="360day")

    ! Longitude 
    call nc_read(filename0,"LON",var2D,missing_value=mv)
    call nc_write(filename,"LON",var2D,dim1=xnm,dim2=ynm,       &
                  long_name="Longitude",units="degrees_east",   &
                  missing_value=mv,grid_mapping=gm)

    ! Latitude 
    call nc_read(filename0,"LAT",var2D,missing_value=mv)
    call nc_write(filename,"LAT",var2D,dim1=xnm,dim2=ynm,       &
                  long_name="Latitude",units="degrees_north",   &
                  missing_value=mv,grid_mapping=gm)

    ! ====== Variables ======

    call nc_read(filename0,"SH",var2D,missing_value=mv)
    call nc_write(filename,"SH",var2D,dim1=xnm,dim2=ynm,        &
                  long_name="Surface Height",units="m",         &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SRF",var2D,missing_value=mv)
    call nc_write(filename,"SRF",var2D,dim1=xnm,dim2=ynm,       &
                  long_name="Surface Type",units="-",           &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SLO",var2D,missing_value=mv)
    call nc_write(filename,"SLO",var2D,dim1=xnm,dim2=ynm,       &
                  long_name="Surface Slope",units="-",          &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"MSK",var2D,missing_value=mv)
    call nc_write(filename,"MSK",var2D,dim1=xnm,dim2=ynm,       &
                  long_name="Ice Sheet Area",units="-",         &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SHSN0",var3D,start=[1,1,1,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"SHSN0",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Ini. old firn/ice thickness",units="m", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SHSN2",var3D,start=[1,1,1,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"SHSN2",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Snow Pack Height above Ice",units="m", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SHSN3",var3D,start=[1,1,2,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"SHSN3",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Snow Pack Height Total",units="m", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SMB",var3D,start=[1,1,1,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"SMB",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Surface Mass Balance (SMB~SF+RF-RU-SU-SW)",units="mmWE/day", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SU",var3D,start=[1,1,1,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"SU",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Sublimation and evaporation",units="mmWE/day", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"ME",var3D,start=[1,1,1,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"ME",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Meltwater production",units="mmWE/day", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"RZ",var3D,start=[1,1,1,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"RZ",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Meltwater refreezing and deposition",units="mmWE/day", &
                  missing_value=mv,grid_mapping=gm)

!     call nc_read(filename0,"SW",var3D,start=[1,1,1,1],count=[nx,ny,1,nt],missing_value=mv)
!     call del29feb(var3D,var3D365)
!     call nc_write(filename,"SW",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
!                   long_name="Surface Water",units="mmWE/day", &
!                   missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SF",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"SF",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Snowfall",units="mmWE/day", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"RF",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"RF",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Rainfall",units="mmWE/day", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"RU",var3D,start=[1,1,1,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"RU",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Run-off of meltwater and rain water",units="mmWE/day", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"UU",var3D,start=[1,1,23,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"UU",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="x-Wind Speed component",units="m/s", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"VV",var3D,start=[1,1,23,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"VV",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="y-Wind Speed component",units="m/s", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"TT",var3D,start=[1,1,23,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"TT",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Temperature",units="C", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"QQ",var3D,start=[1,1,23,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"QQ",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Specific Humidity",units="g/kg", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SP",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"SP",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Surface Pressure",units="hPa", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"RH",var3D,start=[1,1,3,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"RH",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Relative Humidity",units="%", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"TTMIN",var3D,start=[1,1,3,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"TTMIN",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Min Temp",units="C", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"TTMAX",var3D,start=[1,1,3,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"TTMAX",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Max Temp",units="C", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"UV",var3D,start=[1,1,3,1],count=[nx,ny,1,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"UV",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Wind Speed",units="m/s", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SWD",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"SWD",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Short Wave Downward",units="W/m2", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"LWD",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"LWD",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Long  Wave Downward",units="W/m2", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"LWU",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"LWU",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Long  Wave Upward",units="W/m2", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"SHF",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"SHF",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Sensible Heat Flux",units="W/m2", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"LHF",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"LHF",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Latent Heat Flux",units="W/m2", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"AL",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"AL",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Albedo",units="-", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"CC",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"CC",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Cloud Cover",units="-", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"ST",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"ST",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Surface Temperature",units="C", &
                  missing_value=mv,grid_mapping=gm)

    call nc_read(filename0,"PDD",var3D,start=[1,1,1],count=[nx,ny,nt],missing_value=mv)
    call del29feb(var3D,var3D365)
    call nc_write(filename,"PDD",var3D365,dim1=xnm,dim2=ynm,dim3=tnm, &
                  long_name="Postive Degree Day",units="C", &
                  missing_value=mv,grid_mapping=gm)


contains 

    subroutine del29feb(var,var365)

        implicit none 

        real, dimension(:,:,:) :: var, var365
        integer, parameter :: k28 = 59
        integer :: nt 

        nt = size(var,3)

        if (nt .eq. 366) then 
            var365(:,:,1:k28) = var(:,:,1:k28)
            var365(:,:,k28+1:365) = var(:,:,k28+2:366)
        else
            var365 = var 
        end if 

        return

    end subroutine del29feb 

end program rawclean 






