module banderas2018
    
    use gridding_datasets
    use coord
    use ncio 
    use gaussian_filter
    
    implicit none 

    private 
    public :: banderas2018_to_grid

contains 


    subroutine banderas2018_to_grid(outfldr,grid,domain)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !    Banderas et al (2018) data files
        !
        ! =========================================================

        implicit none 

        character(len=*), intent(IN) :: outfldr
        type(grid_class), intent(IN) :: grid 
        character(len=*), intent(IN) :: domain

        ! Local variables 
        integer             :: max_neighbors 
        double precision    :: lat_lim 
        character(len=512)  :: filename 
        character(len=1024) :: desc, ref 

        type(grid_class)    :: grid0 
        character(len=1024) :: filename0  

        type input_type 
            integer :: nx, ny
            double precision, allocatable :: var2D(:,:) 
        end type 

        type(input_type) :: inp

        type(map_class) :: map 
        type(map_scrip_class) :: mps 
        double precision, allocatable :: outvar(:,:)
        integer, allocatable          :: outmask(:,:)
        character(len=12) :: xnm, ynm 
        character(len=56) :: xunits, yunits

        ! Define the output filename 
        write(filename,"(a)") trim(outfldr)//"/"//trim(grid%name)// &
                          "_Banderas2018-MIS3.nc"

        ! ==== Input grid and dataset information ==== 

        filename0 = "/Users/robinson/GoogleDriveUCM/wrk/mypapers/pmip-mis3/wrk/&
                    &data/Banderas2018/Banderas2018_MIS3.nc"

        ! Define input grid
        call grid_init(grid0,name="NH-40KM-B18",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=40.d0,nx=224,dy=40.d0,ny=208, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

        call grid_allocate(grid0,inp%var2D)

        ! Define the input filenames
        desc    = "Ice sheet model output for MIS-3 (index=260) in output file m3out.nc"
        ref     = "Banderas et al. (2018, gmd)"

        ! ==== MAPPING INFORMATION ====

if (.FALSE.) then 
        ! Define input grid in grid description file
        ! Note: this grid definition file was prepared by hand
        ! as climber-x output is not yet prepared to define 
        ! grid properly for cdo.
        call grid_write_cdo_desc_short(grid0,fldr="maps") 
        
        ! Define output grid in grid description file 
        call grid_write_cdo_desc_short(grid,fldr="maps") 
        
        ! Generate SCRIP interpolation weights 
        call map_scrip_init(mps,grid0%name,grid%name,fldr="maps",src_nc=filename0)

else 

        max_neighbors = 10 
        lat_lim       = 1.0d0 
        call map_init(map,grid0,grid,max_neighbors=max_neighbors,lat_lim=lat_lim,fldr="maps",load=.TRUE.)

end if 

        ! Initialize output variable arrays
        call grid_allocate(grid,outvar)
        call grid_allocate(grid,outmask)

        xnm = "xc"
        ynm = "yc" 
        xunits = "kilometers"
        yunits = "kilometers" 

        if (trim(domain) .eq. "Global") then 
            xnm = "lon"
            ynm = "lat"
            xunits = "degrees_east"
            yunits = "degrees_north" 
        end if 

        ! Initialize the output file
        call nc_create(filename)
        call nc_write_dim(filename, xnm,   x=grid%G%x,units=xunits)
        call nc_write_dim(filename, ynm,   x=grid%G%y,units=yunits)
        call grid_write(grid,filename,xnm=xnm,ynm=ynm,create=.FALSE.)
        
        ! Write meta data 
        call nc_write_attr(filename,"Description",desc)
        call nc_write_attr(filename,"Reference",ref)

        ! === Surface elevation ===

        ! Load reference topography in order to adjust temps to sea-level temps  
        call nc_read(filename0,"S",inp%var2D,missing_value=mv) 

        ! Map field
        outvar = mv 
        if (.not. trim(domain) .eq. "Global") then 
            call map_field_conservative_map1(map%map,"z_srf",inp%var2D,outvar, &
                                                            method="mean",missing_value=mv)
        else 
            call map_field(map,"z_srf",inp%var2D,outvar,outmask,"nn",fill=.FALSE., &
                                missing_value=mv,radius=100d3)
        end if 

        ! Write to file
        call nc_write(filename,"z_srf",real(outvar),dim1=xnm,dim2=ynm)

        ! Write variable metadata
        call nc_write_attr(filename,"z_srf","units","m")
        call nc_write_attr(filename,"z_srf","long_name","Surface elevation")
        call nc_write_attr(filename,"z_srf","coordinates","lat2D lon2D")
        
        ! === Ice thickness ===
        
        ! Load reference topography in order to adjust temps to sea-level temps  
        call nc_read(filename0,"H",inp%var2D,missing_value=mv) 
        where(inp%var2D .le. 1.0) inp%var2D = 0.0 

        ! Map field
        outvar = mv 

        if (.not. trim(domain) .eq. "Global") then 
            call map_field_conservative_map1(map%map,"H_ice",inp%var2D,outvar, &
                                                            method="mean",missing_value=mv)
        else 
            call map_field(map,"H_ice",inp%var2D,outvar,outmask,"nn",fill=.FALSE., &
                                missing_value=mv,radius=100d3)
        end if 

        ! Write to file
        call nc_write(filename,"H_ice",real(outvar),dim1=xnm,dim2=ynm)

        ! Write variable metadata
        call nc_write_attr(filename,"H_ice","units","m")
        call nc_write_attr(filename,"H_ice","long_name","Ice thickness")
        call nc_write_attr(filename,"H_ice","coordinates","lat2D lon2D")
        
        ! === Velocity ===
        
        ! Load reference topography in order to adjust temps to sea-level temps  
        call nc_read(filename0,"V",inp%var2D,missing_value=mv) 
        where(inp%var2D .le. 1.0) inp%var2D = 0.0 

        ! Map field
        outvar = mv 

        if (.not. trim(domain) .eq. "Global") then 
            call map_field_conservative_map1(map%map,"uxy",inp%var2D,outvar, &
                                                            method="mean",missing_value=mv)
        else 
            call map_field(map,"uxy",inp%var2D,outvar,outmask,"nn",fill=.FALSE., &
                                missing_value=mv,radius=100d3)
        end if 

        ! Write to file
        call nc_write(filename,"uxy",real(outvar),dim1=xnm,dim2=ynm)

        ! Write variable metadata
        call nc_write_attr(filename,"uxy","units","m")
        call nc_write_attr(filename,"uxy","long_name","Ice velocity")
        call nc_write_attr(filename,"uxy","coordinates","lat2D lon2D")
        

        return 

    end subroutine banderas2018_to_grid


end module banderas2018