module rtopo 

    use gridding_datasets
    use coord
    use ncio 
    
    implicit none 

    private 
    public :: rtopo_latlon_to_grid 
    public :: rtopo_to_grid 

contains 

    subroutine rtopo_latlon_to_grid(outfldr,grid,domain)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       RTOPO DATA
        !       * These data files are enormous (3.5 Gb). Instead of
        !       * using mapping routines, that may crash due to memory
        !       * requirements, the first step will be a quick and dirty
        !       * nearest neighbor solution to a high resolution, local 
        !       * grid. 
        !
        ! =========================================================
        ! PATH: /p/projects/megarun/greenrise/datasets/rtopo-2/
        ! londim = 43201
        ! latdim = 21601
        ! lon, lat 
        ! 
        ! RTopo-2.0.1_30sec_bedrock_topography.nc : bedrock_topography
        ! long_name = "ocean bathymetry; surface topography of continents; bedrock topography under grounded or floating ice"
        ! RTopo-2.0.1_30sec_aux.nc : amask
        ! long_name = "ice ocean rock mask"
        ! RTopo-2.0.1_30sec_ice_base_topography.nc : ice_base_topography
        ! long_name = "ice base topography for the Antarctic and Greenland ice sheets / ice shelves (ice draft for ice shelves and floating glaciers; zero in absence of ice)"
        ! RTopo-2.0.1_30sec_surface_elevation.nc : surface_elevation
        ! long_name = "upper ice surface height for the Antarctic and Greenland ice sheets / ice shelves (bedrock topography for ice-free continent; zero for ocean)"


        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim, grad_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 




        return 

    end subroutine rtopo_latlon_to_grid


    subroutine rtopo_to_grid(outfldr,grid,domain)
        ! Convert the variables to the desired grid format and write to file
        ! =========================================================
        !
        !       RTOPO DATA
        !       * These data files are enormous (3.5 Gb). Instead of
        !       * using mapping routines, that may crash due to memory
        !       * requirements, the first step will be a quick and dirty
        !       * nearest neighbor solution to a high resolution, local 
        !       * grid. 
        !
        ! =========================================================
        ! PATH: /p/projects/megarun/greenrise/datasets/rtopo-2/
        ! londim = 43201
        ! latdim = 21601
        ! lon, lat 
        ! 
        ! RTopo-2.0.1_30sec_bedrock_topography.nc : bedrock_topography
        ! long_name = "ocean bathymetry; surface topography of continents; bedrock topography under grounded or floating ice"
        ! RTopo-2.0.1_30sec_aux.nc : amask
        ! long_name = "ice ocean rock mask"
        ! RTopo-2.0.1_30sec_ice_base_topography.nc : ice_base_topography
        ! long_name = "ice base topography for the Antarctic and Greenland ice sheets / ice shelves (ice draft for ice shelves and floating glaciers; zero in absence of ice)"
        ! RTopo-2.0.1_30sec_surface_elevation.nc : surface_elevation
        ! long_name = "upper ice surface height for the Antarctic and Greenland ice sheets / ice shelves (bedrock topography for ice-free continent; zero for ocean)"


        implicit none 

        character(len=*) :: domain, outfldr 
        type(grid_class) :: grid 
        integer :: max_neighbors 
        double precision :: lat_lim, grad_lim 
        character(len=512) :: filename 
        character(len=1024) :: desc, ref 




        return 

    end subroutine rtopo_to_grid


end module rtopo 
