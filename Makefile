.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
objdir = .obj
libdir = ..

# Command-line options at make call
env   ?= None      # options: manto,eolo,airaki,iplex
debug ?= 0 

ifeq ($(env),manto) ## env=manto

    ## IFORT OPTIONS ##
    FC  = ifort
    INC_NC  = -I/home/jalvarez/work/librairies/netcdflib/include
    LIB_NC  = -L/home/jalvarez/work/librairies/netcdflib/lib -lnetcdf
    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

    FLAGS    = -module $(objdir) -L$(objdir) $(INC_NC)
    LFLAGS   = $(LIB_NC) $(LIB_MKL)

    DFLAGS   = -vec-report0 -O2 -fp-model precise -i_dynamic 
    ifeq ($(debug), 1)
        DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0 -fp-model precise -i_dynamic 
    endif

else ifeq ($(env),eolo) ## env=eolo

    ## IFORT OPTIONS ##
    FC  = ifort
    INC_NC  = -I/home/fispalma22/work/librairies/netcdflib/include
    LIB_NC  = -L/home/fispalma22/work/librairies/netcdflib/lib -lnetcdf
    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

    FLAGS    = -module $(objdir) -L$(objdir) $(INC_NC)
    LFLAGS   = $(LIB_NC) $(LIB_MKL)

    DFLAGS   = -vec-report0 -O2 -fp-model precise
    ifeq ($(debug), 1)
        DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0 -fp-model precise
    endif

else ifeq ($(env),airaki) ## env=airaki

    ## GFORTRAN OPTIONS ##
    FC  = gfortran
    INC_NC  = -I/opt/local/include
    LIB_NC  = -L/opt/local/lib -lnetcdff -lnetcdf

    FLAGS  = -I$(objdir) -J$(objdir) $(INC_NC)
    LFLAGS = $(LIB_NC)

    DFLAGS = -O3
    ifeq ($(debug), 1)  # ,underflow
        DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all
    endif

else ifeq ($(env),iplex) ## env=iplex

    ## IFORT OPTIONS ##
    FC  = ifort
    INC_NC  = -I/home/robinson/apps/netcdf/netcdf/include
    LIB_NC  = -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf
    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

    FLAGS    = -module $(objdir) -L$(objdir) $(INC_NC)
    LFLAGS   = $(LIB_NC)

    DFLAGS   = -vec-report0 -O3
    ifeq ($(debug), 1)
        DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0
    endif

else 
    
    ## None ##
    FC = $(error "Define env")

endif

## Individual libraries or modules ##
$(objdir)/ncio.o: $(libdir)/ncio/ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/planet.o: $(libdir)/coord/planet.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/geodesic.o: $(libdir)/coord/geodesic.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/projection_oblimap2.o: $(libdir)/coord/projection_oblimap2.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/coordinates.o: $(libdir)/coord/coordinates.f90 $(objdir)/projection_oblimap2.o \
	                     $(objdir)/geodesic.o $(objdir)/planet.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(libdir)/coord/interp1D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp2D.o: $(libdir)/coord/interp2D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp_time.o: $(libdir)/coord/interp_time.f90 $(objdir)/interp1D.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/polygons.o: $(libdir)/coord/polygons.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/gridding_datasets.o: gridding_datasets.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

## Complete programs

# Program to test interpolations of CCSM3 data
GRL: $(objdir)/ncio.o $(objdir)/geodesic.o $(objdir)/planet.o \
	         $(objdir)/projection_oblimap2.o $(objdir)/coordinates.o  \
	         $(objdir)/interp1D.o $(objdir)/interp_time.o $(objdir)/interp2D.o \
	         $(objdir)/gridding_datasets.o
	$(FC) $(DFLAGS) $(FLAGS) -o gentopo_GRL.x $^ gentopo_GRL.f90 $(LFLAGS)
	@echo " "
	@echo "    gentopo_GRL.x is ready."
	@echo " "

ANT: $(objdir)/ncio.o $(objdir)/geodesic.o $(objdir)/planet.o \
	         $(objdir)/projection_oblimap2.o $(objdir)/coordinates.o  \
	         $(objdir)/interp1D.o $(objdir)/interp_time.o $(objdir)/interp2D.o \
	         $(objdir)/gridding_datasets.o
	$(FC) $(DFLAGS) $(FLAGS) -o gentopo_ANT.x $^ gentopo_ANT.f90 $(LFLAGS)
	@echo " "
	@echo "    gentopo_ANT.x is ready."
	@echo " "

marclean: $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -o mar_rawclean.x $^ mar_rawclean.f90 $(LFLAGS)
	@echo " "
	@echo "    mar_rawclean.x is ready."
	@echo " "

marmonthly: $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -o mar_calcmonthly.x $^ mar_calcmonthly.f90 $(LFLAGS)
	@echo " "
	@echo "    mar_calcmonthly.x is ready."
	@echo " "

bedmap: $(objdir)/ncio.o $(objdir)/geodesic.o $(objdir)/planet.o \
	         $(objdir)/projection_oblimap2.o $(objdir)/coordinates.o  \
	         $(objdir)/interp1D.o $(objdir)/interp_time.o $(objdir)/interp2D.o
	$(FC) $(DFLAGS) $(FLAGS) -o bedmap2_netcdf.x $^ bedmap2_netcdf.f90 $(LFLAGS)
	@echo " "
	@echo "    bedmap2_netcdf.x is ready."
	@echo " "

nasa-basins: $(objdir)/polygons.o
	$(FC) $(DFLAGS) $(FLAGS) -o nasa_basins.x $^ nasa_basins.f90 $(LFLAGS)
	@echo " "
	@echo "    nasa_basins.x is ready."
	@echo " "

test: $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -o test.x $^ test.f90 $(LFLAGS)
	@echo " "
	@echo "    test.x is ready."
	@echo " "

clean:
	rm -f gentopo_grl.x gentopo_ant.x mar_rawclean.x $(objdir)/*.o $(objdir)/*.mod

# cleanall: cleansico cleanrembo cleansicoX
