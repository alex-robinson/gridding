.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
objdir = .obj
libdir = libs
srcdir = src_datasets

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

#    ## IFORT OPTIONS ##
#    FC  = ifort
#    INC_NC  = -I/home/fispalma22/work/librairies/netcdflib/include
#    LIB_NC  = -L/home/fispalma22/work/librairies/netcdflib/lib -lnetcdf
#    INC_COORD = -I/home/fispalma25/robinson/models/EURICE/coord/.obj
#    LIB_COORD = /home/fispalma25/robinson/models/EURICE/coord/libcoordinates.a
#    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#
#    FLAGS    = -module $(objdir) -L$(objdir) $(INC_COORD) $(INC_NC)
#    LFLAGS   = $(LIB_COORD) $(LIB_NC) $(LIB_MKL)
#
#    DFLAGS   = -vec-report0 -O2 -fp-model precise
#    ifeq ($(debug), 1)
#        DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0 -fp-model precise
#    endif

    ## GFORTRAN OPTIONS ##
    FC  = gfortran
    INC_NC  = -I/home/fispalma25/apps/netcdf/netcdf/include
    LIB_NC  = -L/home/fispalma25/apps/netcdf/netcdf/lib -lnetcdff -lnetcdf
    INC_COORD = -I/home/fispalma25/apps/coordinates/.obj
    LIB_COORD = /home/fispalma25/apps/coordinates/libcoordinates.a

    FLAGS  = -I$(objdir) -J$(objdir) $(INC_COORD) $(INC_NC)
    LFLAGS = $(LIB_COORD) $(LIB_NC)

    DFLAGS = -O3
    ifeq ($(debug), 1)  # ,underflow
        DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all
    endif

else ifeq ($(env),airaki) ## env=airaki

    ## GFORTRAN OPTIONS ##
    FC  = gfortran
    INC_NC  = -I/opt/local/include
    LIB_NC  = -L/opt/local/lib -lnetcdff -lnetcdf
    INC_COORD = -I/Users/robinson/models/EURICE/coord/.obj
    LIB_COORD = /Users/robinson/models/EURICE/coord/libcoordinates.a

    FLAGS  = -I$(objdir) -J$(objdir) $(INC_COORD) $(INC_NC) 
    LFLAGS = $(LIB_COORD) $(LIB_NC)

    DFLAGS = -O3
    ifeq ($(debug), 1)  # ,underflow
        DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all
    endif

else ifeq ($(env),pik) ## env=pik

    ## IFORT OPTIONS ##
    FC  = ifort
    INC_NC  = -I${NETCDF_FORTRANROOT}/include
    LIB_NC  = -L${NETCDF_FORTRANROOT}/lib -lnetcdff -L${NETCDF_CROOT}/lib -lnetcdf
    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
    INC_COORD = -I/p/projects/tumble/robinson/EURICE/coord/.obj
	LIB_COORD = /p/projects/tumble/robinson/EURICE/coord/libcoordinates.a

    FLAGS    = -module $(objdir) -L$(objdir) $(INC_COORD) $(INC_NC) 
    LFLAGS   = $(LIB_COORD) $(LIB_NC)

    DFLAGS   = -O3
    ifeq ($(debug), 1)
        DFLAGS   = -C -g -traceback -ftrapuv -fpe0 -check all
    endif

else 
    
    ## None ##
    FC = $(error "Define env")

endif

## Individual libraries or modules ##
$(objdir)/nml.o: $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/control.o: $(srcdir)/control.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/gridding_datasets.o: $(srcdir)/gridding_datasets.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/regions.o: $(srcdir)/regions.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/generic.o: $(srcdir)/generic.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/AN1CRUST.o: $(srcdir)/AN1CRUST.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/bedmap2.o: $(srcdir)/bedmap2.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/CERES.o: $(srcdir)/CERES.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/climber2.o: $(srcdir)/climber2.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/climber3a.o: $(srcdir)/climber3a.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/davini2015.o: $(srcdir)/davini2015.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/ECMWF.o: $(srcdir)/ECMWF.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/ETOPO.o: $(srcdir)/ETOPO.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/GeothermalHeatFlux.o: $(srcdir)/GeothermalHeatFlux.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/GreenlandVelocity.o: $(srcdir)/GreenlandVelocity.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/grisli_g40.o: $(srcdir)/grisli_g40.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/MAR.o: $(srcdir)/MAR.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/nasaBasins.o: $(srcdir)/nasaBasins.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/RACMO2.o: $(srcdir)/RACMO2.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/Rignot13_BasalMelt.o: $(srcdir)/Rignot13_BasalMelt.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/sediments.o: $(srcdir)/sediments.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/stratigraphy.o: $(srcdir)/stratigraphy.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/topo_reconstructions.o: $(srcdir)/topo_reconstructions.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/topographies_grl.o: $(srcdir)/topographies_grl.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

obj_datasets =     $(objdir)/control.o \
				   $(objdir)/nml.o \
				   $(objdir)/gridding_datasets.o \
				   $(objdir)/generic.o \
				   $(objdir)/regions.o \
			       $(objdir)/AN1CRUST.o \
			       $(objdir)/bedmap2.o \
			       $(objdir)/CERES.o \
			       $(objdir)/climber2.o \
			       $(objdir)/climber3a.o \
			       $(objdir)/davini2015.o \
			       $(objdir)/ECMWF.o \
			       $(objdir)/ETOPO.o \
			       $(objdir)/GeothermalHeatFlux.o \
			       $(objdir)/GreenlandVelocity.o \
			       $(objdir)/grisli_g40.o \
			       $(objdir)/MAR.o \
			       $(objdir)/nasaBasins.o \
				   $(objdir)/RACMO2.o \
				   $(objdir)/Rignot13_BasalMelt.o \
			       $(objdir)/sediments.o \
			       $(objdir)/stratigraphy.o \
			       $(objdir)/topo_reconstructions.o \
			       $(objdir)/topographies_grl.o

obj_datasets_climber = $(objdir)/control.o \
				   $(objdir)/nml.o \
				   $(objdir)/gridding_datasets.o \
			       $(objdir)/climber2.o \
			       $(objdir)/climber3a.o

## Complete programs

gridder: $(obj_datasets) 
	$(FC) $(DFLAGS) $(FLAGS) -o gridder.x $^ gridder.f90 $(LFLAGS)
	@echo " "
	@echo "    gridder.x is ready."
	@echo " "

gridder_help: $(obj_datasets) 
	$(FC) $(DFLAGS) $(FLAGS) -o gridder_help.x $^ gridder_help.f90 $(LFLAGS)
	@echo " "
	@echo "    gridder_help.x is ready."
	@echo " "

gridder_climber: $(obj_datasets_climber) 
	$(FC) $(DFLAGS) $(FLAGS) -o gridder_climber.x $^ gridder_climber.f90 $(LFLAGS)
	@echo " "
	@echo "    gridder_climber.x is ready."
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

bedmap: $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -o bedmap2_netcdf.x $^ bedmap2_netcdf.f90 $(LFLAGS)
	@echo " "
	@echo "    bedmap2_netcdf.x is ready."
	@echo " "

points_to_latlon:
	$(FC) $(DFLAGS) $(FLAGS) -o points_to_latlon.x $^ points_to_latlon.f90 $(LFLAGS)
	@echo " "
	@echo "    points_to_latlon.x is ready."
	@echo " "

clean:
	rm -f *.x $(objdir)/*.o $(objdir)/*.mod

# cleanall: cleansico cleanrembo cleansicoX
