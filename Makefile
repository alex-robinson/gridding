.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
objdir = .obj
libdir = ..
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
    INC_COORD = -I/Users/robinson/models/EURICE/coord/.obj
	LIB_COORD = /Users/robinson/models/EURICE/coord/libcoordinates.a

    FLAGS  = -I$(objdir) -J$(objdir) $(INC_COORD) $(INC_NC) 
    LFLAGS = $(LIB_COORD) $(LIB_NC)

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
    INC_COORD = -I/iplex/01/tumble/robinson/EURICE/coord/.obj
	LIB_COORD = /iplex/01/tumble/robinson/EURICE/coord/libcoordinates.a

    FLAGS    = -module $(objdir) -L$(objdir) $(INC_COORD) $(INC_NC) 
    LFLAGS   = $(LIB_COORD) $(LIB_NC)

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

$(objdir)/gridding_datasets.o: $(srcdir)/gridding_datasets.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/Bamber13.o: $(srcdir)/Bamber13.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/CERES.o: $(srcdir)/CERES.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/ECMWF.o: $(srcdir)/ECMWF.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/MAR.o: $(srcdir)/MAR.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/RACMO2.o: $(srcdir)/RACMO2.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/Rignot13_BasalMelt.o: $(srcdir)/Rignot13_BasalMelt.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/bedmap2.o: $(srcdir)/bedmap2.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/nasaBasins.o: $(srcdir)/nasaBasins.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

grisli_common = $(objdir)/runparam_mod.o $(objdir)/3D-physique-gen_mod.o

obj_datasets_GRL = $(objdir)/gridding_datasets.o \
			   	   $(objdir)/Bamber13.o \
			       $(objdir)/CERES.o \
			       $(objdir)/ECMWF.o \
			       $(objdir)/MAR.o \
			       $(objdir)/RACMO2.o \
			       $(objdir)/nasaBasins.o

obj_datasets_ANT = $(objdir)/gridding_datasets.o \
			       $(objdir)/CERES.o \
			       $(objdir)/ECMWF.o \
			       $(objdir)/MAR.o \
			       $(objdir)/RACMO2.o \
			       $(objdir)/Rignot13_BasalMelt.o \
			       $(objdir)/bedmap2.o \
			       $(objdir)/nasaBasins.o

## Complete programs

ANT: $(obj_datasets_ANT)
	$(FC) $(DFLAGS) $(FLAGS) -o gentopo_ANT.x $^ gentopo_ANT.f90 $(LFLAGS)
	@echo " "
	@echo "    gentopo_ANT.x is ready."
	@echo " "

GRL: $(obj_datasets_GRL) 
	$(FC) $(DFLAGS) $(FLAGS) -o gentopo_GRL.x $^ gentopo_GRL.f90 $(LFLAGS)
	@echo " "
	@echo "    gentopo_GRL.x is ready."
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
	$(FC) $(DFLAGS) $(FLAGS) -o bedmap2_netcdf.x $^ bedmap2_netcdf.f90 ../coord/libcoordinates.a $(LFLAGS)
	@echo " "
	@echo "    bedmap2_netcdf.x is ready."
	@echo " "

test: $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -o test.x $^ test.f90 $(LFLAGS)
	@echo " "
	@echo "    test.x is ready."
	@echo " "

clean:
	rm -f gentopo_grl.x gentopo_ant.x mar_rawclean.x $(objdir)/*.o $(objdir)/*.mod

# cleanall: cleansico cleanrembo cleansicoX
