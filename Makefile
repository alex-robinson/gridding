.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make GRL        : compiles the main program gentopo_GRL.x"
	@echo " make clean      : cleans object files"
	@echo ""

objdir = .obj

ifort ?= 0
debug ?= 0 

ifeq ($(ifort),1)
    FC = ifort 
else
    FC = gfortran
endif 

ifeq ($(ifort),1)
	## IFORT OPTIONS ##
	FLAGS        = -module $(objdir) -L$(objdir) -I/home/robinson/apps/netcdf/netcdf/include
	LFLAGS		 = -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -w -C -traceback -ftrapuv -fpe0 -check all -vec-report0
	else
	    DFLAGS   = -vec-report0 -O3
	endif
else
	## GFORTRAN OPTIONS ##
	FLAGS        = -I$(objdir) -J$(objdir) -I/opt/local/include
	LFLAGS		 = -L/opt/local/lib -lnetcdff -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
	else
	    DFLAGS   = -O3
	endif
endif

## Individual libraries or modules ##
$(objdir)/ncio.o: ../ncio/ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/planet.o: ../coord/planet.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/geodesic.o: ../coord/geodesic.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/projection_oblimap2.o: ../coord/projection_oblimap2.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/coordinates.o: ../coord/coordinates.f90 $(objdir)/projection_oblimap2.o \
	                     $(objdir)/geodesic.o $(objdir)/planet.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp1D.o: ../coord/interp1D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp2D.o: ../coord/interp2D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp_time.o: ../coord/interp_time.f90 $(objdir)/interp1D.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/gridding_datasets.o: gridding_datasets.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

## Complete programs

# Program to test interpolations of CCSM3 data
GRL: $(objdir)/ncio.o $(objdir)/geodesic.o $(objdir)/planet.o \
	         $(objdir)/projection_oblimap2.o $(objdir)/coordinates.o  \
	         $(objdir)/interp1D.o $(objdir)/interp_time.o $(objdir)/gridding_datasets.o
	$(FC) $(DFLAGS) $(FLAGS) -o gentopo_GRL.x $^ gentopo_GRL.f90 $(LFLAGS)
	@echo " "
	@echo "    gentopo_grl.x is ready."
	@echo " "

test: $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -o test.x $^ test.f90 $(LFLAGS)
	@echo " "
	@echo "    test.x is ready."
	@echo " "

clean:
	rm -f gentopo_grl.x $(objdir)/*.o $(objdir)/*.mod

# cleanall: cleansico cleanrembo cleansicoX
