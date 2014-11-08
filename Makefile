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

# PATH options
objdir = .obj
srcdir = src
libdir = ..

# netcdf_inc = /usr/include
# netcdf_lib = /usr/lib
netcdf_inc = /opt/local/include
netcdf_lib = /opt/local/lib
netcdf_inc_ifort = /home/robinson/apps/netcdf/netcdf/include
netcdf_lib_ifort = /home/robinson/apps/netcdf/netcdf/lib

# Command-line options at make call
ifort ?= 0
debug ?= 0 

ifeq ($(ifort),1)
    FC = ifort 
else
    FC = gfortran
endif 

ifeq ($(ifort),1)
	## IFORT OPTIONS ##
	FLAGS        = -module $(objdir) -L$(objdir) -I$(netcdf_inc_ifort)
	LFLAGS		 = -L$(netcdf_lib_ifort) -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0 -heap-arrays
	    # -w 
	else
	    DFLAGS   = -vec-report0 -O3 -heap-arrays
	endif
else
	## GFORTRAN OPTIONS ##
	FLAGS        = -I$(objdir) -J$(objdir) -I$(netcdf_inc)
	LFLAGS		 = -L$(netcdf_lib) -lnetcdff -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow \
	               -fbacktrace -fcheck=all -fbackslash
	else
	    DFLAGS   = -O3 -fbackslash
	endif
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

test: $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -o test.x $^ test.f90 $(LFLAGS)
	@echo " "
	@echo "    test.x is ready."
	@echo " "

clean:
	rm -f gentopo_grl.x gentopo_ant.x mar_rawclean.x $(objdir)/*.o $(objdir)/*.mod

# cleanall: cleansico cleanrembo cleansicoX
