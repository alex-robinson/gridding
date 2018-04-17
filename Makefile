.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = .obj
bindir = ./
libdir = libs
testdir = tests

# Command-line options at make call
debug  ?= 0
openmp ?= 1

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

FC = gfortran

INC_NC  = -I/opt/local/include
LIB_NC  = -L/opt/local/lib -lnetcdff -lnetcdf

COORDROOT = /Users/robinson/models/EURICE/coordinates/libcoordinates
INC_COORD = -I${COORDROOT}/include
LIB_COORD = -L${COORDROOT}/include -lcoordinates

MKLROOT = /opt/intel/mkl
INC_MKL = -I${MKLROOT}/include
# For sequential running (no tbb library)
LIB_MKL =  -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

LISROOT = /Users/robinson/apps/lis/lis
INC_LIS = -I${LISROOT}/include 
LIB_LIS = -L${LISROOT}/lib/ -llis

FFLAGS_DEFAULT = -ffree-line-length-none -fbackslash -I$(objdir) -J$(objdir) $(INC_COORD)
FFLAGS_OPENMP  = $(FFLAGS_DEFAULT) -fopenmp

LFLAGS  = $(LIB_NC) $(LIB_COORD) 
DFLAGS_NODEBUG = -O2
DFLAGS_DEBUG   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
DFLAGS_PROFILE = -O2 -pg


# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

# Determine whether to use openmp flags 
FFLAGS = $(FFLAGS_DEFAULT)
ifeq ($(openmp), 1)
	FFLAGS = $(FFLAGS_OPENMP)
endif

###############################################
##							
## List of rules and source files
##
###############################################

include config/Makefile_gridder.mk

###############################################
##
## Compilation of complete programs
##
###############################################

gridder: $(gridder_libs) $(obj_datasets) 
	$(FC) $(DFLAGS) $(FFLAGS) -o gridder.x $^ gridder.f90 $(LFLAGS)
	@echo " "
	@echo "    gridder.x is ready."
	@echo " "

gridder_help: $(gridder_libs) $(obj_datasets) 
	$(FC) $(DFLAGS) $(FFLAGS) -o gridder_help.x $^ gridder_help.f90 $(LFLAGS)
	@echo " "
	@echo "    gridder_help.x is ready."
	@echo " "

gridder_climber: $(gridder_libs) $(obj_datasets_climber) 
	$(FC) $(DFLAGS) $(FFLAGS) -o gridder_climber.x $^ gridder_climber.f90 $(LFLAGS)
	@echo " "
	@echo "    gridder_climber.x is ready."
	@echo " "

bedmap: 
	$(FC) $(DFLAGS) $(FFLAGS) -o bedmap2_netcdf.x $^ bedmap2_netcdf.f90 $(LFLAGS)
	@echo " "
	@echo "    bedmap2_netcdf.x is ready."
	@echo " "

marclean: $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -o mar_rawclean.x $^ mar_rawclean.f90 $(LFLAGS)
	@echo " "
	@echo "    mar_rawclean.x is ready."
	@echo " "

marmonthly: $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -o mar_calcmonthly.x $^ mar_calcmonthly.f90 $(LFLAGS)
	@echo " "
	@echo "    mar_calcmonthly.x is ready."
	@echo " "

points_to_latlon:
	$(FC) $(DFLAGS) $(FFLAGS) -o points_to_latlon.x $^ points_to_latlon.f90 $(LFLAGS)
	@echo " "
	@echo "    points_to_latlon.x is ready."
	@echo " "

cores_to_grid:
	$(FC) $(DFLAGS) $(FFLAGS) -o cores_to_grid.x $^ cores_to_grid.f90 $(LFLAGS)
	@echo " "
	@echo "    cores_to_grid.x is ready."
	@echo " "

# NOTE: eventually ncio should be extracted from all external libraries (like coord?), and then
# all external libraries should be linked by adding them to the end of the program compilation 
# command, instead of just `$(objdir)/nml.o` as it is now. 

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make gridder : compiles gridder.x"
	@echo " make clean   : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
