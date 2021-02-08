## Individual libraries or modules ##

$(objdir)/ncio.o: $(libdir)/ncio.f90
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(libdir)/nml.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/control.o: $(srcdir)/control.f90 $(objdir)/nml.o
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/gridding_datasets.o: $(srcdir)/gridding_datasets.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/regions.o: $(srcdir)/regions.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/generic.o: $(srcdir)/generic.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/AN1CRUST.o: $(srcdir)/AN1CRUST.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/AntarcticaVelocity.o: $(srcdir)/AntarcticaVelocity.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/badgeley2020.o: $(srcdir)/badgeley2020.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/bedmap2.o: $(srcdir)/bedmap2.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/buizert2018.o: $(srcdir)/buizert2018.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/CERES.o: $(srcdir)/CERES.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/climber2.o: $(srcdir)/climber2.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/climber3a.o: $(srcdir)/climber3a.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/davini2015.o: $(srcdir)/davini2015.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ECMWF.o: $(srcdir)/ECMWF.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ETOPO.o: $(srcdir)/ETOPO.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/GeothermalHeatFlux.o: $(srcdir)/GeothermalHeatFlux.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/GreenlandVelocity.o: $(srcdir)/GreenlandVelocity.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/grisli_g40.o: $(srcdir)/grisli_g40.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/MAR.o: $(srcdir)/MAR.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/Morlighem2013.o: $(srcdir)/Morlighem2013.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/nasaBasins.o: $(srcdir)/nasaBasins.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/pmip3.o: $(srcdir)/pmip3.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/RACMO2.o: $(srcdir)/RACMO2.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/Rignot13_BasalMelt.o: $(srcdir)/Rignot13_BasalMelt.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rtopo.o: $(srcdir)/rtopo.f90 $(objdir)/generic.o
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/sediments.o: $(srcdir)/sediments.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/stratigraphy.o: $(srcdir)/stratigraphy.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/topo_reconstructions.o: $(srcdir)/topo_reconstructions.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/topographies_grl.o: $(srcdir)/topographies_grl.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/vavrus2018.o: $(srcdir)/vavrus2018.f90
		$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

obj_datasets =     $(objdir)/control.o \
									 $(objdir)/gridding_datasets.o \
									 $(objdir)/generic.o \
									 $(objdir)/regions.o \
									 $(objdir)/AN1CRUST.o \
                   					 $(objdir)/AntarcticaVelocity.o \
                   					 $(objdir)/badgeley2020.o \
									 $(objdir)/bedmap2.o \
									 $(objdir)/buizert2018.o \
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
									 $(objdir)/Morlighem2013.o \
									 $(objdir)/nasaBasins.o \
									 $(objdir)/RACMO2.o \
									 $(objdir)/Rignot13_BasalMelt.o \
									 $(objdir)/rtopo.o \
									 $(objdir)/sediments.o \
									 $(objdir)/stratigraphy.o \
									 $(objdir)/topo_reconstructions.o \
									 $(objdir)/topographies_grl.o \
									 $(objdir)/vavrus2018.o

obj_datasets_climber = $(objdir)/control.o \
									 $(objdir)/gridding_datasets.o \
									 $(objdir)/climber2.o \
									 $(objdir)/climber3a.o

obj_datasets_pmip3 = $(objdir)/control.o \
									 $(objdir)/gridding_datasets.o \
									 $(objdir)/pmip3.o

gridder_libs = $(objdir)/ncio.o $(objdir)/nml.o
