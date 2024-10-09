SUBDIRS = shr plume vash dust aero dep seas opt gocart .

lib_LIBRARIES = libgsdchem.a

libgsdchem_a_LIBADD  = shr/*.$(OBJEXT) plume/*.$(OBJEXT) vash/*.$(OBJEXT) dust/*.$(OBJEXT) dep/*.$(OBJEXT) aero/*.$(OBJEXT) seas/*.$(OBJEXT) opt/*.$(OBJEXT) gocart/*.$(OBJEXT)

libgsdchem_a_FCFLAGS  = -I . -I shr -I gocart -I plume -I vash -I dust -I dep -I aero -I seas -I opt

libgsdchem_a_SOURCES = chemCap.F90 chem_methods.F90

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules


libgsdchem_a-chemCap.$(OBJEXT): libgsdchem_a-chem_methods.$(OBJEXT)
libgsdchem_a-chem_methods.$(OBJEXT): $(builddir)/shr/libshr.a
