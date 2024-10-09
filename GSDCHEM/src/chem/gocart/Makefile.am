noinst_LIBRARIES = libgocart.a

libgocart_a_SOURCES  = gocart_mod.F90
libgocart_a_SOURCES += gocart_aerosols_mod.F90
libgocart_a_SOURCES += gocart_chem_mod.F90
libgocart_a_SOURCES += gocart_diag_mod.F90
libgocart_a_SOURCES += gocart_dmsemis_mod.F90
libgocart_a_SOURCES += gocart_model_mod.F90
libgocart_a_SOURCES += gocart_prep_mod.F90
libgocart_a_SOURCES += gocart_settling_mod.F90

libgocart_a_FCFLAGS = -I . -I ../shr -I ../aero -I ../dep -I ../dust -I ../opt -I ../plume -I ../seas -I ../vash

libgocart_a_DEPENDENCIES = ../shr/libshr.a ../dep/libdep.a ../dust/libdust.a ../opt/libopt.a ../plume/libplume.a ../seas/libseas.a ../vash/libvash.a

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules


libgocart_a-gocart_mod.$(OBJEXT) : libgocart_a-gocart_aerosols_mod.$(OBJEXT) libgocart_a-gocart_diag_mod.$(OBJEXT) \
	libgocart_a-gocart_dmsemis_mod.$(OBJEXT) libgocart_a-gocart_chem_mod.$(OBJEXT) libgocart_a-gocart_prep_mod.$(OBJEXT) \
	libgocart_a-gocart_settling_mod.$(OBJEXT)
libgocart_a-gocart_model_mod.$(OBJEXT) : libgocart_a-gocart_mod.$(OBJEXT)
