noinst_LIBRARIES = libplume.a

libplume_a_SOURCES = plume_data_mod.F90 plume_rise_mod.F90 plume_zero_mod.F90 plume_scalar_mod.F90

libplume_a_FCFLAGS = -I ../shr/

libplume_a_DEPENDENCIES = ../shr/libshr.a

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules


libplume_a-plume_rise_mod.$(OBJEXT)   : libplume_a-plume_data_mod.$(OBJEXT) \
	libplume_a-plume_zero_mod.$(OBJEXT) libplume_a-plume_scalar_mod.$(OBJEXT)
libplume_a-plume_scalar_mod.$(OBJEXT) : libplume_a-plume_data_mod.$(OBJEXT) libplume_a-plume_zero_mod.$(OBJEXT)
