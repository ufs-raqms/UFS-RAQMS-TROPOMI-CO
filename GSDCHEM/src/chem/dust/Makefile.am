noinst_LIBRARIES = libdust.a

libdust_a_SOURCES = dust_afwa_mod.F90 dust_data_mod.F90 dust_fengsha_mod.F90 dust_gocart_mod.F90 dust_mod.F90

libdust_a_FCFLAGS = -I ../shr/

libdust_a_DEPENDENCIES = ../shr/libshr.a

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules


libdust_a-dust_afwa_mod.$(OBJEXT): libdust_a-dust_data_mod.$(OBJEXT)
libdust_a-dust_fengsha_mod.$(OBJEXT): libdust_a-dust_data_mod.$(OBJEXT)
libdust_a-dust_gocart_mod.$(OBJEXT): libdust_a-dust_data_mod.$(OBJEXT)
libdust_a-dust_mod.$(OBJEXT): libdust_a-dust_data_mod.$(OBJEXT) \
	libdust_a-dust_afwa_mod.$(OBJEXT) libdust_a-dust_fengsha_mod.$(OBJEXT) libdust_a-dust_gocart_mod.$(OBJEXT)
