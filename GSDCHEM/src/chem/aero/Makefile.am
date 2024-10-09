noinst_LIBRARIES = libaero.a

libaero_a_SOURCES = aero_soa_vbs_mod.F90 aero_soa_vbs_data_mod.F90

libaero_a_FCFLAGS = -I . -I ../shr/ -I ../dust/

libaero_a_DEPENDENCIES = ../shr/libshr.a ../dust/libdust.a

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules


libaero_a-aero_soa_vbs_mod.$(OBJEXT)   : libaero_a-aero_soa_vbs_data_mod.$(OBJEXT)
