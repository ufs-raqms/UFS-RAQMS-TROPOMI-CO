noinst_LIBRARIES = libcedsair.a

libcedsair_a_SOURCES = cedsair_data_mod.F90 cedsair_mod.F90

libcedsair_a_FCFLAGS = -I ../shr/

libcedsair_a_DEPENDENCIES = ../shr/libshr.a

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules


libcedsair_a-cedsair_data_mod.$(OBJEXT)   : libcedsair_a-cedsair_data_mod.$(OBJEXT) 
libcedsair_a-cedsair_mod.$(OBJEXT) : libcedsair_a-cedsair_data_mod.$(OBJEXT) libcedsair_a-cedsair_mod.$(OBJEXT)
