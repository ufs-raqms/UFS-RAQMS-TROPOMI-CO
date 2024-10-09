noinst_LIBRARIES = libvash.a

libvash_a_SOURCES = vash_settling_mod.F90

libvash_a_FCFLAGS = -I ../shr/

libvash_a_DEPENDENCIES = ../shr/libshr.a 

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules
