noinst_LIBRARIES = libdep.a

libdep_a_SOURCES = dep_mod.F90 dep_wet_ls_mod.F90 dep_dry_gocart_mod.F90 dep_dry_mod.F90 dep_simple_mod.F90 \
	dep_vertmx_mod.F90 dep_ctrans_grell_mod.F90 dep_cu_g3_mod.F90

libdep_a_FCFLAGS = -I . -I ../shr/ -I ../aero/ -I ../dust

libdep_a_DEPENDENCIES = ../shr/libshr.a ../aero/libaero.a 

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules


libdep_a-dep_mod.$(OBJEXT)     : libdep_a-dep_wet_ls_mod.$(OBJEXT) libdep_a-dep_dry_mod.$(OBJEXT) libdep_a-dep_ctrans_grell_mod.$(OBJEXT)
libdep_a-dep_dry_mod.$(OBJEXT) : libdep_a-dep_dry_gocart_mod.$(OBJEXT) libdep_a-dep_simple_mod.$(OBJEXT) libdep_a-dep_vertmx_mod.$(OBJEXT)
