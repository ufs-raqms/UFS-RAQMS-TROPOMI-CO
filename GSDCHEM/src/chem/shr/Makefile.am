noinst_LIBRARIES = libshr.a

libshr_a_SOURCES = chem_clock_mod.F90 chem_comm_mod.F90 chem_config_mod.F90 chem_const_mod.F90 \
	chem_data_mod.F90 chem_domain_mod.F90 chem_io_mod.F90 chem_iodata_mod.F90 chem_iolayout_mod.F90 \
	chem_model_mod.F90 chem_rc_mod.F90 chem_species_mod.F90 \
	chem_state_mod.F90 chem_tracers_mod.F90 chem_types_mod.F90 chem_raqms_mod.F90

libshr_a_FCFLAGS = -I .

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules


libshr_a-chem_clock_mod.$(OBJEXT): libshr_a-chem_rc_mod.$(OBJEXT) libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_comm_mod.$(OBJEXT): libshr_a-chem_rc_mod.$(OBJEXT) libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_config_mod.$(OBJEXT): libshr_a-chem_comm_mod.$(OBJEXT) libshr_a-chem_rc_mod.$(OBJEXT) \
	libshr_a-chem_species_mod.$(OBJEXT) libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_const_mod.$(OBJEXT): libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_data_mod.$(OBJEXT): libshr_a-chem_rc_mod.$(OBJEXT) libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_domain_mod.$(OBJEXT): libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_io_mod.$(OBJEXT): libshr_a-chem_comm_mod.$(OBJEXT) libshr_a-chem_model_mod.$(OBJEXT) \
	libshr_a-chem_rc_mod.$(OBJEXT) libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_iodata_mod.$(OBJEXT): libshr_a-chem_comm_mod.$(OBJEXT) libshr_a-chem_config_mod.$(OBJEXT) \
	libshr_a-chem_io_mod.$(OBJEXT) libshr_a-chem_model_mod.$(OBJEXT) libshr_a-chem_rc_mod.$(OBJEXT) \
	libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_model_mod.$(OBJEXT): libshr_a-chem_clock_mod.$(OBJEXT) libshr_a-chem_comm_mod.$(OBJEXT) \
	libshr_a-chem_config_mod.$(OBJEXT) libshr_a-chem_data_mod.$(OBJEXT) libshr_a-chem_domain_mod.$(OBJEXT) \
	libshr_a-chem_iolayout_mod.$(OBJEXT) libshr_a-chem_rc_mod.$(OBJEXT) libshr_a-chem_species_mod.$(OBJEXT) \
	libshr_a-chem_state_mod.$(OBJEXT) libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_state_mod.$(OBJEXT): libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_tracers_mod.$(OBJEXT): libshr_a-chem_rc_mod.$(OBJEXT) libshr_a-chem_config_mod.$(OBJEXT)
