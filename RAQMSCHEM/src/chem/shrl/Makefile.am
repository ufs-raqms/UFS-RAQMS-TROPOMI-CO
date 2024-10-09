noinst_LIBRARIES = libshrl.a

libshrl_a_SOURCES = chem_clock_mod.F90 chem_types_mod.F90 chem_rc_mod.F90 chem_comm_mod.F90 chem_species_mod.F90 \
         chem_config_mod.F90 chem_const_mod.F90 chem_iolayout_mod.F90 chem_io_mod.F90 chem_iodata_mod.F90 \
         chem_model_mod.F90 chem_data_mod.F90  chem_domain_mod.F90 chem_state_mod.F90 

libshrl_a_FCFLAGS = -I .

clean-local:
	-rm -f *.mod  *.o *.a 

libshrl_a-chem_clock_mod.$(OBJEXT): libshrl_a-chem_rc_mod.$(OBJEXT) libshrl_a-chem_types_mod.$(OBJEXT)
libshrl_a_chem_config_mod.$(OBJEXT): libshrl_a-chem_comm_mod.$(OBJEXT) libshrl_a-chem_rc_mod.$(OBJEXT) \
        libshrl_a-chem_species_mod.$(OBJEXT) libshrl_a-chem_types_mod.$(OBJEXT)
libshrl_a_chem_comm_mod.$(OBJEXT): libshrl_a-chem_types_mod.$(OBJEXT) libshrl_a-chem_rc_mod.$(OBJEXT)
libshrl_a-chem_const_mod.$(OBJEXT): libshrl_a-chem_types_mod.$(OBJEXT)
libshrl_a-chem_io_mod.$(OBJEXT): libshrl_a-chem_comm_mod.$(OBJEXT) libshrl_a-chem_model_mod.$(OBJEXT) \
        libshrl_a-chem_rc_mod.$(OBJEXT) libshrl_a-chem_types_mod.$(OBJEXT)
libshrl_a-chem_iodata_mod.$(OBJEXT): libshrl_a-chem_comm_mod.$(OBJEXT) libshrl_a-chem_config_mod.$(OBJEXT) \
        libshrl_a-chem_io_mod.$(OBJEXT) libshrl_a-chem_model_mod.$(OBJEXT) libshrl_a-chem_rc_mod.$(OBJEXT) \
        libshrl_a-chem_types_mod.$(OBJEXT)
libshrl_a-chem_model_mod.$(OBJEXT): libshrl_a-chem_clock_mod.$(OBJEXT) libshrl_a-chem_comm_mod.$(OBJEXT) \
        libshrl_a-chem_config_mod.$(OBJEXT) libshrl_a-chem_data_mod.$(OBJEXT) libshrl_a-chem_domain_mod.$(OBJEXT) \
        libshrl_a-chem_iolayout_mod.$(OBJEXT) libshrl_a-chem_rc_mod.$(OBJEXT) libshrl_a-chem_species_mod.$(OBJEXT) \
        libshrl_a-chem_state_mod.$(OBJEXT) libshrl_a-chem_types_mod.$(OBJEXT)
libshrl_a-chem_data_mod.$(OBJEXT): libshrl_a-chem_rc_mod.$(OBJEXT) libshrl_a-chem_types_mod.$(OBJEXT)
libshrl_a-chem_state_mod.$(OBJEXT): libshrl_a-chem_types_mod.$(OBJEXT)

