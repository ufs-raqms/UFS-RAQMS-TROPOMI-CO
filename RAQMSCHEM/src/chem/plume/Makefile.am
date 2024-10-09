noinst_LIBRARIES = libplume.a

#OPT = -xCore-AVX512 -qopt-zmm-usage=high
 OPT = -O2 -fp-model source -ftz -align array64byte -qno-opt-dynamic-align

F77 = mpiifort $(OPT) -real-size 64 -fPIC -assume realloc_lhs -thread -m64 -mcmodel=small -extend-source  -assume byterecl -convert big_endian -safe-cray-ptr -nowarn -sox -fno-alias -auto -align array64byte -axavx -fp-stack-check -I../shr/ -I../../../../GSDCHEM/src/chem/shr/

libplume_a_SOURCES = raq_plume_data_mod.F90 raq_plume_zero_mod.F90 raq_plume_scalar_mod.F90 raqms_plume_mod.F90

libplume_a_FCFLAGS = -nowarn -I . -I../../../../GSDCHEM/src/chem/shr/ -I../shr/ 
libplume_a_FFFLAGS = $(OPT) -I . -I../../../../GSDCHEM/src/chem/shr/ -I../shr/ 

clean-local:
	-rm -f *.mod  *.o *.a 

libplume_a-raq_plume_data_mod.$(OBJEXT):
libplume_a-raq_plume_zero_mod.$(OBJEXT):
libplume_a-raq_plume_scalar_mod.$(OBJEXT): libplume_a-raq_plume_zero_mod.$(OBJEXT)
libplume_a-raqms_plume_mod.$(OBJEXT): libplume_a-raq_plume_scalar_mod.$(OBJEXT) libplume_a-raq_plume_data_mod.$(OBJEXT)
