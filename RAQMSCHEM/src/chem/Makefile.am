SUBDIRS = shr plume raqms .
#OPT = -axavx -xCore-AVX512 -qopt-zmm-usage=high
OPT = -O2 -fp-model source -traceback -ftz -align array64byte -qno-opt-dynamic-align

lib_LIBRARIES = libraqmschem.a

include $(ESMFMKFILE)

libraqmschem_a_LIBADD  = shr/*.$(OBJEXT) plume/*.$(OBJEXT) raqms/*.$(OBJEXT)


libraqmschem_a_FCFLAGS  = -nowarn -I . -I shr -I plume -I../../../FV3/gfsphysics -I../../../GSDCHEM/src/chem/shr -I raqms
libraqmschem_a_FFFLAGS  = $(OPT) -nowarn -I . -I shr -I plume -I../../../FV3/gfsphysics -I../../../GSDCHEM/src/chem/shr -I raqms
libraqmschem_a_FCFLAGS += $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS)

libraqmschem_a_SOURCES = raqmschemCap.F90 raqmschem_methods.F90

clean-local:
	-rm -f *.mod *.o *.a

libraqmschem_a-raqmschemCap.$(OBJEXT): libraqmschem_a-raqmschem_methods.$(OBJEXT)
libraqmschem_a-raqmschem_methods.$(OBJEXT): $(builddir)/shr/libshr.a $(builddir)/plume/libplume.a
