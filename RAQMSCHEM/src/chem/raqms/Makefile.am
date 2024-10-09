noinst_LIBRARIES = libraqms.a
#OPT = -xCore-AVX512 -qopt-zmm-usage=high
#OPT = -nowarn -C -g -ftrapuv -traceback -xCore-AVX512 -qopt-zmm-usage=high
#OPT = -nowarn -C -g -ftrapuv -traceback 
#OPT = -O
OPT = -O2 -fp-model source -traceback -ftz -align array64byte -qno-opt-dynamic-align
#OPT = -C -g -ftrapuv -traceback -O2 -fp-model source -ftz -align array64byte -qno-opt-dynamic-align

#F77 = mpiifort -real-size 64 -mcmodel=medium -extend-source -O2 -C -g -traceback -ftrapuv -assume byterecl -convert big_endian -I../shr/ -I../../../../GSDCHEM/src/chem/shr/
#F77 = mpiifort -real-size 64 -extend-source -O0 -C -g -traceback -ftrapuv -assume byterecl -convert big_endian -I../shr/ -I../../../../GSDCHEM/src/chem/shr/
#F77 = mpiifort -real-size 64 -fp-model precise -fPIC -assume realloc_lhs -thread -m64 -mcmodel=small -extend-source -O0 -C -g -traceback -ftrapuv -assume byterecl -convert big_endian -safe-cray-ptr -nowarn -sox -fno-alias -auto -align array64byte -axavx -fp-stack-check -I../shr/ -I../../../../GSDCHEM/src/chem/shr/
F77 = mpiifort $(OPT) -real-size 64 -fPIC -assume realloc_lhs -thread -m64 -mcmodel=small -extend-source  -assume byterecl -convert big_endian -safe-cray-ptr -nowarn -sox -fno-alias -auto -align array64byte -axavx  -I../shr/ -I../../../../GSDCHEM/src/chem/shr/
#F77 = mpiifort -real-size 64 -mcmodel=medium -extend-source -O2 -assume byterecl -convert big_endian -I../shr/ -I../../../../GSDCHEM/src/chem/shr/
#F77 = mpiifort -real-size 64 -extend-source -O2 -assume byterecl -convert big_endian -I../shr/ -I../../../../GSDCHEM/src/chem/shr/
#F77 = mpiifort -real-size 64 -extend-source -O0 -C -g -traceback -ftrapuv  -assume byterecl -convert big_endian -I../shr/ -I../../../../GSDCHEM/src/chem/shr/

libraqms_a_SOURCES  = vertmxv_mod.F90 computeaod.F90 calcaod.F90 calcaodgsi.F90 het_kenc.F exchange_o3vmr_mod.F seasalt_mod.F90 \
	raqms_mod.F90 chem_driver_3d_nmhc.F chemgsimod.F90 raqms_model_mod.F90 photcode_larc_47_nmhc_v0.7.F newcol.F \
	vertmx_driver_mod.F90 \
	combined_chemsource_nmhc.F wetdep_mod.F90 wetdep_fam_mod.F90 wetdep_prep_mod.F90 wetdep_alpha_mod.F90  \
	raq_dep_ctrans_grell_mod.F90 raq_dep_wet_ls_mod.F90

libraqms_a_FCFLAGS = -nowarn -I . -I../../../../GSDCHEM/src/chem/shr/ -I../shr/ -I../plume/ -I../../../../FV3/gfsphysics -I../../../../FV3/gfsphysics/physics -I../../../../FMS/fv3gfs
libraqms_a_FFFLAGS = $(OPT) -I . -I../../../../GSDCHEM/src/chem/shr/ -I../shr/ -I../plume/

libraqms_a_DEPENDENCIES = ../../../../GSDCHEM/src/chem/shr/libshr.a ../shr/libshr.a ../plume/libplume.a

clean-local:
	-rm -f *.mod *.o *.a

libraqms_a-computeaod.$(OBJEXT):
libraqms_a-calcaod.$(OBJEXT): libraqms_a-computeaod.$(OBJEXT)
libraqms_a-calcaodgsi.$(OBJEXT): 
libraqms_a-het_kenc.$(OBJEXT):
libraqms_a-vertmxv_mod.$(OBJEXT):
libraqms_a-vermx_driver_mod.$(OBJEXT): libraqms_a-vertmxv_mod.$(OBJEXT)
libraqms_a-wetdep_fam_mod.$(OBJEXT):
libraqms_a-wetdep_mod.$(OBJEXT): libraqms_a-wetdep_prep_mod.$(OBJEXT) libraqms_a-raq_dep_ctrans_grell_mod.$(OBJEXT) libraqms_a-raq_dep_wet_ls_mod.$(OBJEXT) libraqms_a-wetdep_fam_mod.$(OBJEXT) libraqms_a-vertmx_driver_mod.$(OBJEXT)
libraqms_a-wetprep_mod.$(OBJEXT):
libraqms_a-wetdep_alpha_mod.$(OBJEXT):
libraqms_a-raq_dep_ctrans_grell_mod.$(OBJEXT): libraqms_a-wetdep_alpha_mod.$(OBJEXT) libraqms_a-wetdep_fam_mod.$(OBJEXT)
libraqms_a-raq_dep_wet_ls_mod.$(OBJEXT): libraqms_a-wetdep_alpha_mod.$(OBJEXT) libraqms_a-wetdep_fam_mod.$(OBJEXT)
libraqms_a-photcode_larc_47_nmhc_v0.7.$(OBJEXT):
libraqms_a-exchange_o3vmr_mod.$(OBJEXT):
libraqms_a-combined_chemsource_nmhc.$(OBJEXT): libraqms_a-het_kenc.$(OBJEXT)
libraqms_a-chemgsimod.$(OBJEXT) : libraqms_a-chemgsimod.$(OBJEXT) libraqms_a-raqms_mod.$(OBJEXT)
libraqms_a-newcol.$(OBJEXT):
libraqms_a-raqms_model_mod.$(OBJEXT) : libraqms_a-raqms_mod.$(OBJEXT) libraqms_a-chemgsimod.$(OBJEXT)
libraqms_a-chem_driver_3d_nmhc.$(OBJEXT): libraqms_a-exchange_o3vmr_mod.$(OBJEXT) libraqms_a-seasalt_mod.$(OBJEXT)
libraqms_a-raqms_mod.$(OBJEXT) : libraqms_a-newcol.$(OBJEXT) libraqms_a-combined_chemsource_nmhc.$(OBJEXT) libraqms_a-chem_driver_3d_nmhc.$(OBJEXT) libraqms_a-wetdep_mod.$(OBJEXT)
