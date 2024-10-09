SUBDIRS = src
dist_doc_DATA = README
nodist_data_DATA = gsdchem.mk

ACLOCAL_AMFLAGS = -I m4

gsdchem.mk:
	@echo "# ESMF self-describing build dependency makefile fragment" > $(@)
	@echo >> $(@)
	@echo "ESMF_DEP_FRONT     = CHM"                           >> $(@)
	@echo "ESMF_DEP_INCPATH   = $(abs_top_builddir)/src/chem/" >> $(@)
	@echo "ESMF_DEP_CMPL_OBJS = "                              >> $(@)
	@echo "ESMF_DEP_LINK_OBJS = $(libdir)/libgsdchem.a"        >> $(@)
	@echo "ESMF_DEP_SHRD_PATH = "                              >> $(@)
	@echo "ESMF_DEP_SHRD_LIBS = mpp"                           >> $(@)
	@echo
	@echo "Finished generating ESMF self-describing build dependency makefile fragment:" $(@)
	@echo

uninstall-local:
	-rm -f $(datarootdir)/gsdchem.mk

CLEANFILES = gsdchem.mk
