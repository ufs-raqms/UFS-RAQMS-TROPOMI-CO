SUBDIRS = src
dist_doc_DATA = README
EXTRA_DIST = test

install-data-hook: $(datarootdir)/raqmschem.mk

$(datarootdir)/raqmschem.mk:
	@echo "# ESMF self-describing build dependency makefile fragment" > $@
	@echo >> $@
	@echo "ESMF_DEP_FRONT     = RAQMSCHM"                    >> $@
	@echo "ESMF_DEP_INCPATH   = $(abs_top_builddir)/src/chem/" >> $@
	@echo "ESMF_DEP_CMPL_OBJS = "                       >> $@
	@echo "ESMF_DEP_LINK_OBJS = $(libdir)/libraqmschem.a" >> $@
	@echo "ESMF_DEP_SHRD_PATH = "                   >> $@
	@echo "ESMF_DEP_SHRD_LIBS = mpp"                >> $@
	@echo
	@echo "Finished generating ESMF self-describing build dependency makefile fragment:" $@
	@echo
