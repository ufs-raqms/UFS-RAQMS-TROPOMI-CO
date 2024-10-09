# Location of the ESMF makefile fragment for this component:
raqmschem_mk = $(RAQMSCHEM_BINDIR)/raqmschem.mk
all_component_mk_files+=$(raqmschem_mk)

# Location of source code and installation
RAQMSCHEM_SRCDIR?=$(ROOTDIR)/RAQMSCHEM
RAQMSCHEM_BINDIR?=$(ROOTDIR)/RAQMSCHEM_INSTALL

# Make sure the expected directories exist and are non-empty:
$(call require_dir,$(RAQMSCHEM_SRCDIR),RAQMSCHEM source directory)

RAQMSCHEM_ALL_OPTS= \
  COMP_SRCDIR="$(RAQMSCHEM_SRCDIR)" \
  COMP_BINDIR="$(RAQMSCHEM_BINDIR)" \
  FMS_BINDIR="$(FMS_BINDIR)" \
  MACHINE_ID="$(MACHINE_ID)"

########################################################################

# Rule for building this component:

build_RAQMSCHEM: $(raqmschem_mk)

$(raqmschem_mk): configure
	$(MODULE_LOGIC) ; export $(RAQMSCHEM_ALL_OPTS)                  ; \
	set -e                                                        ; \
	cd $(RAQMSCHEM_SRCDIR)                                          ; \
	./configure --prefix=$(RAQMSCHEM_BINDIR)                          \
	  --datarootdir=$(RAQMSCHEM_BINDIR) --libdir=$(RAQMSCHEM_BINDIR)
	+$(MODULE_LOGIC) ; cd $(RAQMSCHEM_SRCDIR) ; exec $(MAKE)             \
	  $(RAQMSCHEM_ALL_OPTS)
	+$(MODULE_LOGIC) ; cd $(RAQMSCHEM_SRCDIR) ; exec $(MAKE)             \
	  $(RAQMSCHEM_ALL_OPTS) install
	test -d "$(RAQMSCHEM_BINDIR)"

########################################################################

# Rule for cleaning the SRCDIR and BINDIR:

clean_RAQMSCHEM:
	+cd $(RAQMSCHEM_SRCDIR) ; test -f Makefile && exec $(MAKE) -k distclean || echo "Nothing to clean up"

distclean_RAQMSCHEM: clean_RAQMSCHEM
	rm -rf $(RAQMSCHEM_BINDIR) $(raqmschem_mk)
