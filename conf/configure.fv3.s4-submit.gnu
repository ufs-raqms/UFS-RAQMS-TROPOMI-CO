## NEMS configuration file
##
## Platform: s4-cardinal
## Compiler: Intel with intelmpi

SHELL=/bin/sh

################################################################################
## Include the common configuration parts

ifdef InNemsMakefile
include $(TOP)/conf/configure.nems.NUOPC
endif

###################### PHYS_MODE ##### CHEM_MODE ###############################
#
#
#

PHYS_MODE       =compile
CHEM_MODE       =compile
ifeq ($(PHYS_MODE),compile)
       PHYS_LIB = $(TOP)/atmos/gsm/gsmphys
       PHYS_INC = $(TOP)/atmos/gsm/gsmphys
       PHYS_DIR = $(TOP)/atmos/gsm/gsmphys
endif
ifeq ($(CHEM_MODE),compile)
        CHEM_LIB = $(TOP)/chem
        CHEM_INC = $(TOP)/chem/gocart/src/Config/
        CHEM_DIR = $(TOP)/chem
       CHEM_MOD = $(TOP)/chem/gocart/${ARCH}/include
      ESMADIR = chem/gocart
endif

############
# commands #
############
FC = mpiifort
CC = mpicc
CXX = mpicxx
LD = mpiifort -mkl=sequential

#########
# flags #
#########
# default is 64-bit OpenMP non-hydrostatic build using AVX2
DEBUG =
#DEBUG = Y
REPRO =
VERBOSE =
OPENMP = N
AVX2 = Y
# for amd nodes
#AVX2 = N
HYDRO = N

include       $(ESMFMKFILE)
ESMF_INC    = $(ESMF_F90COMPILEPATHS)
ifneq ($(NEMSIO_INC),)
NEMSIOINC = -I$(NEMSIO_INC)
endif

NCEPLIBS = $(NEMSIO_LIB) $(BACIO_LIB4) $(SP_LIBd) $(W3EMC_LIBd) $(W3NCO_LIBd)


##############################################
# Need to use at least GNU Make version 3.81 #
##############################################
need := 3.81
ok := $(filter $(need),$(firstword $(sort $(MAKE_VERSION) $(need))))
ifneq ($(need),$(ok))
$(error Need at least make version $(need).  Load module gmake/3.81)
endif
NETCDF_LIB :=
NETCDF_ROOT = $(SSEC_NETCDF4_DIR)
NETCDF4 = $(SSEC_NETCDF4_DIR)
INCLUDE = -I$(NETCDF_ROOT)/include
NETCDF_INC = -I$(NETCDF_ROOT)/include
ifneq ($(findstring netcdf/4,$(LOADEDMODULES)),)
  NETCDF_LIB += -L$(NETCDF4)/lib -lnetcdff -lnetcdf
else
  NETCDF_LIB += -L$(NETCDF4)/lib -lnetcdff -lnetcdf
endif

FPPFLAGS := -fpp -Wp,-w $(INCLUDE)
CFLAGS := $(INCLUDE)

FFLAGS := $(INCLUDE) -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte

ifeq ($(HYDRO),Y)
CPPDEFS += -Duse_libMPI -Duse_netCDF -DSPMD -DUSE_LOG_DIAG_FIELD_INFO -Duse_LARGEFILE -DUSE_GFSL63 -DGFS_PHYS
else
CPPDEFS += -Duse_libMPI -Duse_netCDF -DSPMD -DUSE_LOG_DIAG_FIELD_INFO -Duse_LARGEFILE -DUSE_GFSL63 -DGFS_PHYS -DMOIST_CAPPA -DUSE_COND
endif

CPPDEFS += -DNEW_TAUCTMAX -DINTERNAL_FILE_NML
CPPDEFS += -Duse_WRTCOMP

ifeq ($(32BIT),Y)
CPPDEFS += -DOVERLOAD_R4
FFLAGS += -i4 -real-size 32
else
FFLAGS += -i4 -real-size 64 -no-prec-div -no-prec-sqrt
endif
# for cardinal turn off ajl
FFLAGS += -axavx 
CFLAGS += -axavx 
#  FFLAGS += -xCore-AVX512 -qopt-zmm-usage=high
# FFLAGS += -xCore-AVX2
#FFLAGS += -mcmodel=medium
#CFLAGS += -mcmodel=medium

FFLAGS_OPT = -O2 -traceback -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3
FFLAGS_REPRO = -O2 -debug minimal -fp-model source -qoverride-limits -g -traceback
FFLAGS_DEBUG = -g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -traceback -ftrapuv
FFLAGS_DEBUGS = -traceback

TRANSCENDENTALS := -fast-transcendentals
#FFLAGS_OPENMP = -qopenmp
FFLAGS_VERBOSE = -v -V -what

CFLAGS += -D__IFC -sox -fp-model source

CFLAGS_OPT = -O2 -debug minimal
CFLAGS_REPRO = -O2 -debug minimal
#CFLAGS_OPENMP = -qopenmp
CFLAGS_DEBUG = -O0 -g -ftrapuv -traceback

# Optional Testing compile flags.  Mutually exclusive from DEBUG, REPRO, and OPT
# *_TEST will match the production if no new option(s) is(are) to be tested.
#FFLAGS_TEST = -O3 -debug minimal -fp-model source -qoverride-limits
FFLAGS_TEST = -O2 -debug minimal -fp-model source -qoverride-limits
CFLAGS_TEST = -O2

LDFLAGS :=
LDFLAGS_OPENMP := -qopenmp
LDFLAGS_VERBOSE := -Wl,-V,--verbose,-cref,-M

# start with blank LIBS
LIBS :=

ifneq ($(REPRO),)
CFLAGS += $(CFLAGS_REPRO)
FFLAGS += $(FFLAGS_REPRO)
FAST :=
else ifneq ($(DEBUG),)
CFLAGS += $(CFLAGS_DEBUG)
FFLAGS += $(FFLAGS_DEBUG)
FAST :=
else ifneq ($(TEST),)
CFLAGS += $(CFLAGS_TEST)
FFLAGS += $(FFLAGS_TEST)
FAST :=
else
CFLAGS += $(CFLAGS_OPT)
FFLAGS += $(FFLAGS_OPT)
FAST := $(TRANSCENDENTALS)
endif

ifneq ($(OPENMP),)
CFLAGS += $(CFLAGS_OPENMP)
FFLAGS += $(FFLAGS_OPENMP)
LDFLAGS += $(LDFLAGS_OPENMP)
endif

ifneq ($(VERBOSE),)
CFLAGS += $(CFLAGS_VERBOSE)
FFLAGS += $(FFLAGS_VERBOSE)
LDFLAGS += $(LDFLAGS_VERBOSE)
endif

LDFLAGS += $(LIBS)

ifdef InNemsMakefile
FFLAGS += $(ESMF_INC)
CPPFLAGS += -traditional
EXTLIBS = $(NCEPLIBS) $(ESMF_LIB) $(LDFLAGS) $(NETCDF_LIB)
endif
