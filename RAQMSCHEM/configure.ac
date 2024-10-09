#                                               -*- Autoconf -*- 
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([raqmschem], [0.1] , [allenl@ssec.wisc.edu])
AM_INIT_AUTOMAKE
#AM_EXTRA_RECURSIVE_TARGETS([fortran nuopc])

# Checks for programs.
AC_PROG_GREP
AC_PROG_SED
AC_PATH_PROG(CUT,cut)

AC_ARG_WITH([esmfmk],
  [AS_HELP_STRING([--with-esmfmk=<path-to-esmf.mk>],
    [set full path to makefile fragment from ESMF library.])],
  [], 
  [with_esmfmk=no])

AS_IF([test "x$with_esmfmk" != xno],
  AC_CHECK_FILE(["$with_esmfmk"],
    [AC_SUBST(ESMFMKFILE,["$with_esmfmk"])
     with_esmfmk=yes],
    [AC_MSG_WARN([File $with_esmfmk not found. Carry on.])
     with_esmfmk=no]))

AS_IF([test "x$with_esmfmk" = xno],
  [AC_CHECKING([environment for ESMFMKFILE])
   AC_CHECK_FILE(["${ESMFMKFILE}"],
      [AC_SUBST(ESMFMKFILE,["${ESMFMKFILE}"])
       AC_MSG_RESULT([getting ESMF Makefile fragment from environment variable ESMFMKFILE])
       with_esmfmk=yes],
      [AC_MSG_WARN([Makefile fragment from ESMF library not found.])])])


AC_PROG_CC


AS_IF([test "x$with_esmfmk" = xyes],
  [FCOMPILER=`$GREP ESMF_F90COMPILER ${ESMFMKFILE} | $SED -e 's/ESMF_F90COMPILER *= *//g'`
   for v in ESMF_F90COMPILEOPTS ESMF_F90COMPILEPATHS ESMF_F90COMPILECPPFLAGS
   do  
     FCFLAGS="${FCFLAGS} `$GREP $v ${ESMFMKFILE} 2>/dev/null | $CUT -d= -f2-`"
   done
   FCLIBS="${LDFLAGS}"
   for v in ESMF_F90LINKOPTS ESMF_F90LINKPATHS ESMF_F90LINKRPATHS ESMF_F90ESMFLINKLIBS
   do  
     FCLIBS="${FCLIBS} `$GREP $v ${ESMFMKFILE} 2>/dev/null | $CUT -d= -f2-`"
   done
   AC_PROG_FC([$FCOMPILER])
   AC_SUBST([FCLIBS])],
   [AC_PROG_FC([mpiifort])
    ="${FCFLAGS} -O2"])

AC_PROG_F77([mpiifort])
FFLAGS="-O2"

AS_IF([test "x$FC" = x], 
  [AC_MSG_ERROR([Cannot find a suitable Fortran compiler.])])

# default
AC_ARG_ENABLE([safe-fp],
  [AS_HELP_STRING([--enable-safe-fp],[Enable value-safe floating-point optimization for Intel compilers (default: yes).])],
  [],[])

AS_IF([test "x$enable_safe_fp" != xno],
 [AC_CHECKING([compiler for -fp-model option])
  save_FCFLAGS=${FCFLAGS}
  FCFLAGS="-fp-model source ${FCFLAGS}"
  FCFLAGS="-traceback -ftz -align array64byte -qno-opt-dynamic-align -fno-alias -nowarn -safe-cray-ptr ${FCFLAGS}" 
  AC_LANG_PUSH(Fortran)
  AC_COMPILE_IFELSE(AC_LANG_PROGRAM([]),
    [AC_MSG_RESULT([adding -fp-model precise option for results consistency])],
    [FCFLAGS=${save_FCFLAGS}
     AC_MSG_WARN([-fp-model option not recognized. Please add corresponding option for value-safe floating-point optimization.])])
  AC_LANG_POP()],[])

AC_PROG_INSTALL
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_PID_T

# Checks for library functions.
AC_CHECK_FUNCS([bzero clock_gettime strerror])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                src/Makefile
                src/chem/Makefile
                src/chem/shr/Makefile
                src/chem/plume/Makefile
                src/chem/raqms/Makefile])
AC_OUTPUT


AS_ECHO("")
AS_ECHO("${PACKAGE_NAME} version ${PACKAGE_VERSION} successfully configured.")
AS_ECHO("")
AS_ECHO("Configuration settings:")
AS_ECHO("=======================")
AS_ECHO("  CC        = ${CC}")
AS_ECHO("  CFLAGS    = ${CFLAGS}")
AS_ECHO("")
AS_ECHO("  CPP       = ${CPP}")
AS_ECHO("  CPPFLAGS  = ${CPPFLAGS}")
AS_ECHO("")
AS_ECHO("  FC        = ${FC}")
AS_ECHO("  FCFLAGS   = ${FCFLAGS}")
AS_ECHO("")
AS_ECHO("  LDFLAGS   = ${LDFLAGS}")
AS_ECHO("  LIBS      = ${LIBS}")
AS_ECHO("  FCLIBS    = ${FCLIBS}")
AS_ECHO("=======================")
AS_IF([test "x$with_esmfmk" != xno],
  [AS_ECHO("NOTE: Settings imported from ESMF Makefile fragment: ${ESMFMKFILE}")],[])
AS_ECHO("")
AS_ECHO("Type 'make' to build.")
