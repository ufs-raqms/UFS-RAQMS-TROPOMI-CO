#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([gsdchem], [0.9.0], [raffaele.montuoro@noaa.gov])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

# Checks for programs.
AC_PROG_GREP
AC_PROG_SED
AC_PATH_PROG(CUT,cut)

# Checks for required ESMF library
AC_ARG_WITH([esmfmk],
  [AS_HELP_STRING([--with-esmfmk=<path-to-esmf.mk>],
    [set full path to ESMF makefile fragment.])],
  [],[with_esmfmk=no])

AS_CASE(["x$with_esmfmk"],
  [x|xauto|xcheck|xyes|xno], [esmfmk=${ESMFMKFILE}],
  [esmfmk=$withval]
)

AC_MSG_CHECKING([for ESMF library])
AS_IF([test "x$esmfmk" != x], [dnl
  AC_CHECK_FILE([$esmfmk], [AC_SUBST(ESMFMKFILE,["$esmfmk"])], [esmfmk=])],
  [])

AS_IF([test "x$esmfmk" = x], [AC_MSG_ERROR([Unable to find required ESMF library])],
  [dnl Retrieve ESMF compiler settings
   ESMF_FC=`$GREP ESMF_F90COMPILER ${ESMFMKFILE} | $SED -e 's/ESMF_F90COMPILER *= *//g'`
   ESMF_FC_INC=`$GREP ESMF_F90COMPILEPATHS ${ESMFMKFILE} 2>/dev/null | $CUT -d= -f2-`
  ]
)

# Looks for Fortran compiler
# Assume we are cross compiling by default to prevent issues on WCOSS-Cray
save_cross_compiling=${cross_compiling}
cross_compiling=yes
# Silence resulting warning
save_ac_tool_warned=${ac_tool_warned}
ac_tool_warned=yes
AC_PROG_FC([$FC $ESMF_FC mpiifort mpifort mpif90 mpixlf90_r mpxlf90_r])
# Restore original cross compilation settings
cross_compiling=${save_cross_compiling}
ac_tool_warned=${save_ac_tool_warned}

# Checks if Fortran compiler supports MPI
AC_MSG_CHECKING([whether MPI is supported])
AC_LANG_PUSH([Fortran])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[dnl
       use mpi
       integer, parameter :: comm = MPI_COMM_WORLD
  ])],
  [AC_MSG_RESULT(yes)],
  [dnl No MPI found. Disable compiler.
   FC=
   AC_MSG_RESULT(no)
   AC_MSG_ERROR([Cannot find a suitable MPI Fortran compiler.])]
)

AC_FC_PP_SRCEXT([F90])
AX_COMPILER_VENDOR

# Add configure options
# - safe floating point arithmetic
AC_ARG_ENABLE([safe-fp],
  [AS_HELP_STRING([--enable-safe-fp],[Enable value-safe floating-point arithmetic (default: yes).])],
  [],[])

# - debug
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],[Enable debug mode (default: no).])],
  [enable_debug=${enableval}], [enable_debug=no])

# Checks for debug/optimization options
AS_IF([test "x$enable_debug" = xyes],[
  AC_MSG_CHECKING([compiler for debug options])
  save_FCFLAGS=${FCFLAGS}
  FCFLAGS=`echo $FCFLAGS | $SED -e 's/-O@<:@0-9s@:>@*//g'`
  FCFLAGS=`echo $FCFLAGS | $SED -e 's/-g@<:@0-9@:>@*//g'`
  debug_opts=
  AS_CASE(["$ax_cv_fc_compiler_vendor"],
    [gnu|clang], [debug_opts="-Wall -Wno-unused -Wno-unused-dummy-argument -fbacktrace"],
    [ibm],       [debug_opts="-qcheck -qfullpath"],
    [intel],     [debug_opts="-traceback -check pointer,stack,uninit"],
    [portland],  [debug_opts="-traceback"],
    []            dnl no debug options available for other platforms
  )
  AS_IF([test "x$ac_cv_prog_fc_g" = xyes], [debug_opts="-g ${debug_opts}"], [])
  FCFLAGS="-O0 ${debug_opts} ${FCFLAGS}"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [],
    [FCFLAGS=${save_FCFLAGS}
     debug_opts="unknown. Please add debug flags to FCFLAGS"]
  )
  AC_MSG_RESULT([${debug_opts}])
  ], [dnl default optimization settings, if none present
  AC_MSG_CHECKING([compiler for optimization options])
  fcflags_noopt=
  perf_opts=
  for fc_flag in $FCFLAGS
  do
    AS_CASE([$fc_flag],
      [-O@<:@0-9s@:>@*|-Ofast], [perf_opts=$fc_flag], dnl find optimization settings
      [fcflags_noopt="$fcflags_noopt $fc_flag"]
    )
  done
  AS_IF([test "x$perf_opts" = x], [dnl add default optimization settings
    AS_CASE(["$ax_cv_fc_compiler_vendor"],
      [gnu|clang], [dnl
                    perf_opts="-O2"
                    AS_IF([test "x$enable_safe_fp" != xno],
                           [perf_opts="${perf_opts} -fno-unsafe-math-optimizations"],[])],
      [ibm],       [perf_opts="-O2"],
      [intel],     [dnl
                    perf_opts="-O2"
                    AS_IF([test "x$enable_safe_fp" != xno],
                          [perf_opts="${perf_opts} -fp-model source -ftz -align array64byte -xCORE-AVX2 -qno-opt-dynamic-align"],[])],
      [portland],  [dnl
                    perf_opts="-O2"
                    AS_IF([test "x$enable_safe_fp" != xno],
                          [perf_opts="${perf_opts} -Kieee"],[])],
      []            dnl no debug options available for other platforms
    )], [])
  AS_IF([test "x$perf_opts" = x], [dnl optization flags unknown
    perf_opts=unknown],
    [dnl test optimization flags
    save_FCFLAGS=${fcflags_noopt}
    FCFLAGS="${fcflags_noopt} ${perf_opts}"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
      [],
      [FCFLAGS=${save_FCFLAGS}
       perf_opts=unknown]
    )
    ]
  )
  AC_MSG_RESULT([$perf_opts])
  ]
)

# Checks for compiler characteristics
AC_FC_LINE_LENGTH([unlimited])
AC_FC_MODULE_EXTENSION
AC_FC_MODULE_FLAG

# Adds ESMF include paths
fc_inc=
for incpath in ${ESMF_FC_INC}
do
  incpath=`echo "${incpath}" | $SED -e's/^-@<:@A-Z,a-z@:>@*\//\//g'`
  fc_inc="${fc_inc}${FC_MODINC}${incpath} "
done

AC_MSG_CHECKING([whether ESMF include flags work])
FCFLAGS="${fc_inc}${FCFLAGS}"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
  [AC_MSG_RESULT([yes])],
  [dnl ESMF-derived include flags may be corrupted. Abort.
   AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unrecognized include flags])]
)

# Checks whether ESMF is compatible with Fortran compiler
AC_MSG_CHECKING([whether ESMF works])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[dnl
       use esmf
  ])],
  [AC_MSG_RESULT(yes)],
  [dnl ESMF is incompatible with selected Fortran compiler or not found. Abort.
   AC_MSG_RESULT(no)
   AC_MSG_FAILURE([ESMF library incompatible or not found])]
)

AC_LANG_POP

# Checks for auxiliary build tools
AC_PROG_INSTALL
AC_PROG_RANLIB
AM_PROG_AR

# Output Makefiles
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/chem/Makefile
                 src/chem/shr/Makefile
                 src/chem/gocart/Makefile
                 src/chem/plume/Makefile
                 src/chem/vash/Makefile
                 src/chem/seas/Makefile
                 src/chem/dust/Makefile
                 src/chem/aero/Makefile
                 src/chem/dep/Makefile
                 src/chem/opt/Makefile])
AC_OUTPUT


AS_ECHO("")
AS_ECHO("${PACKAGE_NAME} version ${PACKAGE_VERSION} successfully configured.")
AS_ECHO("")
AS_ECHO("Compiler settings:")
AS_ECHO("=======================")
AS_ECHO("  FC      = ${FC}")
AS_ECHO("  FCFLAGS = ${FCFLAGS}")
AS_ECHO("=======================")
AS_ECHO("NOTE: Settings imported from ESMF Makefile fragment: ${ESMFMKFILE}")
AS_ECHO("")
AS_ECHO("Type 'make' to build.")
