AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-mesa_r15140], [1.0])

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Set CPU_COUNT
AX_COUNT_CPUS()
AC_SUBST(CPU_COUNT)

# Find the compiler(s)
AC_PROG_CC()

AC_PROG_CXX()

AC_PROG_FC()
AX_FC_ISO_C_BINDING()

# Find tools for creating static libraries
AC_CHECK_TARGET_TOOL([AR], [ar], [:])
AS_IF([test $AR = ":"], [AC_MSG_ERROR([ar command not found.])])

AC_PROG_RANLIB()
AS_IF([test $RANLIB = ":"], [AC_MSG_ERROR([ranlib command not found.])])

# Find something to download things with
AMUSE_DOWNLOAD()

# Find AMUSE libraries
AMUSE_LIB_FORSOCKETS()

# Find external dependencies
AC_LANG_PUSH([C])

AX_MPI()

AC_OPENMP()
AC_SUBST(OPENMP_CFLAGS)

AX_LAPACK()

AC_LANG_POP([C])


AC_LANG_PUSH([C++])

AX_MPI()

AC_OPENMP()
AC_SUBST(OPENMP_CFLAGS)

dnf AX_LAPACK()

AC_LANG_POP([C++])


AC_LANG_PUSH([Fortran])

AX_MPI()

AC_OPENMP()
AC_SUBST(OPENMP_FCFLAGS)

AX_LAPACK()
AX_LIB_HDF5()

AC_LANG_POP([Fortran])


# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT

