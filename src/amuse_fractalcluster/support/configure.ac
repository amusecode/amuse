AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-fractalcluster], [1.0])

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Find tools for creating static libraries
AC_CHECK_TARGET_TOOL([AR], [ar], [:])
AS_IF([test $AR = ":"], [AC_MSG_ERROR([ar command not found.])])

AC_PROG_RANLIB()
AS_IF([test $RANLIB = ":"], [AC_MSG_ERROR([ranlib command not found.])])


# Find the compiler
AC_PROG_FC()
AX_FC_ISO_C_BINDING()

AC_LANG_PUSH([Fortran])

# Find MPI
AX_MPI()

# Find AMUSE libraries
AMUSE_LIB_FORSOCKETS()

AC_LANG_POP([Fortran])


# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT

