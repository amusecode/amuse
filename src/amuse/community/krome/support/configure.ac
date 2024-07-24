AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-krome], [1.0])

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Set CPU_COUNT
AX_COUNT_CPUS()
AC_SUBST(CPU_COUNT)

# Find the compiler(s)
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
AC_LANG_PUSH([Fortran])

AX_MPI()

# TODO: it looks like Krome can use Intel MKL, and conda-forge packages it

AC_LANG_POP([Fortran])

# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT

