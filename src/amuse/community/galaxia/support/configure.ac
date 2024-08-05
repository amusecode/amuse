AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-galaxia], [1.0])

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Find the compiler
AC_PROG_FC()
AX_FC_ISO_C_BINDING()

# Find AMUSE libraries
AMUSE_LIB_FORSOCKETS()

AC_LANG_PUSH([Fortran])

# Find MPI
AX_MPI()

AX_LAPACK()     # calls AX_BLAS automatically

AC_LANG_POP([Fortran])

# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT
