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

# The AX_LAPACK macro (which also calls AX_BLAS) is designed to be called by a C wrapper
# So, we need to do this is order to detect blas/lapack even if we have a pure fortran code
AC_LANG_PUSH([C])
AX_LAPACK()
AC_LANG_POP([C])

AC_LANG_PUSH([Fortran])

# Find MPI
AX_MPI()

AC_LANG_POP([Fortran])

# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT
