AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-mosse], [1.0])

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Set the worker language
AC_LANG([Fortran])

# Set CPU_COUNT
AX_COUNT_CPUS()
AC_SUBST(CPU_COUNT)

# Find the compiler
AC_PROG_FC()
AX_FC_ISO_C_BINDING()

# Find AMUSE libraries
AMUSE_LIB_FORSOCKETS()

# Find external dependencies
AC_LANG_PUSH([Fortran])
AX_MPI()
AC_LANG_POP([Fortran])

# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT