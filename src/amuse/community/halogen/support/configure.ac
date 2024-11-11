AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-halogen], [1.0])

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Set CPU_COUNT
AX_COUNT_CPUS()
AC_SUBST(CPU_COUNT)

# Set the worker language
AC_PROG_CC()
AC_PROG_CXX()

# Find AMUSE libraries
AMUSE_LIB_AMUSE_MPI()

# Find external dependencies
AC_LANG_PUSH([C])
AX_MPI()
AC_LANG_POP([C])

AC_LANG_PUSH([C++])
AX_MPI()
AC_LANG_POP([C++])


# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT

