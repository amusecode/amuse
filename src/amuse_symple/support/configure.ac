AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-symple], [1.0])

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Set the worker language
AC_PROG_CC()
AC_PROG_CXX()

# Find tools for creating static libraries
AC_CHECK_TARGET_TOOL([AR], [ar], [:])
AS_IF([test $AR = ":"], [AC_MSG_ERROR([ar command not found.])])

AC_PROG_RANLIB()
AS_IF([test $RANLIB = ":"], [AC_MSG_ERROR([ranlib command not found.])])


# Find AMUSE libraries
AMUSE_LIB_STOPCONDMPI()
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

