AC_CONFIG_MACRO_DIRS([shared/m4])
AC_INIT([amuse-lib-amuse_mpi], [1.0])

AMUSE_DETECT_OS()

AC_PROG_CC()

AC_CHECK_TARGET_TOOL([AR], [ar], [:])
AS_IF([test $AR = ":"], [AC_MSG_ERROR([ar command not found.])])

AC_PROG_RANLIB()
AS_IF([test $RANLIB = ":"], [AC_MSG_ERROR([ranlib command not found.])])

AC_CONFIG_AUX_DIR([shared])
AC_PROG_INSTALL()

AC_LANG_PUSH(C)
AX_MPI()
AC_LANG_POP(C)

AC_CONFIG_FILES([config.mk])
AC_OUTPUT

