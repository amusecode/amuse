AC_CONFIG_MACRO_DIRS([shared/m4])
AC_INIT([amuse-lib-simple_hash], [1.0])

AMUSE_DETECT_OS()

AC_PROG_CC()

AC_CHECK_TARGET_TOOL([AR], [ar], [:])
AS_IF([test $AR = ":"], [AC_MSG_ERROR([ar command not found.])])

AC_PROG_RANLIB()
AS_IF([test $RANLIB = ":"], [AC_MSG_ERROR([ranlib command not found.])])

AC_CHECK_TOOL([PYTHON], [python3], [:])
AS_IF([test $PYTHON = ":"], [AC_MSG_WARN([Python 3 not found, cannot run tests.])])

AC_CONFIG_AUX_DIR([shared])
AC_PROG_INSTALL()

AC_CONFIG_FILES([config.mk])
AC_OUTPUT

