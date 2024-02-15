AC_CONFIG_MACRO_DIRS([shared/m4])
AC_INIT([amuse-lib-forsockets], [1.0])

AC_PROG_CC()
AC_PROG_FC()
AX_FC_ISO_C_BINDING()

AC_CHECK_TARGET_TOOL([AR], [ar], [:])
AS_IF([test $AR = ":"], [AC_MSG_ERROR([ar command not found.])])

AC_PROG_RANLIB()
AS_IF([test $RANLIB = ":"], [AC_MSG_ERROR([ranlib command not found.])])

AC_CONFIG_AUX_DIR([shared])
AC_PROG_INSTALL()

AC_CONFIG_FILES([config.mk])
AC_OUTPUT

