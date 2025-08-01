AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-tupan], [1.0])

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Set CPU_COUNT
AX_COUNT_CPUS()
AC_SUBST(CPU_COUNT)

# Find the compiler(s)
AC_PROG_CC()

# Find tools to download and unpack with
AMUSE_DOWNLOAD()
AC_CHECK_TOOL(TAR, tar)
AC_CHECK_TOOL(BUNZIP2, bunzip2)

# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT

