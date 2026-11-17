AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-metisse], [1.0])

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Set the worker language
AC_LANG([Fortran])

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


# Find tools to download and unpack with
AMUSE_DOWNLOAD()
AC_CHECK_TOOL(TAR, tar)
AC_CHECK_TOOL(UNZIP, unzip)
AC_CHECK_TOOL(GUNZIP, gunzip)
AC_CHECK_TOOL(BUNZIP2, bunzip2)
AC_CHECK_TOOL(UNXZ, unxz)
AC_CHECK_TOOL(PERL, perl)


# Find AMUSE libraries
AMUSE_LIB_STOPCOND()
AMUSE_LIB_STOPCONDMPI()
AMUSE_LIB_AMUSE_MPI()
AMUSE_LIB_FORSOCKETS()
AMUSE_LIB_SIMPLE_HASH()
AMUSE_LIB_G6LIB()
AMUSE_LIB_SAPPORO_LIGHT()


# Find external dependencies
##### These select a language to use to detect a library #####
##### Wrap the library tests below in these depending on #####
##### from which language the library will be used by the code #####
AC_LANG_PUSH([Fortran])
AC_LANG_POP([Fortran])

##### Wrap these into language macros as needed #####
##### MPI is always needed to build the worker #####
AX_MPI()


# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT
