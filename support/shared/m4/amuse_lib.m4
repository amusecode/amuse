# Helper macros for detecting AMUSE libraries.
#
# See https://stackoverflow.com/questions/10220946/pkg-check-modules-considered-harmful
# for some background on why we're not using PKG_CHECK_MODULES directly. Do note that
# PKG_CHECK_MODULES checks for pkg-config and gives a useful error these days, so that
# that page is partially outdated.

# AMUSE_LIB(prefix, module, library, function)
#
# Searches for an AMUSE library and sets ${prefix}_CFLAGS, ${prefix}_LIBS to the
# appropriate values if it is found.
#
# prefix: prefix for the variables to be set
# module: name of the pkg-config module to seach for if needed
# library: name of the library to be searched for
# function: name of a function in the library to use for the link check
AC_DEFUN([AMUSE_LIB], [
    amuse_save_LIBS="$LIBS"

    # Search for the library, first directly then fall back to pkg-config
    AC_SEARCH_LIBS([$4], [$3], [
        $1_LIBS="$LIBS"
        $1_CFLAGS=""
    ], [
        PKG_CHECK_MODULES([$1], [$2])
    ])

    LIBS="$amuse_save_LIBS"

    AC_SUBST([$1][_CFLAGS])
    AC_SUBST([$1][_LIBS])
])


# AMUSE_LIB_STOPCOND()
#
# Searches for the AMUSE stopping conditions library and sets STOPCOND_CFLAGS
# and STOPCOND_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_STOPCOND], [
    AMUSE_LIB([STOPCOND], [stopcond], [stopcond], [is_condition_enabled])
])


# AMUSE_LIB_STOPCONDMPI()
#
# Searches for the AMUSE stopping conditions library and sets STOPCONDMPI_CFLAGS
# and STOPCONDMPI_LIBS to the appropriate values if it is found.
#
# We use a different function here, to avoid getting a cached value from
# AMUSE_LIB_STOPCOND if both are used. Seems like it indexes by function name.
AC_DEFUN([AMUSE_LIB_STOPCONDMPI], [
    AMUSE_LIB([STOPCONDMPI], [stopcondmpi], [stopcondmpi], [get_set_conditions_])
])


# AMUSE_LIB_AMUSE_MPI()
#
# Searches for the AMUSE MPI helper library and sets AMUSE_MPI_CFLAGS and
# AMUSE_MPI_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_AMUSE_MPI], [
    AMUSE_LIB([AMUSE_MPI], [amuse_mpi], [amuse_mpi], [get_comm_world])
])


# AMUSE_LIB_FORSOCKETS()
#
# Searches for the AMUSE forsockets library and sets FORSOCKETS_CFLAGS and
# FORSOCKETS_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_FORSOCKETS], [
    AMUSE_LIB([FORSOCKETS], [forsockets], [forsockets], [forsockets_close])
])


# AMUSE_LIB_SIMPLE_HASH()
#
# Searches for the AMUSE simple hash library and sets SIMPLE_HASH_CFLAGS and
# SIMPLE_HASH_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_SIMPLE_HASH], [
    AMUSE_LIB([SIMPLE_HASH], [simple_hash], [simple_hash], [init_hash])
])


# AMUSE_LIB_G6LIB()
#
# Searches for the g6 library and sets G6LIB_CFLAGS and G6LIB_LIBS to
# the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_G6LIB], [
    AMUSE_LIB([G6LIB], [g6lib], [g6], [g6_npipes])
])


# AMUSE_LIB_SAPPORO_LIGHT()
#
# Searches for the Sapporo light library and sets SAPPORO_LIGHT_CFLAGS and
# SAPPORO_LIGHT_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_SAPPORO_LIGHT], [
    AMUSE_LIB([SAPPORO_LIGHT], [sapporo_light], [sapporo], [get_device_count])
])

