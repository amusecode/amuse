# Helper macro for installing into Conda envs and virtualenvs.

# AMUSE_VENV()
#
# Detect the active Conda env or virtualenv and extend LDFLAGS and
# PKG_CONFIG_PATH to point to it. This then makes it possible to detect the AMUSE
# libraries as installed in the environment and link to them.
#
AC_DEFUN([AMUSE_VENV], [
    AS_IF([test "x$VIRTUAL_ENV" != x], [
        LDFLAGS="$LDFLAGS -L${VIRTUAL_ENV} -Wl,-rpath ${VIRTUAL_ENV}/lib"
        PKG_CONFIG_PATH="$VIRTUAL_ENV/lib/pkgconfig:$PKG_CONFIG_PATH"
    ])

    AS_IF([test "x$CONDA_PREFIX" != x], [
        LDFLAGS="$LDFLAGS -L${CONDA_PREFIX} -Wl,-rpath ${CONDA_PREFIX}/lib"
        # Conda pkg-config includes this already, but in case we have one from
        # the system...
        PKG_CONFIG_PATH="$PKG_CONFIG_PATH:${CONDA_PREFIX}/lib/pkgconfig"
    ])
    # Needs to be exported or the PKG_CHECK_MODULES macro won't see it
    export PKG_CONFIG_PATH
])

