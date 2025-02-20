#!/bin/sh

PACKAGE="$1"

if test -z "${PACKAGE}" ; then
    echo "  Incorrect usage, please specify a package name as the first argument."
    exit 1
fi

if test -n "${VIRTUAL_ENV}${CONDA_DEFAULT_ENV}" ; then
    pip_package_line=$(pip list | grep "${PACKAGE} ")
    if [ "a${pip_package_line}" = "a" ] ; then
        printf '%s\n' "  Package ${PACKAGE} was not installed, not uninstalling."
        exit 0
    fi
fi

if test -n "${VIRTUAL_ENV}" ; then
    # grep -q stops reading when a match is found, which then crashes pip list, so we
    # redirect instead. With conda list we can use -q as usual.
    if [ "a${pip_package_line}" != "a" ] ; then
        pip uninstall -y ${PACKAGE}
    fi
fi

if test -n "${CONDA_DEFAULT_ENV}" ; then
    conda_package_line=$(conda list | grep "^${PACKAGE} ")
    if printf '%s' "${conda_package_line}" | grep "<develop>" >/dev/null 2>&1 ; then
        # Conda is showing a pip-installed develop package. However, there may be a
        # conda package hidden underneath, so we're going to try to conda uninstall
        # anyway to fix that if needed. This will fail if everything is as it should
        # be, so we make sure to inhibit the error message to reduce confusion.
        EXPECTING_FAIL=true
    fi

    if printf '%s' "${conda_package_line}" | grep -v "pypi$" >/dev/null 2>&1 ; then
        TMPOUT=$(mktemp)
        conda remove -y ${PACKAGE} >${TMPOUT} 2>&1
        err="$?"
        # If it failed and we were not expecting it to, print the output
        if test "$err" != 0 -a -z "${EXPECTING_FAIL}" ; then
            cat ${TMPOUT}
        fi
        rm ${TMPOUT}
    fi

    if [ "a${pip_package_line}" != "a" ] ; then
        pip uninstall -y ${PACKAGE}
    fi
fi

