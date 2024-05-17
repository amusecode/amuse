#!/bin/sh

PACKAGE="$1"

if test -z ${PACKAGE} ; then
    echo "Incorrect usage, please specify a package name as the first argument."
    exit 1
fi

if test -n ${VIRTUAL_ENV} ; then
    # grep -q stops reading when a match is found, which then crashes pip list, so we
    # redirect instead. With conda list we can use -q as usual.
    if pip list | grep "^${PACKAGE}" >/dev/null 2>&1 ; then
        pip uninstall -y ${PACKAGE}
    fi
fi

if test -n ${CONDA_DEFAULT_ENV} ; then
    if conda list | grep "^{PACKAGE}" | grep -q "<develop>" ; then
        # Conda is showing a pip-installed develop package. However, there may be a
        # conda package hidden underneath, so we're going to try to conda uninstall
        # anyway to fix that if needed. This will fail if everything is as it should
        # be, so we make sure to output an error message to reduce confusion.
        EXPECTING_FAIL=true
    fi

    if conda list | grep -q "^${PACKAGE}" ; then
        TMPOUT=$(mktemp)
        conda uninstall -y ${PACKAGE} >${TMPOUT} 2>&1
        # If it failed and we were not expecting it to, print the output
        if test $? != 0 -a -z ${EXPECTING_FAIL} ; then
            cat ${TMPOUT}
        fi
        rm ${TMPOUT}
    fi

    if pip list | grep "^${PACKAGE}" >/dev/null 2>&1 ; then
        pip uninstall -y ${PACKAGE}
    fi
fi

