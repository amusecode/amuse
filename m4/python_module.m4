# AC_CHECK_PYTHON_MODULE Autoconf macro, revision 0
#
# Copyright (C) 2008 Stephan Peijnik
#
# Copying and distribution of this file, with or without modification,         
# are permitted in any medium without royalty provided the copyright notice    
# and this notice are preserved.  
#
# Sources:
#   http://blog.sp.or.at/2008/08/31/autoconf-and-python-checking-for-modules/
#

AC_DEFUN([AC_CHECK_PYTHON_MODULE],[
    # AC_CHECK_PYTHON_MODULE(MODULE_NAME [,VERSION_VARIABLE])

    # the python module name
    MODULE_NAME=$1
    # the python variable that contains the module's version
    # If this is not set the version will not be retrieved from the module.
    # Example: __version__.
    VERSION_VARIABLE=$1.$2

    # check for the python binary defined in $PYTHON
    # fall back to "python"
    if test -z $PYTHON;
    then
        PYTHON="python"
    fi

    AC_MSG_CHECKING(for python module $MODULE_NAME)

    if test -z "$2"
    then
      $PYTHON -c "import $MODULE_NAME" 2>/dev/null
      if test $? -eq 0
      then
        eval PYTHON_${MODULE_NAME}=1
        AC_MSG_RESULT(found)
       else
        eval PYTHON_${MODULE_NAME}=0
        AC_MSG_RESULT(not found)
      fi
      AC_SUBST(PYTHON_$1)
    else
      VERSION=`$PYTHON -c "import $MODULE_NAME; print ($VERSION_VARIABLE)" 2>/dev/null`
      if test $? -eq 0
      then
        eval PYTHON_${MODULE_NAME}_VERSION=$VERSION
        eval PYTHON_${MODULE_NAME}=1

        AC_MSG_RESULT([found ($VERSION)])
      else
        eval PYTHON_${MODULE_NAME}=0
        eval PYTHON_${MODULE_NAME}_VERSION=0

        AC_MSG_RESULT(not found)
      fi
      AC_SUBST(PYTHON_$1_VERSION)
      AC_SUBST(PYTHON_$1)
    fi
])
