
AC_DEFUN([AC_CHECK_PYTHON_DEV],[
    PYTHON_DEV='yes'
    
    AC_ARG_VAR([PYTHONCONFIG], [Python config script to determine python module compile flags (python-config)])
    AS_IF([test "x$PYTHONCONFIG" = "x"], [
	 PYTHONCONFIG=python-config
   	 AC_PATH_PROGS([PYTHONCONFIG], [$PYTHONCONFIG], ['no'])])
    AS_IF([test "x$PYTHONCONFIG" = "xno"], [PYTHON_DEV='no'])
    AC_ARG_VAR([CYTHON], [cython script to compile code in the Cython language])
    AC_ARG_VAR([PYTHONDEV_CFLAGS], [CFLAGS for python compilation])
    AC_ARG_VAR([PYTHONDEV_LDFLAGS], [LDFLAGS for python language])
    AS_IF([test "x$CYTHON" = "x"], [
	CYTHON=cython
    	AC_PATH_PROGS([CYTHON], [$CYTHON], [])])
    AS_IF([test "x$PYTHON_DEV" = "xyes"], [
       	AS_IF([test "x$PYTHONDEV_CFLAGS" = "x"],[ 
        	PYTHONDEV_CFLAGS=`$PYTHONCONFIG --cflags`
        	PYTHONDEV_LDFLAGS=`$PYTHONCONFIG --ldflags`
        	save_CFLAGS="$CFLAGS"
        	save_CPPFLAGS="$CPPFLAGS"
        	CFLAGS="$PYTHONDEV_CFLAGS $save_CFLAGS"
        	CPPFLAGS="$PYTHONDEV_CFLAGS $save_CPPFLAGS"
        	AC_CHECK_HEADER(
            	[Python.h],
            	[],
            	[
            	PYTHON_DEV="no"
            	PYTHONDEV_CFLAGS=''
            	PYTHONDEV_LDFLAGS=''
            	
            	AC_MSG_WARN([Cannot find headers (Python.h)])]
        	) 
        	AS_IF([test "x$PYTHON_DEV" = "xyes"], [
        	
            	save_LDFLAGS="$LDFLAGS"
            	save_LIBS="$LIBS"
            	LIBS=`$PYTHONCONFIG --libs`
            	LDFLAGS="$PYTHONDEV_LDFLAGS $save_LDFLAGS"
            	AC_MSG_CHECKING([if possible to embed python])
            	AC_RUN_IFELSE(
                	[AC_LANG_PROGRAM([[#include "Python.h"]],
                   	[[
                       	Py_Initialize();
                       	PyRun_SimpleString("print 'embedded python',");
                       	Py_Finalize();
                   	]])],
                	[AC_MSG_RESULT([yes])],
                	[PYTHON_DEV="no"
                	PYTHONDEV_CFLAGS=''
                	PYTHONDEV_LDFLAGS=''
                	
                	AC_MSG_RESULT([no])]
            	)
            	LDFLAGS="$save_LDFLAGS"
            	LIBS="$save_LIBS"
        	])
        	CFLAGS="$save_CFLAGS"
        	CPPFLAGS="$save_CPPFLAGS"
	],[])
    ], [])
    
    AC_SUBST(PYTHON_DEV)
    AC_SUBST(CYTHON)
    AC_SUBST(PYTHONDEV_CFLAGS)
    AC_SUBST(PYTHONDEV_LDFLAGS)
])
