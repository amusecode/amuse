
AC_DEFUN([AX_GFORTRAN_OPTION], [
if test "x$GCC" = "xyes"; then
	AC_MSG_CHECKING([if gfortran accepts $2 option])
   	if AC_TRY_COMMAND($FC $2) >/dev/null 2>&1; then
   		$1=$3
   	        AC_MSG_RESULT([yes])
   	else
   		$1=$4
   		AC_MSG_RESULT([no])
   	fi
else
	unset $1
        AC_MSG_RESULT([sorry, no gcc available])
fi
])



AC_DEFUN([AX_GFORTRAN_VERSION], [
  GFORTRAN_VERSION=""
  AX_GFORTRAN_OPTION(ax_gcc_version_option, [-dumpversion],
     [yes],
     [no])
  AS_IF([test "x$GCC" = "xyes"],[
  
    AS_IF([test "x$ax_gcc_version_option" != "xno"],[
      
	    AC_MSG_CHECKING([checking gfortran version])
        ax_cv_gcc_version="`$FC -v 2>&1 |  grep gcc\ version | cut -d\  -f3`"
        AS_IF([test "x$ax_cv_gcc_version" = "x"],[
          ax_cv_gcc_version=""
          AC_MSG_RESULT([could not determine version])
        ], [
          AC_MSG_RESULT([$ax_cv_gcc_version])])
      
      GFORTRAN_VERSION=$ax_cv_gcc_version
    ])
  ])
  AC_SUBST([GFORTRAN_VERSION])
])


#
# AC_FIND_LIB(LIBRARY, FUNCTION, LIST-OF-DIRECTORIES [, ACTION-IF-FOUND
#		[, ACTION-IF-NOT-FOUND [, OTHER-LIBRARIES ]]])
AC_DEFUN(AC_FIND_LIB,
    [AC_MSG_CHECKING([directories for -l$1])
    #
    # A lot of this is taken from AC_CHECK_LIB.  Note that we always check
    # the "no directory" case first.
    #
    ac_lib_var=`echo $1['_']$2 | tr './+\055' '__p_'`
    ac_save_LIBS="$LIBS"
    AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
    [for dir in "" $3
    do
        ac_cache_save_LIBS="$LIBS"
        if test "X$dir" = "X"; then
            LIBS="$LIBS -l$1 $6"
        else
            LIBS="$LIBS -L$dir -l$1 $6"
        fi
        AC_TRY_LINK(
    ifelse([$2], [main], , # Avoid conflicting decl of main.
    [/* Override any gcc2 internal prototype to avoid an error.  */
    ]
    [/* We use char because int might match the return type of a gcc2
        builtin and then its argument prototype would still apply.  */
    char $2();
    ]),
        [$2()],
        if test "X$dir" = "X"; then
            eval "ac_cv_lib_$ac_lib_var=yes"
        else
            eval "ac_cv_lib_$ac_lib_var=$dir"
        fi
        break,
        eval "ac_cv_lib_$ac_lib_var=no")
        LIBS="$ac_cache_save_LIBS"
    done
    ])#
    LIBS="$ac_save_LIBS"
    if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" != no"; then
      eval "dir=\"`echo '$ac_cv_lib_'$ac_lib_var`\""
      if test "$dir" = "yes"; then
        AC_MSG_RESULT([found])
      else
        AC_MSG_RESULT([found in $dir])
      fi
        if test "$dir" = "yes"; then
          LIBS="$LIBS -l$1"
        else 
          LIBS="$LIBS -L$dir -l$1"
        fi
        $4
    else
      AC_MSG_RESULT(not found)
      $5
    fi
])
