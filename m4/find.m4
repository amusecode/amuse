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
