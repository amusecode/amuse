AC_DEFUN([AX_FC_WORKS],[
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([whether the Fortran 90 compiler ($FC $FCFLAGS $LDFLAGS) works])
AC_LINK_IFELSE([
    AC_LANG_SOURCE([
        program conftest
        integer, dimension(10) :: n
        end
    ])
],[
    ax_cv_prog_fc_works="yes"
    AC_MSG_RESULT([$ax_cv_prog_fc_works])
],[
    ax_cv_prog_fc_works="no"
    AC_MSG_WARN([installation or configuration problem: Fortran 90 compiler cannot create executables.])
])
# The intel compiler sometimes generates these work.pc and .pcl files
rm -f work.pc work.pcl
AC_LANG_POP(Fortran)
])


AC_DEFUN([AX_FC_ISO_C_BINDING],[
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([if the fortran compiler supports iso c binding])
AC_LINK_IFELSE([
    AC_LANG_SOURCE([
        program conftest
        use ISO_C_BINDING
        integer, dimension(10) :: n
        end
    ])
],[
    ax_cv_fc_iso_c_bindings="yes"
    AC_MSG_RESULT([$ax_cv_fc_iso_c_bindings])
],[
    ax_cv_fc_iso_c_bindings="no"
    AC_MSG_RESULT([no, fortran codes and sockets or embedding will not work.])
])
# The intel compiler sometimes generates these work.pc and .pcl files
rm -f work.pc work.pcl
AC_LANG_POP(Fortran)
])
