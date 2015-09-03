# AX_NETCDF
# ----------------------------------------------------------
# set up for NETCDF
#
AC_DEFUN([AX_NETCDF],[
        AC_ARG_WITH(netcdf,
             AC_HELP_STRING([--with-netcdf=PFX], [Prefix where NETCDF has been installed] ),
             [
                test "$withval" = no && AC_MSG_WARN([netcdf is a required package for some modules])
                test "$withval" = yes || netcdf_prefix="$withval" 
                with_netcdf=yes ],
             [ with_netcdf=yes ] 
        )
    
        AS_IF([test "x$with_netcdf" != xno ],
        [
            #user override
            AS_IF([test "x$NETCDF_LIBS" != x && test "x$NETCDF_FLAGS" != x ],
            [
                have_netcdf=yes
                FOUND_NETCDF="yes"
            ],
            [
                saved_LIBS="$LIBS"
                saved_CXXFLAGS="$CXXFLAGS"
                NETCDF_PREFIX=""
                NETCDF_LIBS=""
                NETCDF_FLAGS=""
                FOUND_NETCDF="no"
                if test x$netcdf_prefix == x; then
                   if test x$PREFIX != x; then
		                netcdf_prefix=$PREFIX 
                   fi  
                fi
                if test x$netcdf_prefix != x; then
                    ac_NETCDF_CFLAGS="-I$netcdf_prefix/include"
                    ac_NETCDF_LDOPTS="-L$netcdf_prefix/lib"
                    
                    
                    save_CFLAGS="$CFLAGS"
                    save_CPPFLAGS="$CPPFLAGS"
                    CFLAGS="$ac_NETCDF_CFLAGS $save_CFLAGS"
                    CPPFLAGS="$ac_NETCDF_CFLAGS $save_CPPFLAGS"
                    AC_CHECK_HEADER(
                        [netcdf.h],
                        [NETCDF_FLAGS="$ac_NETCDF_CFLAGS"
                        NETCDF_PREFIX="$netcdf_prefix"
                        FOUND_NETCDF="yes"],
                        [AC_MSG_WARN([Cannot find headers (netcdf.h) of the library NETCDF in $netcdf_prefix/include.])]
                    )
                    CFLAGS="$save_CFLAGS"
                    CPPFLAGS="$save_CPPFLAGS"
                    
                    save_LIBS="$LIBS"
                    LIBS="$ac_NETCDF_LDOPTS -lnetcdf  $save_LIBS"
                    AC_CHECK_LIB([netcdf], [nc_inq_libvers],
                            [NETCDF_LIBS="$ac_NETCDF_LDOPTS -lnetcdf"],
                            [AC_MSG_WARN([libnetcdf : library missing. (Cannot find symbol nc_inq_libvers) in $netcdf_prefix. Check if libnetcdf3 is installed and if the version is correct])]
                            )
                    LIBS="$save_LIBS"
	fi
	if test x$FOUND_NETCDF != xyes; then
                    AC_CHECK_HEADER(
                        [netcdf.h],
                        [NETCDF_FLAGS=""
                        NETCDF_PREFIX=""
                        FOUND_NETCDF="yes"],
                        [AC_MSG_WARN([Cannot find headers (netcdf.h) of the library NETCDF in $netcdf_prefix/include.])]
                    )
                    AC_CHECK_LIB(
                        [netcdf],
                        [nc_inq_libvers],
                        [NETCDF_LIBS="-lnetcdf"],
                        [
                        FOUND_NETCDF="no"
                        AC_MSG_WARN([libnetcdf : library missing. (Cannot find symbol nc_inq_libvers). Check if libnetcdf is installed and if the version is correct])]
                    )
		fi
                
               
            ])
        ])

        AC_SUBST(FOUND_NETCDF)
        AC_SUBST(NETCDF_FLAGS)
        AC_SUBST(NETCDF_LIBS)
        AC_SUBST(NETCDF_PREFIX)
    ]
)
