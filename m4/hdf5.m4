# AX_HDF5
# ----------------------------------------------------------
# set up for HDF5
#
AC_DEFUN([AX_HDF5],[
        AC_ARG_WITH(hdf5,
             AC_HELP_STRING([--with-hdf5=PFX], [Prefix where HDF5 has been installed] ),
             [
                test "$withval" = no && AC_MSG_WARN([hdf5 is a required package for some modules])
                test "$withval" = yes || hdf5_prefix="$withval" 
                with_hdf5=yes ],
             [ with_hdf5=yes ] 
        )
    
        AS_IF([test "x$with_hdf5" != xno ],
        [
            #user override
            AS_IF([test "x$HDF5_LIBS" != x && test "x$HDF5_FLAGS" != x ],
            [
                have_hdf5=yes
                FOUND_HDF5="yes"
            ],
            [
                saved_LIBS="$LIBS"
                saved_CXXFLAGS="$CXXFLAGS"
                HDF5_LIBS=""
                HDF5_FLAGS=""
                FOUND_HDF5="no"
                if test x$hdf5_prefix == x; then
                   if test x$PREFIX != x; then
		                hdf5_prefix=$PREFIX 
                   fi  
                fi
                if test x$hdf5_prefix != x; then
                    ac_HDF5_CFLAGS="-I$hdf5_prefix/include"
                    ac_HDF5_LDOPTS="-L$hdf5_prefix/lib"
                    
                    
                    save_CFLAGS="$CFLAGS"
                    save_CPPFLAGS="$CPPFLAGS"
                    CFLAGS="$ac_HDF5_CFLAGS $save_CFLAGS"
                    CPPFLAGS="$ac_HDF5_CFLAGS $save_CPPFLAGS"
                    AC_CHECK_HEADER(
                        [hdf5.h],
                        [HDF5_FLAGS="$ac_HDF5_CFLAGS"
                        FOUND_HDF5="yes"],
                        [AC_MSG_WARN([Cannot find headers (hdf5.h) of the library HDF5 in $hdf5_prefix/include.])]
                    )
                    CFLAGS="$save_CFLAGS"
                    CPPFLAGS="$save_CPPFLAGS"
                    
                    save_LIBS="$LIBS"
                    LIBS="$ac_HDF5_LDOPTS -lhdf5  $save_LIBS"
                    AC_CHECK_LIB([hdf5], [H5_init_library],
                            [HDF5_LIBS="$ac_HDF5_LDOPTS -lhdf5"],
                            [AC_MSG_WARN([libhdf5 : library missing. (Cannot find symbol H5_init_library) in $hdf5_prefix. Check if libhdf53 is installed and if the version is correct])]
                            )
                    LIBS="$save_LIBS"
	fi
	if test x$FOUND_HDF5 != xyes; then
                    AC_CHECK_HEADER(
                        [hdf5.h],
                        [HDF5_FLAGS=""
                        FOUND_HDF5="yes"],
                        [AC_MSG_WARN([Cannot find headers (hdf5.h) of the library HDF5 in $hdf5_prefix/include.])]
                    )
                    AC_CHECK_LIB(
                        [hdf5],
                        [H5_init_library],
                        [HDF5_LIBS="-lhdf5"],
                        [
                        FOUND_HDF5="no"
                        AC_MSG_WARN([libhdf5 : library missing. (Cannot find symbol H5_init_library). Check if libhdf5 is installed and if the version is correct])]
                    )
		fi
                
               
            ])
        ])

        AC_SUBST(FOUND_HDF5)
        AC_SUBST(HDF5_FLAGS)
        AC_SUBST(HDF5_LIBS)
    ]
)
