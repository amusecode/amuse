AC_DEFUN(FIND_MPIEXEC, [


	AC_ARG_ENABLE(mpiexec,
	    [AS_HELP_STRING([--enable-mpiexec],[Enable MPI exec usage by distributed AMUSE, by default enabled])],
	    [WITH_MPIEXEC=no
	    AS_IF([test "x$enable_mpiexec" != "xno"], [
		    WITH_MPIEXEC=yes
	    ])
	    ],
	    [WITH_MPIEXEC=yes]
	)

	AC_SUBST(WITH_MPIEXEC)


	AC_ARG_VAR([MPIEXEC], [mpiexec or mpirun location])

	AS_IF([test x"$WITH_MPIEXEC" != xno],[
		AC_PATH_PROGS(MPIEXEC, $MPIEXEC mpiexec mpirun)
		AC_SUBST(MPIEXEC)
	], [
		MPIEXEC=""
		AC_SUBST(MPIEXEC)
	])
])
