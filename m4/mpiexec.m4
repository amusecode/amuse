AC_DEFUN(FIND_MPIEXEC, [
	AC_ARG_VAR([MPIEXEC], [mpiexec or mpirun location])
	AC_CHECK_PROGS(MPIEXEC, mpiexec mpirun)
	AC_SUBST(MPIEXEC)
])
