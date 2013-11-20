!************************************************************************
SUBROUTINE rt_read_hydro_params(nml_ok)

! 'Interface' for reading any additional parameters from .nml file.
! This routine can be overridden by any patch to include more parameters.
!------------------------------------------------------------------------
  implicit none
  logical::nml_ok
!------------------------------------------------------------------------
  call read_rt_params(nml_ok)

END SUBROUTINE rt_read_hydro_params




