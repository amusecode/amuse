subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! scale_d converts mass density from user units into g/cc
  scale_d = 1.0

  ! scale_t converts time from user units into seconds
  scale_t = 1.0

  ! scale_l converts distance from user units into cm
  scale_l = 1.0

  ! scale_v convert velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = scale_d

end subroutine units
