subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  use poisson_commons
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------
  real(dp)::v200,c

  v200 = gravity_params(5) ! Virial velocity of the halo (km.s-1)
  c = gravity_params(6)    ! Concentration of the halo
  
  ! scale_d converts mass density from user units into g/cc
  scale_d = rhoc*200./3.0*c*c*c/(log(1d0+c)-c/(1d0+c))

  ! scale_t converts time from user units into seconds
  scale_t = 1.0/sqrt(6.67d-8*scale_d)

  ! scale_l converts distance from user units into cm
  scale_l = v200*1d5 * scale_t / sqrt(4d0*ACOS(-1d0)*(log(1d0+c)-c/(1d0+c))/c)

  ! scale_v convert velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

end subroutine units
