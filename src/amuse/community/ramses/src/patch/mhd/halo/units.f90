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
  real(dp)::v200,c,r200,h

  v200 = gravity_params(1) ! Virial velocity of the halo (km.s-1)
  c = gravity_params(2)    ! Concentration of the halo
  
  ! NFW scale radius
  h = 0.7
  r200 = v200*1d5/sqrt(2d0*twopi*6.67d-8/3.*200.*rhoc*h**2)

  ! scale_d converts mass density from user units into g/cc
  scale_d = rhoc*h**2*200./3.0*c*c*c/(log(1d0+c)-c/(1d0+c))

  ! scale_t converts time from user units into seconds
  scale_t = 1.0/sqrt(6.67d-8*scale_d)

  ! scale_l converts distance from user units into cm
  scale_l = r200/c

  ! scale_v convert velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

end subroutine units
