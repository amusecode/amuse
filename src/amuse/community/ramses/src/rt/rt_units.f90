!************************************************************************
subroutine rt_units(scale_Np, scale_Fp)

! Conversion factors from user units into cgs units for RT variables.
!
! scale_np   <=  photon number density scale [# cm-3]
! scale_pf   <=  photon flux density scale [# cm-2 s-1]
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  implicit none
  real(dp)::scale_Np, scale_Fp
  real(dp)::scale_nH, scale_T2, scale_t, scale_v, scale_d, scale_l
!------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  scale_Np = 1.0
  if(cosmo) scale_Np = 1./aexp**3         ! or maybe it should be aexp**4

  scale_Fp = scale_Np  * scale_v

end subroutine rt_units
