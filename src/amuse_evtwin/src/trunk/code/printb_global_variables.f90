module printb_global_variables
   ! Module for exchanging variables between subroutines within PRINTB only
   use real_kind

   implicit none
   real(double) ::  dmtr, perc, rcz, drcz, dmt, tet, RcoreH, RcoreHe
   real(double) :: dmsw, oms, fr, fl, sdc, sdm, sds, stc, stm, sts, psurf, pcntr, vk2
   real(double) :: raf, dphi, ephi, f1, df, htot, aj, dl, tkh, mex(12), ecc

   ! Meshpoint were maximum temperature is reached
   integer :: im_tmax

   ! Pressures at the edge of the H and He exhausted cores
   real(double) :: ph
   real(double) :: phe
end module printb_global_variables






