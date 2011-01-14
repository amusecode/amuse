module neutrinos

   use real_kind
   
   implicit none
   ! Neutrino-loss rates from Itoh et al (1983-1992)
   real(double) :: cnu(60,41,2)
   
contains

   subroutine get_neutrino_rate(tf, frho, abund, en)
      use real_kind
      use eostate_types
      
      implicit none
      real(double), intent(in) :: tf, frho
      type(abundance), intent(in) :: abund
      real(double), intent(out) :: en
      real(double) :: tt, tu, rn, ru, enp, enb
      integer :: it, ir
      en = 0.0d0
      if ( tf >= 7.0d0 .and. frho >= 0.0d0 ) then
         tt = 20.0d0*(tf - 6.95d0)
         it = max0(1, min0(59, int(tt)))
         tt = tt - it
         tu = 1.0d0 - tt
         rn = 4.0d0*(frho + 0.25d0)
         ir = max0(1, min0(40, int(rn)))
         rn = rn - ir
         ru = 1.0d0 - rn
         enp = tt*(rn*cnu(it + 1, ir + 1, 1) + ru*cnu(it + 1, ir, 1))  &
              + tu*(rn*cnu(it, ir + 1, 1) + ru*cnu(it, ir, 1))
         enb = tt*(rn*cnu(it + 1, ir + 1, 2) + ru*cnu(it + 1, ir, 2))  &
              + tu*(rn*cnu(it, ir + 1, 2) + ru*cnu(it, ir, 2))
         en = -10.0d0**enp - abund%nzz*10.0d0**enb
      end if
   end subroutine get_neutrino_rate
   
end module neutrinos
