! Various subroutines for calculating quantities in Roche geometry
module roche
   use real_kind

contains

   ! Calculate the scaled (a = 1) x coordinate of the binary component Jstar.
   ! Internally, all calculations assume m1 > m2, q = m2/m1 < 1 and x1 > 0.
   ! This means that when the mass ratio inverts, the coordinate system changes.
   ! This can be handled entirely automatically, as long as the arguments passed
   ! in are consistent and all functions perform the required coordinate
   ! transformation.
   elemental function calc_scaled_x(m1, m2, Jstar)
   implicit none
   real(double), intent(in) :: m1, m2
   integer, intent(in) :: Jstar
   real(double) :: calc_scaled_x
   real(double) :: q
   integer :: n

   q = m2 / m1
   n = Jstar
   if (q > 1.0d0) then
      q = m1 / m2
      n = 3 - Jstar
   end if

   if (n == 1) then
      calc_scaled_x = q / (1.0d0 + q)
   else
      calc_scaled_x = -1.0d0 / (1.0d0 + q)
   end if
   end function calc_scaled_x

   elemental function calc_scaled_potential(m1, m2, x, y, z)
   implicit none
   real(double), intent(in) :: m1, m2, x, y, z
   real(double) :: calc_scaled_potential
   real(double) :: q, xx, qrat, yz2

   q = m2 / m1
   xx = x
   if (q > 1.0d0) then
      q = 1.0d0 / q
      xx = 1.0d0 - x
   end if

   qrat = 1.0d0 / (1.0d0 + q)
   yz2 = y**2 + z**2
   calc_scaled_potential = qrat * (1.0d0/sqrt((x - q*qrat)**2+yz2) + q / sqrt((x+qrat)**2+yz2)) + 0.5d0*(x**2+y**2)
   end function calc_scaled_potential

   elemental function calc_scaled_xl1(m1, m2)
   implicit none
   real(double), intent(in) :: m1, m2
   real(double) :: calc_scaled_xl1
   real(double) :: q,q1,q2,q3,q4,q5,q6,xl1
   integer :: ij

   q = m1 / m2
   if ( q < 1.0d0 ) q = 1.0d0/q
   q5 = 1.0d0 + q
   xl1 = (q/(3.0d0*q5))**(1.0d0 / 3.0d0)
   q3 = 3.0d0 + q
   q2 = q + q
   q4 = 3.0d0 + q2
   q1 = 2.0d0 + q
   q6 = 1.0d0/q5
  ! Newton-Raphson iteration for L1
   do ij = 1, 4
     xl1 = xl1 + &
          (q - xl1*(q2 - xl1*(q - xl1*(q3 - xl1*(q4 - xl1*q5)))))&
          /(q2 - xl1*(q2 - xl1*(3.0d0*q3 - xl1*(4.0d0*q4 - xl1*5.0d0*q5))))
   end do
   calc_scaled_xl1 = q6 - xl1
   end function calc_scaled_xl1

   elemental function calc_scaled_xl2(m1, m2)
   implicit none
   real(double), intent(in) :: m1, m2
   real(double) :: calc_scaled_xl2
   real(double) :: q,q1,q2,q3,q4,q5,q6,xl1,xl2
   integer :: ij

   q = m1 / m2
   if ( q < 1.0d0 ) q = 1.0d0/q
   q5 = 1.0d0 + q
   xl1 = (q/(3.0d0*q5))**(1.0d0 / 3.0d0)
   xl2 = 2.0d0 - q*xl1/q5
   q3 = 3.0d0 + q
   q2 = q + q
   q4 = 3.0d0 + q2
   q1 = 2.0d0 + q
   q6 = 1.0d0/q5
  ! Newton-Raphson iteration for L2
   do ij = 1, 4
     xl2 = xl2 + &
          (q - xl2*(q2 - xl2*(q1 - xl2*(q3 - xl2*(q4 - xl2*q5)))))&
          /(q2 - xl2*(2.0d0*q1 - xl2*(3.0d0*q3 - xl2*(4.0d0*q4 - xl2*5.0d0*q5))))
   end do
   calc_scaled_xl2 = q6 - xl2
   end function calc_scaled_xl2

   ! Calculate the correction factor for the cross section of the mass flow at level surface xx (dimensionless, L1=0, L2=1)
   elemental function stream_cross_section_primary(q, xx)
      use real_kind
      implicit none
      real(double), intent(in) :: q, xx
      real(double) :: stream_cross_section_primary
      real(double) :: a, b, c, d

      a = 4.50d-4 + (3.68d-4 + 3.15d-4 * q) * q
      if (q < 0.3d0) then
         b = 2.32d0 + (-1.31d0 + 2.87d0 * q) * q
      else
         b = 2.21d0 + 0.238d0 * q
      end if
      c = q * (1.22d0 + q*(-1.22d0 + q*0.564d0))
      d = 0.323d0 * q**0.230d0 - 7.39d-2

      stream_cross_section_primary = (c * xx + a)**b + d * xx
   end function stream_cross_section_primary

   ! Calculate the correction factor for the cross section of the mass flow at level surface xx (dimensionless, L1=0, L2=1)
   elemental function stream_cross_section_secondary(q, xx)
      use real_kind
      implicit none
      real(double), intent(in) :: q, xx
      real(double) :: stream_cross_section_secondary
      real(double) :: a, b, c, d

      a = 9.53d-4 * q**0.501d0 + 1.65d-4
      b = 3.25d-1 * q**0.326d0 + 2.03d0
      c = 6.33d-1 * q**0.271d0 - 7.18d-2
      d = 2.30d-1 * q**0.361d0 + 1.89d-2

      stream_cross_section_secondary = (c * xx + a)**b + d * xx
   end function stream_cross_section_secondary

   !> \brief Compute the Roche-lobe radius Rl/a
   !!
   !! \param x  q^(1/3)
   elemental function rlobe(x)
      use real_kind
      use constants

      implicit none
      real(double), intent(in) :: x
      real(double) :: rlobe,x2

      x2 = x**2
      rlobe = 0.49_dbl*x2 / (0.6_dbl*x2 + log(1.0_dbl + x))
   end function rlobe

end module
