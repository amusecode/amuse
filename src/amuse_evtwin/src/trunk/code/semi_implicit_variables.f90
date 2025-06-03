module semi_implicit_variables
   use real_kind
   implicit none

   ! Expansion coefficients for the radius in spherical harmonics, for rotation and tides
   real(double), allocatable, save :: a_rot2(:,:), a_tide(:, :, :)

   ! Ratio between the average radius r0 and the volume radius rv
   real(double), allocatable, save :: r0_over_rv(:, :)

   ! fp and ft correction terms for distorted stars
   real(double), allocatable, save :: fp_star(:, :)
   real(double), allocatable, save :: ft_star(:, :)

   ! Vertical component of the meridional circulation
   real(double), allocatable, save :: Vmc2_star(:, :)

   ! Difference in rotation rate, ~domega / dk
   real(double), allocatable, save :: diff_omega(:, :)

   ! Total mass loss rate due to RLOF
   real(double), save :: mdot_rlof(2)
   real(double), save :: mdot_rlof0(2)
end module semi_implicit_variables

