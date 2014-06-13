module cuberoot

contains

!> \brief Compute the cubic root
!!
!! \param x  input variable

elemental function cbrt(x)
   use real_kind
   use constants

   implicit none
   real(double), intent(in) :: x
   real(double) :: cbrt

   cbrt = abs(x)**c3rd
end function cbrt

end module cuberoot
