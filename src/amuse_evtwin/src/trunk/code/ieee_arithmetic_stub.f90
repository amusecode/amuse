! Module providing ieee_is_nan, to be used as a work-around if ieee_arithmetic
! is not provided by the compiler using some pre-processor trickery.
module ieee_arithmetic_stub

#ifdef ieee_arithmetic

contains

   logical elemental function ieee_is_nan(x)
      use real_kind
      implicit none
      real(double), intent(in) :: x
      ieee_is_nan = isnan(x)
   end function ieee_is_nan
#endif

end module ieee_arithmetic_stub
