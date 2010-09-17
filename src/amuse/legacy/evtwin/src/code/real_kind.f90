!
!> module real_kind:
!! 
!! Contains the integers double and dbl, which shall be used in (almost) all routines
!! to provide the kind of a (currently double-precision) real variable type.
!! 
!! Variables can be declared using "real(double) :: "; constants can be defined as 
!! e.g. "x = 3.0_dbl".
!< 

module real_kind
   implicit none
   integer, parameter :: double = selected_real_kind(15,307)
   integer, parameter :: dbl = selected_real_kind(15,307)

   ! Single precision, should not normally be used
   integer, parameter :: single = selected_real_kind(6, 37)
end module real_kind


