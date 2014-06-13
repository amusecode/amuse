module step_functions

contains

pure function step (f, d)
   ! Smoothed step function (=0 for F<=-D, =1 for F>>D)
   ! D is the smoothing length
   use real_kind
   implicit none
   real(double), intent(in) :: f, d
   real(double) :: step
   step = smf (1, f, d)
end function step



pure function pstv (f, d)
   ! Smoothed triangle function (=0 for F<=-D, =F for F>>D)
   ! D is the smoothing length
   use real_kind
   implicit none
   real(double), intent(in) :: f, d
   real(double) :: pstv
   pstv = smf (2, f, d)
end function pstv



pure function smf (iht, f, d)
   ! Smoothed step function (IHT=1) or triangle function (IHT=2)
   use real_kind
   implicit none
   real(double), intent(in) :: f, d
   integer, intent(in) :: iht
   real(double) :: smf
   real(double) :: x, y
   x = f + d  
   if ( x.le.1.0d-150 ) then  ! Make sure that X**2>0 numerically
      smf = 0.0d0
   else
      if(x < 1.d130) then    ! safety net for big numbers !SdM
         y = x*x             ! moved to location where it is needed !SdM
         smf = y/(d*d + y)
      else
         smf = 1.d0
      end if
   end if
   if ( iht == 2 ) smf = x*smf           
end function smf

end module step_functions
