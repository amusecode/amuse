      module compare_floats
      contains
      
! Floating point number comparisons can be rather unreliable and tend to
!  break at unexpected times. 
! Ideas here are inspired by this page:
!  http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
      logical function feq(f1, f2)
      implicit none
      double precision :: f1, f2
      double precision, parameter :: abserr = 1d-300;
      double precision, parameter :: relerr = 1d-24;
      !logical :: t1, t2, t3;
      !
      !t1 = f1 == f2;
      !t2 = dabs(f1-f2) < abserr;
      !t3 = dabs( (f1-f2)/max(f1,f2) ) < relerr;
      !if (t1 /= t2 .or. t1 /= t3 .or. t2/=t3) then
      !   print *, t1, t2, t3, f1-f2
      !end if
      
      ! Check absolute difference
      if (f1 == f2 .or. abs(f1-f2) < abserr) then
        feq = .true.;
        return
      end if
       
      ! Check relative difference
      if ( abs( (f1-f2)/max(f1,f2) ) < relerr) then
        feq = .true.;
        return
      end if

      ! Difference too large - return          
      feq = .false.;
      return
       
      end function
      
      logical function nfeq(f1, f2)
      implicit none
      double precision :: f1, f2;
      
      nfeq = .not. feq(f1, f2);
      return;
      end function
      
      end module
