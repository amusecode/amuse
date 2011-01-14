! Module containing functions useful to construct and calculate
! interpolation functions based on a stellar model.
module interpolate
   use real_kind
   use mesh
   
   implicit none
   ! Do we want to use spline or linear interpolation?
   logical :: use_spline_interpolation = .false.
   
   
contains
   
   
   ! ------------------------------------------------------------------------------
   !  SPLINE_INIT
   !   Taken from Make Me A Star, because the spline code already in the evolution
   !   code is too integrated into the opacity table (sigh...)
   !   This routine is a very slightly modified version of one from the book
   !   Computer Methods for Mathematical Computations, by George Forsythe, Mike
   !   Malcolm, and Cleve Moler. The original copy was obtained from the Netlib
   !   mathematical software file server
   !   (see http://www2.ucsc.edu/cats/sc/software/netlib/) with the command
   !   "send spline from fmm".
   !
   !  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
   !  for a cubic interpolating mmas_spline
   !
   !    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
   !
   !    for  x(i) .le. x .le. x(i+1)
   !
   !  input..
   !
   !    n = the number of data points or knots (n.ge.2)
   !    x = the abscissas of the knots in strictly increasing order
   !    y = the ordinates of the knots
   !
   !  output..
   !
   !    b, c, d  = arrays of spline coefficients as defined above.
   !
   !  using  p  to denote differentiation,
   !
   !    y(i) = s(x(i))
   !    b(i) = sp(x(i))
   !    c(i) = spp(x(i))/2
   !    d(i) = sppp(x(i))/6  (derivative from the right)
   !
   !  the accompanying function subprogram  seval  can be used
   !  to evaluate the mmas_spline.
   ! ------------------------------------------------------------------------------
   !  Input:
   !        N  - The number of grid points (>=2)
   !        X  - The independent variable (usually the mass coordinate)
   !        Y  - The dependent variable
   !  Output:
   !        B  - Linear term in the interpolation function
   !        C  - Quadratic term in the interpolation function
   !        D  - Cubic term in the interpolation function
   ! ------------------------------------------------------------------------------
   subroutine spline_init (n, x, y, b, c, d)
      use real_kind
      
      implicit none
      integer, intent (in) ::  n
      real(double), intent(in) :: x(n), y(n)
      real(double), intent(out) :: b(n), c(n), d(n)
      integer :: nm1, ib, i
      real(double) :: t
      !
      nm1 = n-1
      if ( n < 2 ) return
      if ( n < 3 ) then
         b(1) = (y(2)-y(1))/(x(2)-x(1))
         c(1) = 0.d0
         d(1) = 0.d0
         b(2) = b(1)
         c(2) = 0.d0
         d(2) = 0.d0
         return
      end if
      !
      !  set up tridiagonal system
      !
      !  b = diagonal, d = offdiagonal, c = right hand side.
      !
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.d0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
      end do
      !
      !  end conditions.  third derivatives at  x(1)  and  x(n)
      !  obtained from divided differences
      !
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.d0
      c(n) = 0.d0
      if ( .not. ( n .eq. 3 ) ) then
         c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
         c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
         c(1) = c(1)*d(1)**2/(x(4)-x(1))
         c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      end if
      !
      !  forward elimination
      !
      do i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
      end do
      !
      !  back substitution
      !
      c(n) = c(n)/b(n)
      do ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
      !
      !  c(i) is now the sigma(i) of the text
      !
      !  compute polynomial coefficients
      !
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.d0*c(n))
      do i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.d0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.d0*c(i)
      end do
      c(n) = 3.d0*c(n)
      d(n) = d(n-1)
   end subroutine spline_init
   
   
   
   ! ------------------------------------------------------------------------------
   ! LINEAR_INIT
   ! Set up linear interpolation tables in a way that is compatible with
   ! iptable_eval.
   !  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
   !  for a cubic interpolating mmas_spline
   !
   !    s(x) = y(i) + b(i)*(x-x(i))
   !    c(i) = d(i) = 0.0
   !
   !    for  x(i) .le. x .le. x(i+1)
   !
   ! ------------------------------------------------------------------------------
   !  Input:
   !        N  - The number of grid points (>=2)
   !        X  - The independent variable (usually the mass coordinate)
   !        Y  - The dependent variable
   !  Output:
   !        B  - Linear term in the interpolation function
   !        C  - Quadratic term in the interpolation function (0)
   !        D  - Cubic term in the interpolation function (0)
   ! ------------------------------------------------------------------------------
   subroutine linear_init (n, x, y, b, c, d)
      use real_kind
      
      implicit none
      integer, intent (in) ::  n
      real(double), intent(in) :: x(n), y(n)
      real(double), intent(out) :: b(n), c(n), d(n)
      integer :: nm1, i
      
      nm1 = n-1
      if ( n .lt. 2 ) return
      
      b(:) = 0
      c(:) = 0
      d(:) = 0
      
      do i=2, nm1
         ! Average forward and backward derivatives
         !b(i) = 0.5 * ( (y(i-1) - y(i)) / (x(i-1) - x(i)) +  (y(i) - y(i+1)) / (x(i) - x(i+1)))
         b(i) = (y(i) - y(i+1)) / (x(i) - x(i+1))
      end do
      b(1) = (y(1) - y(2)) / (x(1) - x(2))
      b(n) = (y(n-1) - y(n)) / (x(n-1) - x(n))
   end subroutine linear_init
   
   
   
   ! ------------------------------------------------------------------------------
   ! IPTABLE_EVAL
   !  Again taken from Make Me A Star.
   ! This routine is from the book Computer Methods for Mathematical
   ! Computations, by George Forsythe, Mike Malcolm, and Cleve Moler.
   ! This copy was obtained from the Netlib mathematical software file server
   ! (see http://www2.ucsc.edu/cats/sc/software/netlib/) with the command
   ! "send seval from fmm".
   !
   !  this subroutine evaluates the cubic mmas_spline function
   !
   !    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
   !
   !    where  x(i) .lt. u .lt. x(i+1), using horner's rule
   !
   !  if  u .lt. x(1) then  i = 1  is used.
   !  if  u .ge. x(n) then  i = n  is used.
   !
   !  input..
   !
   !    n = the number of data points
   !    u = the abscissa at which the mmas_spline is to be evaluated
   !    x,y = the arrays of data abscissas and ordinates
   !    b,c,d = arrays of mmas_spline coefficients computed by mmas_spline
   !
   !  if  u  is not in the same interval as the previous call, then a
   !  binary search is performed to determine the proper interval.
   !
   !  This code will also accept the linear interpolation tables calculated by
   !  linear_init() in this module.
   ! ------------------------------------------------------------------------------
   !  Input:
   !        N  - The number of grid points (>=2)
   !        U  - Wanted value for the independent variable
   !        X  - The independent variable (usually the mass coordinate)
   !        Y  - The dependent variable
   !        B  - Linear term in the interpolation function
   !        C  - Quadratic term in the interpolation function (0)
   !        D  - Cubic term in the interpolation function (0)
   !  Return value:
   !        The interpolated value of Y for X = U
   ! ------------------------------------------------------------------------------
   function iptable_eval(n, u, x, y, b, c, d)
      use real_kind
      
      implicit none
      real(double)  :: iptable_eval
      integer, intent(in) :: n
      real(double) , intent(in) :: u, x(n), y(n), b(n), c(n), d(n)
      integer :: i
      integer :: j, k
      real(double) :: dx
      
      !
      !  binary search - NB: Assumes that the array x is stored largest-smallest
      !
      i = 1;
      j = n+1
      do while (j > i+1)
         k = (i+j)/2
         if ( u < x(k) ) then
            i = k 
         else 
            j = k
         end if
      end do
      !
      !  evaluate mmas_spline
      !
      dx = u - x(i)
      iptable_eval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
   end function iptable_eval
   
   
   
   ! ------------------------------------------------------------------------------
   ! IPTABLE_INIT
   !  Prepare data structures for either linear or cubic spline
   !  interpolation, depending on settings.
   ! ------------------------------------------------------------------------------
   !  Input:
   !        N  - The number of grid points (>=2)
   !        X  - The independent variable (usually the mass coordinate)
   !        Y  - The dependent variable
   !  Output:
   !        B  - Linear term in the interpolation function
   !        C  - Quadratic term in the interpolation function (or 0)
   !        D  - Cubic term in the interpolation function (or 0)
   ! ------------------------------------------------------------------------------
   subroutine iptable_init (n, x, y, b, c, d)
      use real_kind
      
      implicit none
      integer, intent (in) ::  n
      real(double), intent(in) :: x(n), y(n)
      real(double), intent(out) :: b(n), c(n), d(n)
      if (use_spline_interpolation) then
         call spline_init (n, x, y, b, c, d)
      else
         call linear_init (n, x, y, b, c, d)
      end if
   end subroutine iptable_init
   
end module interpolate


