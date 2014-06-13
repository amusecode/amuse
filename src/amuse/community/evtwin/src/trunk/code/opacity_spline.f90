module opacity_spline


contains

! Calculate the coefficients of a 1-D cubic spline:
! Forsythe, Malcolm, Moler, Computer Methods for Mathematical
! Computations, Prentice-Hall, 1977, p.76
pure subroutine spline ( k, x, f )
   use real_kind

   implicit none
   integer, intent(in) :: k
   real(double), intent(in) :: x(*)
   real(double), intent(out) :: f(4,*)

   integer :: i,ib
   real(double) :: t

   f(2:4,k) = 0.0d0
   ! Set up a tridiagonal system for A*y=B where y(i) are the second
   ! derivatives at the knots.
   ! f(2,i) are the diagonal elements of A
   ! f(4,i) are the off-diagonal elements of A
   ! f(3,i) are the B elements/3, and will become c/3 upon solution
   f(4,1) = x(2)-x(1)
   f(3,2) = (f(1,2) - f(1,1))/f(4,1)
   do i = 2, k - 1
      f(4,i) = x(i+1) - x(i)
      f(2,i) = 2.0d0*(f(4,i-1) + f(4,i))
      f(3,i+1) = (f(1,i+1) - f(1,i))/f(4,i)
      f(3,i) = f(3,i+1) - f(3,i)
   end do
   ! Boundaries.
   f(2,2) = f(4,1) + 2.0d0*f(4,2)
   f(3,2) = f(3,2)*f(4,2)/(f(4,1) + f(4,2))
   f(2,k-1) = 2.0d0*f(4,k-2) + f(4,k-1)
   f(3,k-1) = f(3,k-1)*f(4,k-2)/(f(4,k-1) + f(4,k-2))
   ! Forward elimination.
   t = f(4,2)/f(2,2)
   f(2,3) = f(2,3) - t*(f(4,2) - f(4,1))
   f(3,3) = f(3,3) - t*f(3,2)
   do i = 4, k - 2
      t = f(4,i-1)/f(2,i-1)
      f(2,i) = f(2,i)-t*f(4,i-1)
      f(3,i) = f(3,i)-t*f(3,i-1)
   end do
   t = (f(4,k-2) - f(4,k-1))/f(2,k-2)
   f(2,k-1) = f(2,k-1) - t*f(4,k-2)
   f(3,k-1) = f(3,k-1) - t*f(3,k-2)
   ! Back substitution.
   f(3,k-1) = f(3,k-1)/f(2,k-1)
   do ib = 1, k - 4
      i = k - 1 - ib
      f(3,i) = (f(3,i) - f(4,i)*f(3,i+1))/f(2,i)
   end do
   f(3,2) = (f(3,2) - (f(4,2) - f(4,1))*f(3,3))/f(2,2)
   ! Reset d array to step size.
   f(4,1) = x(2) - x(1)
   f(4,k-1) = x(k) - x(k-1)
   ! Set f(3,1) for not-a-knot.
   f(3,1) = (f(3,2)*(f(4,1) + f(4,2)) - f(3,3)*f(4,1))/f(4,2)
   f(3,k) = f(3,k-1) + (f(3,k-1) - f(3,k-2))*f(4,k-1)/f(4,k-2)
   ! Compute the polynomial coefficients.
   do i = 1, k - 1
      f(2,i) = (f(1,i+1) - f(1,i))/f(4,i) - f(4,i)*(f(3,i+1) &
           + 2.0d0*f(3,i))
      f(4,i) = (f(3,i+1) - f(3,i))/f(4,i)
      f(3,i) = 3.0d0*f(3,i)
      f(4,i) = f(4,i)
   end do
end subroutine spline

end module opacity_spline
