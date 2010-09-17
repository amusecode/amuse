!  Helper function: calculate degeneracy parameter xi given lnf
function calculate_degeneracy_parameter(lnf)
   use real_kind
   
   implicit none
   real(double) :: calculate_degeneracy_parameter
   real(double), intent(in) :: lnf
   real(double) :: uf, wf, xi
   
   !   if (auf < -100.0) then
   !      xi = auf
   !   else if (auf > 100.0) then
   !      uf = dexp(auf)
   !      xi = 2.0d0*sqrt(uf)
   !   else
   uf = dexp(lnf)
   wf = dsqrt(1.0d0 + uf)
   xi = dlog(uf) + 2.0d0*(wf - dlog(1.0d0 + wf))
   !   end if
   
   calculate_degeneracy_parameter = xi
end function calculate_degeneracy_parameter



!  Find the parameter "log f" from the positron degeneracy parameter XI
function find_positron_f(xi)
   use real_kind
   
   implicit none
   real(double), intent(in) :: xi
   real(double) :: find_positron_f
   real(double), external :: calculate_degeneracy_parameter
   real(double), external :: solve_f_interval
   
   !  Asymptotic expansions
   if (xi > 15.0) then
      find_positron_f = log(0.25*xi**2 + 1)
      return
   else if (xi < -5) then
      find_positron_f = xi - 2.0 + log(4.0)
      return
   end if
   !  If logf is below -100.0 we really don't care enough to determine it
   !  more accurately; all Fermi functions go to 0 as exp(xi) in this region
   !  anyway.
   find_positron_f = solve_f_interval(calculate_degeneracy_parameter, xi,   -5.d0, 15.d0) 
end function find_positron_f



!  Calculate Fermi-Dirac integrals for positrons, passing in the appropriate
!  values of F and G for Positrons
subroutine positron_fdirac ( f, g, dp )
   use real_kind
   
   implicit none
   real(double), intent (in) :: f, g
   real(double), intent (out) :: dp(16)
   
   !  If the degeneracy parameter for positrons is very small, we don't care
   !  because all contributions are going to be essentially 0 anyway.
   if (f < 1.0d-300) then
      dp(:) = 0.0d0
      return
   end if
   
   !  Compute positron FD functions
   call fdirac(f, g, dp)
end subroutine positron_fdirac



!  Solve function: find the point x in the interval (x1, x2) where the
!  function value equals y, using bisection.
!  Returns x, or x1 - 1 if the function did not converge.
!  Code based on the Numerical Recipes routine RTSEC
function solve_f_interval(func,y,x1,x2)
   use real_kind
   
   implicit none
   real(double), parameter :: eps = 1.0d-8
   integer, parameter :: maxiter = 10
   real(double), intent(in) :: y, x1, x2
   real(double), external :: func
   real(double) :: solve_f_interval
   real(double) :: dx, fh, fl, xl, xh
   integer :: i
   
   fl = func(x1) - y
   fh = func(x2) - y
   if (abs(fl) < abs(fh)) then
      !     Swap highest/lowest variables
      xh = x1
      xl = x2
      !     Swap function values
      dx = fl
      fl = fh
      fh = dx
   else
      xl = x1
      xh = x2
   end if
   !  Loop until converged to the required accuracy
   do i = 1, maxiter
      !     Calculate correction
      dx = (xl - xh) * fh / (fh - fl)
      !     On each iteration, flip boundaries
      xl = xh
      fl = fh
      xh = xh + dx
      fh = func(xh) - y
      !     Converged?
      if (dabs(dx) < eps .or. fh == 0.0d0) exit
   end do
   solve_f_interval = xh
end function solve_f_interval



