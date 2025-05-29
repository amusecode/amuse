      module poly32
      use real_kind
!     Module for an n=3/2 polytrope

!     Number of data points used for polytrope interpolation
      integer, parameter :: n_poly_mesh = 7347

!     Range of xi and theta values
      real(double) :: xi32(n_poly_mesh)        ! Scaled coordinate
      real(double) :: theta32(n_poly_mesh)     ! Dimensionless density
      real(double) :: dthetadxi32(n_poly_mesh) ! -xi^2 dtheta/dxi

      contains

      subroutine read_poly32
      use filenames
      implicit none
      integer :: f, n

      f = get_free_file_unit()
      open (file = trim(evpath)//'/input/poly32.dat', unit = f)
      do n=1, n_poly_mesh
         read (f, *) xi32(n), theta32(n), dthetadxi32(n)
      end do
      close(f)
      end subroutine

!     Interpolate on the table for arbitrary values of xi
!     We use linear interpolation to avoid spurious oscillations from
!     splines
      subroutine interpolate_poly32(xi, theta, dthetadxi)
      implicit none
      real(double), intent(in) :: xi
      real(double), intent(out) :: theta, dthetadxi
      real(double) :: dxi, dxii
      integer :: n, nh, nl

      if (xi<0 .or. xi > xi32(n_poly_mesh)) return

      ! Binary search
      nh = n_poly_mesh+1
      nl = 1
      do while (nh > nl+1)
         n = nl + (nh - nl) / 2
         if (xi32(n) > xi) then     ! n is new upper bound
            nh = n
         else                       ! n is new lower bound
            nl = n
         end if
      end do 

      dxi = (xi - xi32(nl)) / (xi32(nh) - xi32(nl))
      dxii = 1.0d0 - dxi

      theta = dxii * theta32(nl) + dxi * theta32(nh)
      dthetadxi = dxii * dthetadxi32(nl) + dxi * dthetadxi32(nh)

      !print *, xi,        xi32(nl),         xi32(nh)
      !print *, theta,     theta32(nl),      theta32(nh)
      !print *, dthetadxi, dthetadxi32(nl),  dthetadxi32(nh)
      end subroutine

      end module
