! Mesh dependent artificial energy term
! This file contains:
! Functions for setting up linear/spline interpolation tables
! Functions for evaluating those tables, either as a function of
! meshpoint (not particularly useful) or as a function of mass
! coordinate (more useful)
      module mesh_enc
      use mesh
      implicit none

      logical :: mutate = .false.           ! TRUE if we want to mutate rather than evolve a star
      logical :: usemenc = .false.          ! TRUE if we want to use the mesh-dependent term
      logical :: adj_mea = .false.          ! TRUE if we want to adjust mea to luminosity
      logical :: adj_comp = .false.         ! TRUE if we want to adjust the composition
      logical :: start_model_only = .false. ! TRUE if we only want to construct a starting model

      ! Mutation modes, how to check for convergence and what variables
      ! to use for matching the target entropy profile
      integer, parameter :: MMODE_EV = -1    ! Use logT and logf to calculate S
      integer, parameter :: MMODE_SONLY = -2 ! Use S from input file directly
      integer :: mutation_mode = MMODE_EV    ! Default mutation mode

      ! Do we want to use spline or linear interpolation?
      logical :: use_spline_interpolation = .false.

      ! Input file for target profiles
      integer :: tprof_in = 62

      ! init.dat file for mutations
      integer, parameter :: mutate_dat = 63

      ! Target profile and interpolation coefficients b, c, d
      double precision, save :: TH(NVAR, NM)
      double precision, save :: THb(NVAR, NM)
      double precision, save :: THc(NVAR, NM)
      double precision, save :: THd(NVAR, NM)

      ! Entropy
      double precision, save :: TSa(NM)
      double precision, save :: TSb(NM)
      double precision, save :: TSc(NM)
      double precision, save :: TSd(NM)

      ! Number of meshpoints in interpolation profiles
      integer, save :: ip_mesh_size

      ! Convergence monitor
      double precision, save :: MUTANT_H(NVAR, NM)
      double precision, save :: best_diffsqr
      double precision, save :: last_diffsqr
      double precision, save :: curr_diffsqr
      double precision, save :: curr_maxdiffsqr
      integer, save :: best_mod

      ! Output file for artificial energy term, corresponding to fort.13
      ! menc_out+1 corresponds to fort.14, while menc_out+2 contains a
      ! sequence of models that corresponds to the output from fort.15
      integer, parameter :: menc_out = 35

      ! Variable numbers for menc, mea and met
      integer, parameter :: NMENC = 22
      integer, parameter :: NMEA = 23
      integer, parameter :: NMET = 24

      ! Equation numbers for menc, mea, met
      integer, parameter :: EMENC = 16
      integer, parameter :: EMEA = 18
      integer, parameter :: EMET = 22

      contains
      
c subroutine spline_init (n, x, y, b, c, d)
c  Taken from Make Me A Star, because the spline code already in the evolution
c  code is too integrated into the opacity table (sigh...)
C This routine is a very slightly modified version of one from the book
C Computer Methods for Mathematical Computations, by George Forsythe, Mike
C Malcolm, and Cleve Moler. The original copy was obtained from the Netlib
C mathematical software file server
C (see http://www2.ucsc.edu/cats/sc/software/netlib/) with the command
C "send mmas_spline from fmm".
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating mmas_spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the mmas_spline.
c
c
      subroutine spline_init (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.d0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
      end do
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
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
c
c  forward elimination
c
      do i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
      end do
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.d0*c(n))
      do i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.d0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.d0*c(i)
      end do
      c(n) = 3.d0*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.d0
      d(1) = 0.d0
      b(2) = b(1)
      c(2) = 0.d0
      d(2) = 0.d0
      return
      end subroutine

c subroutine linear_init (n, x, y, b, c, d)
c Set up linear interpolation tables in a way that is compatible with
c iptable_eval.
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating mmas_spline
c
c    s(x) = y(i) + b(i)*(x-x(i))
c    c(i) = d(i) = 0.0
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
      subroutine linear_init (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
      integer nm1, ib, i
      double precision t

      nm1 = n-1
      if ( n .lt. 2 ) return

      b(:) = 0
      c(:) = 0
      d(:) = 0

      do i=2, nm1
         ! Average forward and backward derivatives
         !b(i) = 0.5 * ( (y(i-1) - y(i)) / (x(i-1) - x(i)) +  (y(i) - y(i+1)) / (x(i) - x(i+1)))
         b(i) = (y(i) - y(i+1)) / (x(i) - x(i+1))
      enddo
      b(1) = (y(1) - y(2)) / (x(1) - x(2))
      b(n) = (y(n-1) - y(n)) / (x(n-1) - x(n))
      return
      end subroutine

      subroutine iptable_init (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
      if (use_spline_interpolation) then
         call spline_init (n, x, y, b, c, d)
      else
         call linear_init (n, x, y, b, c, d)
      end if
      end subroutine

c function iptable_eval(n, u, x, y, b, c, d)
c  Again taken from Make Me A Star.
C This routine is from the book Computer Methods for Mathematical
C Computations, by George Forsythe, Mike Malcolm, and Cleve Moler.
C This copy was obtained from the Netlib mathematical software file server
C (see http://www2.ucsc.edu/cats/sc/software/netlib/) with the command
C "send seval from fmm".
c
c  this subroutine evaluates the cubic mmas_spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the mmas_spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of mmas_spline coefficients computed by mmas_spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
c  This code will also accept the linear interpolation tables calculated by
c  linear_init() in this module.
c
      function iptable_eval(n, u, x, y, b, c, d)
      implicit none
      double precision  :: iptable_eval
      integer, intent(in) :: n
      double precision , intent(in) :: u, x(n), y(n), b(n), c(n), d(n)
      integer :: i;
      integer :: j, k
      double precision :: dx

      i = 1;
      if ( i >= n ) i = 1
      if ( u < x(i) ) go to 10
      if ( u >= x(i+1) ) go to 30
c
c  binary search - NB: Assumes that the array x is stored largest-smallest
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u < x(k) ) then
         i = k 
      else 
         j = k
      end if
         
      !if ( u >= x(k) ) j = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate mmas_spline
c
   30 dx = u - x(i)
      iptable_eval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end function
      
      ! Read interpolated variable number j at mass coordinate m
      function get_iph(m, j)
      double precision :: get_iph;
      double precision, intent(in) :: m
      integer, intent(in) :: j
      integer :: i = 50;
      
      get_iph = iptable_eval(ip_mesh_size, m, TH(4,:), TH(j,:), THb(j,:), THc(j,:), THd(j,:))
      end function

      ! Read interpolated entropy as a function of the mass coordinate m
      function get_s(m)
      double precision :: get_s;
      double precision, intent(in) :: m
      integer :: i = 50;
      
      get_s = iptable_eval(ip_mesh_size, m, TH(4,:), TSa(:), TSb(:), TSc(:), TSd(:))
      end function

      ! Compute quantities from the equation of state for mass
      ! coordinate m
      subroutine calc_equation_of_state(m)
         use settings
         implicit none
         double precision, intent(in) :: m;
         ! Import common blocks from STARS code
         double precision :: XA(9), NA(9)
         double precision :: NEO, NIO, NZZ, AVM, NE;
         double precision :: AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA
         double precision :: SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR(23)
         double precision :: DELTA, PHI, EXT, FKT, FKR
         double precision :: AT, AF
         COMMON /ABUND / XA, NA, NEO, NIO, NZZ, AVM, NE
         COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA,
     &    SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR, DELTA, PHI, EXT, FKT, FKR

         ! Call equation of state to calculate P and rho
         XA(1) = get_iph(m, 5);
         XA(2) = get_iph(m, 9);
         XA(3) = get_iph(m, 10);
         XA(4) = get_iph(m, 16);
         XA(5) = get_iph(m, 3);
         XA(6) = get_iph(m, 11);
         XA(7) = CMG*CZS;
         XA(8) = CSI*CZS;
         XA(9) = CFE*CZS;
         XA(7) = max(1.0 - XA(1) - XA(2) - XA(3) - XA(4) - XA(5) - XA(6) - XA(8) - XA(9), 0.0d0);

         AF = get_iph(m, 1);
         AT = get_iph(m, 2);

         call statef(AF, AT);
      end subroutine

      ! Calculate the target entropy at masscoordinate m
      function calc_s(m)
         use settings
         implicit none
         double precision, intent(in) :: m;
         double precision :: calc_s;
         double precision :: AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA
         double precision :: SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR(23);
         double precision :: DELTA, PHI, EXT, FKT, FKR
         double precision :: AT, AF
         COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA,
     &    SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR, DELTA, PHI, EXT, FKT, FKR

         call calc_equation_of_state(m)
         calc_s = S
      end function

      ! Read the target model for interpolation; this is stored in the
      !  same format as the output models from the evolution code.
      subroutine read_target_model()
      USE MESH
      double precision :: sm, dty, age, per, bms, ecc, p1, enc;
      integer :: kh, kp, jmod, jb, jin, jf;
      integer :: ik, ij;
      integer :: ioerror;
      double precision :: XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW(14)
      double precision :: AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA,
     : SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP,  
     : RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO, ROO, RGNE, 
     : RGMG, RCCG, RPNG, EX, ENX, WMU
      double precision :: DELTA, PHI, EXT, FKT, FKR
      double precision :: LOLEDD, DG, EG, GRADT, ETH, EGR, R, QQ, QMU, 
     : WL, WCV, HP, TW, PHIMU, GMR, SEP, M3,PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW
      COMMON /VBLES / LOLEDD, DG, EG, GRADT, ETH, EGR, R, QQ, QMU, 
     : WL, WCV, HP, TW, PHIMU, GMR, SEP, M3,PX, SX,QA
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, 
     : SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP,  
     : RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO, ROO, RGNE, 
     : RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHI, EXT, FKT, FKR
     
      if (use_spline_interpolation) then
         print *, 'Using spline interpolation of target profile.'
      else 
         print *, 'Using linear interpolation of target profile.'
      end if
      ! Special: read m, S(m) and X(m) from input file
      if (mutation_mode == MMODE_SONLY) then
         print *, 'Reading S-only model'
         ! Read profile information
         read  (tprof_in, *, iostat=ioerror) sm, ip_mesh_size
         kh = ip_mesh_size
         if (ioerror /= 0) then
            close (tprof_in);
            print *, "Error reading target profile header"
            return;
         end if

         ! Read profile
         print *, "Reading target profile"
         do ik = 1, kh
            read  (tprof_in, *, iostat=ioerror) TH(4, ik), TSa(ik),
     &         TH(5, ik), TH(9, ik), TH(10, ik), TH(16, ik), TH(3, ik), TH(11, ik)
         end do
         close (tprof_in);
         
         ! Set up interpolation tables
         print *, "Calculating interpolation coefficients"
         call iptable_init (kh, TH(4,:), TH(5,:), THb(5,:), THc(5,:), THd(5,:))
         call iptable_init (kh, TH(4,:), TH(9,:), THb(9,:), THc(9,:), THd(9,:))
         call iptable_init (kh, TH(4,:), TH(10,:), THb(10,:), THc(10,:), THd(10,:))
         call iptable_init (kh, TH(4,:), TH(16,:), THb(16,:), THc(16,:), THd(16,:))
         call iptable_init (kh, TH(4,:), TH(3,:), THb(3,:), THc(3,:), THd(3,:))
         call iptable_init (kh, TH(4,:), TH(11,:), THb(11,:), THc(11,:), THd(11,:))
         call iptable_init (kh, TH(4,:), TSa(:), TSb(:), TSc(:), TSd(:));
         print *, "Interpolation initialisation done."
      
         ! Call EOS routine to find reaction rates/entropy (we'll need those)
         call check_convergence_to_target_structure();
         print *, "Target profile initialisation completed."
         return
      end if


      ! Read input parameters for the target profile
      read  (tprof_in, *, iostat=ioerror)
     & sm, dty, age, per, bms, ecc, p1, enc, kh, kp, jmod, jb, jin, jf      
     
      if (ioerror /= 0) then
         close (tprof_in);
         print *, "Error reading target profile header"
         return;
      end if
      
      ip_mesh_size = kh;
      
      ! Read profile
      print *, "Reading target profile"
      do ik = 1, kh
         read  (tprof_in, *, iostat=ioerror) (TH(ij,ik), ij = 1, jin)
      end do
      close (tprof_in);
      
      if (ioerror /= 0) then
         print *, "Error reading target profile data"
         return;
      end if
      
      ! Initialise interpolation coefficients  
      print *, "Calculating interpolation coefficients"
      do ij = 1, jin
         call iptable_init (kh, TH(4,:), TH(ij,:), THb(ij,:), THc(ij,:), THd(ij,:))
      end do

      do ik = 1, kh
         XH = TH(5, ik);
         XHE = TH(9, ik);
         XC = TH(10, ik);
         XN = TH(16, ik);
         XO = TH(3, ik);
         XNE = TH(11, ik);
         call statef(TH(1, ik), TH(2, ik));
         !write (1024, '(2E11.3)') GRADA, DG
         TSa(ik) = S;
      end do
      call iptable_init (kh, TH(4,:), TSa(:), TSb(:), TSc(:), TSd(:));
      print *, "Interpolation initialisation done."
      
      print *, "Target profile initialisation completed."
      
      end subroutine
      
      
      
      ! Compute the relative squared deviation from teh target model for
      ! meshpoint ik.
      subroutine compute_entropy_difference(ik, curr_maxdiffsqr, dssqr)
      USE MESH
      implicit none
      integer, intent(in) :: ik;             ! Meshpoint to work with
      double precision, intent(inout) :: curr_maxdiffsqr;
      double precision, intent(inout) :: dssqr;
      double precision :: sm, dty, age, per, bms, ecc, p1, enc;
      integer :: GKH, KTW, KW(260)
      double precision :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      double precision :: LOM, m, dm
      double precision :: XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW(14)
      double precision :: AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA,
     : SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP,  
     : RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO, ROO, RGNE, 
     : RGMG, RCCG, RPNG, EX, ENX, WMU
      double precision :: DELTA, PHI, EXT, FKT, FKR
      double precision :: s_target, f, ds, logA, logA_target;
      integer :: kh, kp, jmod, jb, jin, jf;
      integer :: ij, ikh;
      COMMON H, DH, EPS, DEL, DH0, GKH, KTW, KW
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, 
     : SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP,  
     : RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO, ROO, RGNE, 
     : RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHI, EXT, FKT, FKR

      m = H(4,ik);

      ! If we are outside the mass bounds for the interpolation data, do nothing
      if (m < TH(4, ip_mesh_size) .or. m > TH(4, 1)) then
         return;
      end if

      !f = 1.0/gkh;
      ! Relative importance is weighted with the fractional mass of the shell
      if (ik>1 .and. ik<gkh) then
         dm = H(4,ik-1) - H(4,ik+1);
      else if (ik == 1) then
         dm = 0;
      else
         dm = H(4,gkh);
      end if
      f = dm/H(4, 1)

      ! Get desired entropy at current meshpoint
      select case (mutation_mode)
         case (MMODE_EV)
            s_target = calc_s(m);
         case (MMODE_SONLY)
            s_target = get_s(m);
      end select

      ! Get actual entropy at current meshpoint from the equation of state
      ! Pass composition to the equation of state through the common block
      XH = H(5, ik);
      XHE = H(9, ik);
      XC = H(10, ik);
      XN = H(16, ik);
      XO = H(3, ik);
      XNE = H(11, ik);
      call statef(H(1, ik), H(2, ik));
      call nucrat(H(2, ik));

      ! Get relative deviation for this mass shell
      ds = f*(1.0D0 - S/s_target)**2

      dssqr = dssqr + ds;
      curr_maxdiffsqr = MAX(curr_maxdiffsqr, ds)
      end subroutine


      ! Quantify the difference between the current structure model and
      ! the target structure model.
      ! The precise method for doing this depends on the method used to
      ! construct the target model because this determines what
      ! variables are available.
      subroutine check_convergence_to_target_structure()
      USE MESH
      implicit none
      integer :: GKH, KTW, KW(260)
      double precision :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      double precision :: dssqr, tau;
      integer :: ik;
      COMMON H, DH, EPS, DEL, DH0, GKH, KTW, KW

      dssqr = 0.0D0;
      curr_maxdiffsqr = 0.0D0;
      
      !do ik = 40, gkh
      do ik = 1, gkh
         call compute_entropy_difference(ik, curr_maxdiffsqr, dssqr)
      end do
      
      ! Convergence monitor
      last_diffsqr = curr_diffsqr;
      curr_diffsqr = dssqr;
      
      end subroutine

      end module
