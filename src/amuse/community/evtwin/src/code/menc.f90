! Mesh dependent artificial energy term
module mesh_enc
   use real_kind
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

   ! Input file for target profiles
   integer :: tprof_in = 62

   ! init.dat file for mutations
   integer, parameter :: mutate_dat = 63

   ! Target profile and interpolation coefficients b, c, d
   real(double), save :: TH(NVAR, NM)
   real(double), save :: THb(NVAR, NM)
   real(double), save :: THc(NVAR, NM)
   real(double), save :: THd(NVAR, NM)

   ! Entropy
   real(double), save :: TSa(NM)
   real(double), save :: TSb(NM)
   real(double), save :: TSc(NM)
   real(double), save :: TSd(NM)

   ! Nucleosynthesis data for merger remnant
   real(double), allocatable, save, private :: THnuc(:, :)
   real(double), allocatable, save, private :: THnucb(:, :)
   real(double), allocatable, save, private :: THnucc(:, :)
   real(double), allocatable, save, private :: THnucd(:, :)
   logical, save, private :: have_nucleosynthesis_data

   ! Number of meshpoints in interpolation profiles
   integer, save :: ip_mesh_size

   ! Convergence monitor
   real(double), save :: MUTANT_H(NVAR, NM)
   real(double), save :: best_diffsqr
   real(double), save :: last_diffsqr
   real(double), save :: curr_diffsqr
   real(double), save :: curr_maxdiffsqr
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



   ! Read interpolated variable number j at mass coordinate m
   function get_iph(m, j)
      use real_kind
      use interpolate

      implicit none
      real(double) :: get_iph
      real(double), intent(in) :: m
      integer, intent(in) :: j

      get_iph = iptable_eval(ip_mesh_size, m, TH(4,:), TH(j,:), THb(j,:), THc(j,:), THd(j,:))
   end function get_iph



   ! Get interpolated abundance for isotope j at mass coordinate m
   function get_ipnucabund(m, j)
      use real_kind
      use interpolate
      implicit none
      real(double) :: get_ipnucabund
      real(double), intent(in) :: m
      integer, intent(in) :: j
      integer :: nms

      if (.not. have_nucleosynthesis_data) then
         get_ipnucabund = 0.0d0
         return
      end if
      nms = ip_mesh_size
      get_ipnucabund = iptable_eval(nms, m, TH(4,1:nms), THnuc(j,1:nms), THnucb(j,1:nms), THnucc(j,1:nms), THnucd(j,1:nms))
   end function get_ipnucabund




   ! Return TRUE or FALSE depending on whether the loaded interpolation model has nucleosynthesis data or not
   function interpolation_has_nucabund()
      implicit none
      logical :: interpolation_has_nucabund
      interpolation_has_nucabund = have_nucleosynthesis_data
   end function



   ! Read interpolated entropy as a function of the mass coordinate m
   function get_s(m)
      use real_kind
      use interpolate

      implicit none
      real(double) :: get_s;
      real(double), intent(in) :: m

      get_s = iptable_eval(ip_mesh_size, m, TH(4,:), TSa(:), TSb(:), TSc(:), TSd(:))
   end function get_s

   ! Compute quantities from the equation of state for mass
   ! coordinate m
   subroutine calc_equation_of_state(m, eos)
      use real_kind
      use settings
      use eostate_types

      implicit none
      real(double), intent(in) :: m;
      type(eostate), intent(out) :: eos
      type(abundance) :: abund
      real(double) :: XA(9)
      real(double) :: AT, AF

      ! Call equation of state to calculate P and rho
      XA(1) = get_iph(m, 5)
      XA(2) = get_iph(m, 9)
      XA(3) = get_iph(m, 10)
      XA(4) = get_iph(m, 16)
      XA(5) = get_iph(m, 3)
      XA(6) = get_iph(m, 11)
      XA(7) = CMG*CZS
      XA(8) = CSI*CZS
      XA(9) = CFE*CZS
      XA(7) = max(1.0 - XA(1) - XA(2) - XA(3) - XA(4) - XA(5) - XA(6) - XA(8) - XA(9), 0.0d0)

      AF = get_iph(m, 1)
      AT = get_iph(m, 2)

      call statef ( af, at, xa, abund, eos )
   end subroutine calc_equation_of_state



   ! Calculate the target entropy at masscoordinate m
   function calc_s(m)
      use real_kind
      use settings
      use eostate_types

      implicit none
      real(double), intent(in) :: m
      real(double) :: calc_s
      type(eostate) :: eos

      call calc_equation_of_state(m, eos)
      calc_s = eos%S
   end function calc_s



   ! Read the target model for interpolation; this is stored in the
   !  same format as the output models from the evolution code.
   subroutine read_target_model()
      use real_kind
      use mesh, only: NVAR
      use settings
      use interpolate
      use eostate_types
      use structure_variables

      implicit none
      real(double) :: sm, dty, age, per, bms, ecc, p1, enc
      integer :: kp, jmod, jb, jin, jf
      integer :: ik, ij;
      integer :: ioerror;
      real(double) :: dh1(NVAR)
      real(double) :: XA(9)
      type(eostate) :: eos
      type(abundance) :: abund


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
         if (ioerror /= 0) then
            close (tprof_in);
            print *, "Error reading target profile header"
            return;
         end if

         ! Read profile
         print *, "Reading target profile"
         do ik = 1, ip_mesh_size
            read  (tprof_in, *, iostat=ioerror) TH(4, ik), TSa(ik),&
                 TH(5, ik), TH(9, ik), TH(10, ik), TH(16, ik), TH(3, ik), TH(11, ik)
         end do
         close (tprof_in);

         ! Set up interpolation tables
         print *, "Calculating interpolation coefficients"
         call iptable_init (ip_mesh_size, TH(4,:), TH(5,:), THb(5,:), THc(5,:), THd(5,:))
         call iptable_init (ip_mesh_size, TH(4,:), TH(9,:), THb(9,:), THc(9,:), THd(9,:))
         call iptable_init (ip_mesh_size, TH(4,:), TH(10,:), THb(10,:), THc(10,:), THd(10,:))
         call iptable_init (ip_mesh_size, TH(4,:), TH(16,:), THb(16,:), THc(16,:), THd(16,:))
         call iptable_init (ip_mesh_size, TH(4,:), TH(3,:), THb(3,:), THc(3,:), THd(3,:))
         call iptable_init (ip_mesh_size, TH(4,:), TH(11,:), THb(11,:), THc(11,:), THd(11,:))
         call iptable_init (ip_mesh_size, TH(4,:), TSa(:), TSb(:), TSc(:), TSd(:));
         print *, "Interpolation initialisation done."

         ! Call EOS routine to find reaction rates/entropy (we'll need those)
         call check_conv_to_target_structure()
         print *, "Target profile initialisation completed."
         return
      end if


      ! Read input parameters for the target profile
      read  (tprof_in, *, iostat=ioerror)&
           sm, dty, age, per, bms, ecc, p1, enc, ip_mesh_size, kp, jmod, jb, jin, jf

      if (ioerror /= 0) then
         close (tprof_in);
         print *, "Error reading target profile header"
         return;
      end if

      if (ip_mesh_size > NM) then
         write (*, '(A,I4,A,I4,A)') '*** Error: model has more than ', NM,' mesh points (need', ip_mesh_size,').'
         write (*, *) 'Please increase the number of meshpoints in mesh.f and recompile'
         stop
      end if

      ! Read profile
      print *, "Reading target profile"
      do ik = 1, ip_mesh_size
         read  (tprof_in, *, iostat=ioerror) (TH(ij,ik), ij = 1, jin)
      end do

      if (ioerror /= 0) then
         print *, "Error reading target profile data"
         close (tprof_in);
         return;
      end if

      ! Skip DH if present
      if (iand(jf, 4) == 4) then
         print *, 'Skipping DH data'
         do ik = 1, ip_mesh_size
            read (tprof_in, *, iostat=ioerror) dh1(1:jin)
         end do
      end if

      if (ioerror /= 0) then
         print *, "Error reading target profile data"
         close (tprof_in);
         return;
      end if

      ! Read nucleosynthesis data, if we have it
      have_nucleosynthesis_data = (iand(jf, 8) == 8)
      if (have_nucleosynthesis_data) then
         print *, 'Reading nucleosynthesis data'
         allocate(THnuc(50, ip_mesh_size))
         allocate(THnucb(50, ip_mesh_size))
         allocate(THnucc(50, ip_mesh_size))
         allocate(THnucd(50, ip_mesh_size))
         do ik = 1, ip_mesh_size
            read (tprof_in, *, iostat=ioerror) THnuc(1:50, ik)
            if (ioerror /= 0) exit
         end do
      end if

      if (ioerror /= 0) then
         print *, "Error reading target profile data"
         close (tprof_in);
         return;
      end if

      ! Close the file now that we are done with it
      close (tprof_in);

      ! Initialise interpolation coefficients
      print *, "Calculating interpolation coefficients"
      do ij = 1, jin
         call iptable_init (ip_mesh_size, TH(4,:), TH(ij,:), THb(ij,:), THc(ij,:), THd(ij,:))
      end do

      ! Ditto abundances for nucleosynthesis code
      if (have_nucleosynthesis_data) then
         do ij = 1, 50
            ik = ip_mesh_size
            call iptable_init (ip_mesh_size, TH(4,1:ik), THnuc(ij,1:ik), THnucb(ij,1:ik), THnucc(ij,1:ik), THnucd(ij,1:ik))
         end do
      end if

      do ik = 1, ip_mesh_size
         XA(1) = TH(5, ik)
         XA(2) = TH(9, ik)
         XA(3) = TH(10, ik)
         XA(4) = TH(16, ik)
         XA(5) = TH(3, ik)
         XA(6) = TH(11, ik)
         XA(7) = CMG*CZS
         XA(8) = CSI*CZS
         XA(9) = CFE*CZS
         call statef ( TH(1, IK), TH(2, IK), xa, abund, eos )
         TSa(ik) = eos%S
      end do
      call iptable_init (ip_mesh_size, TH(4,:), TSa(:), TSb(:), TSc(:), TSd(:))
      print *, "Interpolation initialisation done."
      print *, "Target profile initialisation completed."
   end subroutine read_target_model



   ! Compute the relative squared deviation from teh target model for
   ! meshpoint ik.
   subroutine compute_entropy_difference(ik, curr_maxdiffsqr, dssqr)
      use real_kind
      use mesh
      use settings
      use eostate_types

      implicit none
      integer, intent(in) :: ik;             ! Meshpoint to work with
      real(double), intent(inout) :: curr_maxdiffsqr
      real(double), intent(inout) :: dssqr
      real(double) :: m, dm
      real(double) :: XA(9)
      real(double) :: s_target, f, ds
      type(eostate) :: eos
      type(abundance) :: abund

      m = H(4,ik)

      ! If we are outside the mass bounds for the interpolation data, do nothing
      if (m < TH(4, ip_mesh_size) .or. m > TH(4, 1)) then
         return
      end if

      !f = 1.0/gkh;
      ! Relative importance is weighted with the fractional mass of the shell
      if (ik>1 .and. ik<kh) then
         dm = H(4,ik-1) - H(4,ik+1)
      else if (ik == 1) then
         dm = 0;
      else
         dm = H(4,kh)
      end if
      f = dm/H(4, 1)

      ! Get desired entropy at current meshpoint
      s_target = 0
      select case (mutation_mode)
      case (MMODE_EV)
         s_target = calc_s(m);
      case (MMODE_SONLY)
         s_target = get_s(m);
      end select

      ! Get actual entropy at current meshpoint from the equation of state
      ! Pass composition to the equation of state through the common block
      XA(1) = H(5, ik)
      XA(2) = H(9, ik)
      XA(3) = H(10, ik)
      XA(4) = H(16, ik)
      XA(5) = H(3, ik)
      XA(6) = H(11, ik)
      XA(7) = CMG*CZS
      XA(8) = CSI*CZS
      XA(9) = CFE*CZS
      call statef ( H(1, IK), H(2, IK), xa, abund, eos )

      ! Get relative deviation for this mass shell
      ds = f*(1.0D0 - eos%S/s_target)**2

      dssqr = dssqr + ds
      curr_maxdiffsqr = MAX(curr_maxdiffsqr, ds)
   end subroutine compute_entropy_difference


   ! Quantify the difference between the current structure model and
   ! the target structure model.
   ! The precise method for doing this depends on the method used to
   ! construct the target model because this determines what
   ! variables are available.
   subroutine check_conv_to_target_structure()
      use real_kind
      use mesh

      implicit none
      real(double) :: dssqr
      integer :: ik;

      dssqr = 0.0D0;
      curr_maxdiffsqr = 0.0D0;

      !do ik = 40, gkh
      do ik = 1, kh
         call compute_entropy_difference(ik, curr_maxdiffsqr, dssqr)
      end do

      ! Convergence monitor
      last_diffsqr = curr_diffsqr;
      curr_diffsqr = dssqr;

   end subroutine check_conv_to_target_structure

end module mesh_enc
