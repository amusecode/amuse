! This module is designed to work along the MUSE library.
! It takes an array with pressure, density, mass coordinate and
! composition for each mesh point (to be extended with angular momentum)
! and constructs a new stellar model, in the vein of mkmergermod.
! This facility is forked off to its own module because (for the moment)
! it makes use of numerical recipes code that shouldn't taint the main
! library file.
! Order of variables stored at each meshpoint:
!  Mass coordinate [Msun], Radius [Rsun], log density [cgs],
!  log pressure [cgs], XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE
! Mespoints should be passed in from *surface* to *centre*
      function import_stellar_merger(nmesh, numvar, model, age_tag)
      use twin_library
      use fudge_control
      use mesh_enc
      use constants
      use settings
      use init_dat
      implicit none
      integer :: import_stellar_merger
      integer, intent(in) :: nmesh, numvar
      double precision, intent(in) :: model(numvar, nmesh)
      double precision, intent(in) :: age_tag
      double precision :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      integer :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
      double precision :: XA(9), NA(9)
      double precision :: NEO, NIO, NZZ, AVM, NE
      COMMON /ABUND / XA, NA, NEO, NIO, NZZ, AVM, NE
      double precision :: AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA
      double precision :: SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR(23)
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR
      double precision :: CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9), KZN(9)
      COMMON /ATDATA/ CH2, CHI, COM, CAN, CBN, KZN
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
      DOUBLE PRECISION :: HPR(NVAR,NM), HT(60025+4*NM)
      COMMON /STORE / HPR, HT

      type(init_dat_settings) :: initdat
      integer :: KH2, KR1, KR2, KSV, KT5, JCH
      integer :: n, i, id, iter, status
      double precision :: mass, entropy_max, target_diffsqr, vma, tkh
      character*500 :: outputfilename, basename

      basename = "star";

      call push_init_dat(initdat, KH2, KR1, KR2, KSV, KT5, JCH)
      KR2 = get_number_of_iterations();
      KX = 0
      KY = 0
      KZ = 0
      ! FIXME: there are more mixing and mass loss options that can be
      ! set in the init.dat file, these need to be stored/restored as
      ! well!
      CTH = 0.0
      CRD = 0.0
      CMR = 0.0
      CMJ = 0.0
      CMV = 0.0
      CMK = 0.0
      CMNL = 0.0
      mass = model(1, 1)

      ! Invert equation of state
      print *, 'Inverting EoS...'
      do n=1, nmesh
         ! Convert mass fractions back baryon number fractions
         NA(1:9) = model(5:13, n)
         do i=1, 9
            XA(i) = NA(i) * CBN(i)/CAN(i)
         end do
         vma = sum(XA(1:9))
         XA(1:9) = XA(1:9) / vma
         TH(4, n) = model(1, n) * CMSN    ! Mass coordinate
         TH(5, n) = XA(1)                 ! Hydrogen abundance
         TH(9, n) = XA(2)
         TH(10, n) = XA(3)
         TH(16, n) = XA(4)
         TH(3, n) = XA(5)
         TH(11, n) = XA(6)
         call prtoft (model(4, n), model(3, n), TH(1, n), TH(2, n))
      end do
      ip_mesh_size = nmesh

      ! Construct interpolation tables
      print *, 'Constructing interpolation tables'
      do n = 1, 16
         call iptable_init (nmesh, TH(4,:),TH(n,:),THb(n,:),THc(n,:),THd(n,:))
      end do

      print *, 'Loading ZAMS star of mass', mass
      call flush_star
      id = load_zams_star(mass, 0.0D0)
      call select_star(id)

      ! Stage 1: match composition
      adj_comp = .true.
      impose_composition_factor = 0.0d0
      do iter=1, 100
         print *, 'Composition adjustment factor =', impose_composition_factor
         status = twin_evolve()
         call flush_star()
         if (status /= 0) then
            print *, '*** Failed to converge on timestep', iter, 'with code', status
            stop
         end if

         if (impose_composition_factor>=1.0D0) exit
         if (impose_composition_factor<=1.0D-2) then
            impose_composition_factor = min(1.5d0*impose_composition_factor, 1.0d0);
         else
            impose_composition_factor = min(1.2d0*impose_composition_factor, 1.0d0);
         end if
         impose_composition_factor = max(impose_composition_factor, 1.0d-4);

         call flush(6)
      end do
      
      ! Store output
      call flush_star()      
      outputfilename = trim(basename)//'.comp_mod'
      print *, 'Writing output to ', trim(outputfilename)
      call dump_twin_model(id, outputfilename);

      ! Stage 2: adjust entropy profile, keep composition fixed
      usemenc = .true.
      impose_entropy_factor = 0.0d0
      entropy_max = 1.0d2
      curr_diffsqr = 1.0d3
      best_diffsqr = 1.0d3
      target_diffsqr = EPS
      target_diffsqr = 1.0e-4
      call set_number_of_iterations(20)
      do iter=1, 100
         AGE = 0.0
         print *, 'Entropy adjustment factor =', impose_entropy_factor
         ! Construct the next stellar model in the pseudo-evolution
         ! sequence. Make sure the timestep is always close to the
         ! thermal time scale for the best accuracy. Make sure the
         ! solution is stable at this timestep by iterating on it while
         ! only changing the timestep.
         status = twin_evolve()
         TKH = 1.0D22*CG*H(4, 1)*2/(exp(H(7, 1))*H(8, 1)*CSY)
         do while (status == 0 .and. DT < 10.0*CSY)
            !print *, 'Grow timestep'
            DT = 1.01*DT
            status = twin_evolve()
         end do
         do while (status == 0 .and. DT > TKH * CSY)
            !print *, 'Shrink timestep'
            DT = DT*0.8
            status = twin_evolve()
         end do
         !print *, DT/CSY, TKH
         DT = TKH * CSY
         call flush_star()
         if (status /= 0) then
            print *, '*** Failed to converge on timestep', iter, 'with code', status
            if (impose_entropy_factor >= 1.0d0) exit;
            stop
         end if
         
         ! Check convergence, adjust entropy matching factor
         call check_convergence_to_target_structure()
         if ( curr_diffsqr < best_diffsqr ) then
            best_diffsqr = curr_diffsqr
            best_mod = iter
            MUTANT_H(1:24, 1:KH) = H(1:24, 1:KH)
         end if
         !WRITE (6, '(1P, 3D16.9, I6)'), CURR_MAXDIFFSQR, CURR_DIFFSQR, BEST_DIFFSQR, BEST_MOD
         print *, 'Converged to',BEST_DIFFSQR,'at',BEST_MOD,'now', iter
         

         if ( ( impose_entropy_factor>=entropy_max .and. best_diffsqr>curr_diffsqr ) .or. best_diffsqr<target_diffsqr ) exit
         if (impose_entropy_factor < 1.0) then
            impose_entropy_factor = min(1.5d0*impose_entropy_factor, entropy_max);
         elseif (impose_entropy_factor < 1.0d2) then
            impose_entropy_factor = min(1.25d0*impose_entropy_factor, entropy_max);
         elseif (impose_entropy_factor < 1.0d3) then
            impose_entropy_factor = min(1.05d0*impose_entropy_factor, entropy_max);
         else
            impose_entropy_factor = min(1.01d0*impose_entropy_factor, entropy_max);
         end if
         impose_entropy_factor = max(impose_entropy_factor, 1.0d-8);

         call flush(6)
      end do

      call pop_init_dat(initdat, KH2, KR1, KR2, KSV, KT5, JCH)
      H(1:24, 1:KH) = MUTANT_H(1:24, 1:KH)
      HPR(1:24, 1:KH) = MUTANT_H(1:24, 1:KH)
      adj_comp = .false.
      usemenc = .false.
      impose_entropy_factor = 0.0d0
      impose_composition_factor = 0.0d0
      ! Set timestep
      !DT = 1.0D3 * CSY
      DT = TKH * CSY
      ! Set age: make sure this is not reset when we go back one model
      AGE = age_tag
      QR(10) = AGE
      PPR(10) = AGE
      call flush_star()
      call set_star_iter_parameters( id, 10, 20, 0 )

      ! Store output
      outputfilename = trim(basename)//'.pmutate'
      print *, 'Writing output to ', trim(outputfilename)
      call dump_twin_model(id, outputfilename);

      import_stellar_merger = id
      end function
