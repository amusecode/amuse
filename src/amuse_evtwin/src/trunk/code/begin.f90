#include "assert.h"

subroutine star12_loop ( jo, ksv, kt5, kp, it, dty )
   use real_kind
   use init_dat
   use mesh
   use mesh_enc
   use control
   use settings
   use extrapolate_dh
   use explicit_functions
   use nucleosynthesis
   use test_variables
   use current_model_properties
   use constants
   use indices

   implicit none
   integer, intent(in) :: ksv, kt5, kp, it
   integer, intent(out) :: jo
   real(double), intent(inout) :: dty
   integer :: kr,jf
   integer :: iter,ik,Jstar,Jstop,itry
   real(double) :: rz,prz
   real(double) :: cdd_1, cdd_2
   real(double) :: dx12, dx23
   integer :: jo_flash = 0

   call printb ( jo, 1, it)
   if ( jo == 1 ) return
   if ( ktw == 2 ) then
      ! Average change wanted by *1
      cdd_1 = cdd  !zq(30)
      call printb ( jo, 2, it )
      ! Average change wanted by *2
      cdd_2 = cdd  !zq(30)
      ! Set average change to minimum of the two components
      cdd = min(cdd_1, cdd_2)  !zq(30)
   end if
   if ( kp.eq.0 ) return
   call nextdt ( dty, jo, it )
   if ( jo.eq.3 ) return
   iter = kr1

   ! Begin evolutionary loop of KP time steps:
   jnn = 1
   jter = 0
   prz = 0
   do_kp: do
      jo = 0
      ! Solve for structure, mesh, and major composition variables
      joc = 1

      call smart_solver ( iter, id, kt5, jo )
      call relax_timestep ( iter, id, kt5, jo, dty )

      ! If no convergence, restart from 2 steps back, DT decreased substantially
5     continue
      if (jo /= 0) then
         call backup ( dty, jo )
         if (use_quadratic_predictions) then
            call restore_previous_model()
            call predict_dh(dty, dh)
         end if
         if ( jo.eq.2 ) exit do_kp    ! Abort if timestep below limit
         goto 4
      end if

      ! If required, solve for minor composition variables; mesh, structure fixed.
      if (ktw > 1) call update_accretion_abundances2()
      !call smooth_spikes()
      if (allocated(dhnuc)) dhnuc(1:ktw, :, 1:kh_nuc) = 0.0d0
      call update_concentration_gradients()
      call update_radative_accelerations2()
      do_jstar: do Jstar = 1, ktw
         joc = Jstar + 1
         call solver ( kr_nucsyn, ie, kt5, jo )

         if (jo == 0) then
            call update_concentration_gradients()
            call update_radative_accelerations2()
            call solver ( kr_nucsyn, ie, kt5, jo )
         end if

         ! Set small abundances to 0 if no convergence
         if (jo /= 0) then
            jo = 0
            dhnuc(Jstar, :, 1:kh_nuc) = 0.0d0
            where (hnuc < 1.0d-20) hnuc = 0.0d0
            do ik = 1, kh_nuc
               if (h(idx_for_star(VAR_H1, jstar), ik) == 0.0d0) hnuc(jstar, 41, ik) = 0.0d0
            end do
            call solver ( kr_nucsyn, ie, kt5, jo )
         end if

         if (jo /= 0) then
            ! Smooth "spikes", if they're likely to cause the problem
            !call smooth_spikes()
            !jo = 0
            !dhnuc(Jstar, :, 1:kh_nuc) = 0.0d0
            !call solver ( kr_nucsyn, ie, kt5, jo )
            do_itry: do itry = 1, 5
               if (point_max_err > 1 .and. point_max_err < kh) then
                  dx12 = hnuc(jstar, var_max_err, point_max_err-1) - hnuc(jstar, var_max_err, point_max_err  )
                  dx23 = hnuc(jstar, var_max_err, point_max_err  ) - hnuc(jstar, var_max_err, point_max_err+1)
                  if ( dx12*dx23 < 0.0_dbl .and. H(idx_for_star(VAR_LNT, Jstar), point_max_err) < 6*CLN ) then
                     call mix_gridpoint(Jstar, point_max_err)
                     jo = 0
                     dhnuc(Jstar, :, 1:kh_nuc) = 0.0d0
                     call solver ( kr_nucsyn, ie, kt5, jo )
                     if (jo == 0) exit do_itry
                  else
                     exit do_itry
                  end if
               else
                  exit do_itry
               end if
            end do do_itry
         end if

         if (jo /= 0) then
            jo = 17
            exit do_jstar
         end if
      end do do_jstar  ! do Jstar = 1, ktw

      ! If no convergence, restart from 2 steps back, DT decreased substantially
      if (jo /= 0) then
         call backup ( dty, jo )
         if (use_quadratic_predictions) then
            call restore_previous_model()
            call predict_dh(dty, dh)
         end if

         if ( jo.eq.2 ) exit do_kp    ! Abort if timestep below limit
         goto 4
      end if

      ! If mass < 0, probably NaN, so exit with code 99
      ! (works for Lahey F95, but probably compiler-dependent!)
      if(h(VAR_MASS,1).lt.0.d0) jo = 99

      ! If model didn't converge, give up
      if ( jo >= 1 ) exit do_kp

      prz = rz
      call printb ( jo, 1, it )
      rz = rlf(1)
      if ( ktw.eq.2 ) then
        !average change wanted by *1
         cdd_1 = cdd  !zq(30)
         call printb ( jo, 2, it )
         rz = rlf(2)
        !average change wanted by *2
         cdd_2 = cdd  !zq(30)
        !set average change to minimum of the two components
         cdd = min(cdd_1, cdd_2)  !zq(30)
      end if

      ! If the helium or carbon flash code was triggered, remember this but do not break the evolution run yet
      jo_flash = 0
      if (jo == 10 .or. jo == 9 .or. jo == 8 .or. jo == 6) then
         jo_flash = jo
         jo = 0
      end if

      if ( jo.eq.14 ) goto 5

      ! Check for manual termination
      read (11, *) jstop
      if ( jstop.ne.0 ) jo = 15
      close (11)
      write (11, *) 0
      close (11)

      if ( jo >= 2 ) exit do_kp
      call update ( dty )
      call update2
      if (use_quadratic_predictions) call store_current_model()

      ! 2nd condition saves every model, doesn't work with mod:
      if ( mod(jmod, ksv) == 1 .or. ksv == 1 ) call output ( kp, 15, jo, jf )

4     continue
      call nextdt ( dty, jo, it )
      if (use_quadratic_predictions) call predict_dh(dty, dh)
      if ( jo == 3 ) exit do_kp
      if ( jo == 5 ) exit do_kp       ! Reached age limit
      if ( jnn >= 4 ) iter = kr2
      jnn = jnn + 1

      ! Don't break the run just as RLOF starts!
      kr = int(abs(rz)/(abs(rz - prz) + 1.0d-6))

      if (kp > 0 .and. jnn >= kp) exit do_kp
      if (.not.( (jnn<=kp .or. kp<0) .or. jhold<3 .or. jter>2 .or. kr<=10 )) exit do_kp
   end do do_kp     ! Evolutionary loop of KP time steps:  do - ~150 lines up

   ! Did we stop because of the He flash or a similar condition?
   if (jo_flash /= 0) jo = jo_flash
end subroutine star12_loop

!> Read first model from JIP (InPut) and write to JOP (OutPut)
!! Read init.dat from IT
!<
subroutine star12 ( jo, kh2, jch, jop, jip, ksv, kt5, it )
   use real_kind
   use init_dat
   use mesh
   use mesh_enc
   use control
   use settings
   use extrapolate_dh
   use explicit_functions
   use nucleosynthesis
   use test_variables
   use current_model_properties
   use constants
   use indices
   use resolve_helium_flash, only: make_post_flash_model

   implicit none
   integer, intent(in) :: kh2, ksv, kt5, jop, jip, it
   integer, intent(inout) :: jch
   integer, intent(out) :: jo
   integer :: kp,jf
   real(double) :: dty, tm, oa, ecc, bms, p1, omega(kh2), vi

   if(debug_level.ge.2) write(6,'(/,A20,6I12)')'star12: ',jo, jop, jip, ksv, it

   ! Read an initial model:
   jnn = 0
   expl_var(:, :, :) = 0.0d0

   ! Load models and prepare for the first timestep
   call beginn ( jip, kh2, jch, dty, kp, it, jo, jf )

   if ( jo /= -2 ) then
      call star12_loop ( jo, ksv, kt5, kp, it, dty )
      if (jo == 8 .and. construct_zahb_model) then
         tm  = h(VAR_MASS,  1)
         oa  = h(VAR_HORB,  1)
         ecc = h(VAR_ECC,   1)
         bms = h(VAR_BMASS, 1)
         omega = h(VAR_OMEGA, 1:kh2)
         vi  = h(VAR_INERT, 1)
         p1 = 2.0*cpi/(h(VAR_OMEGA, 1) * csday)
         write (jb, *) 'He flash (construct new model)'
         call make_post_flash_model(jo)
         tm  = h(VAR_MASS,  1)
         call remesh ( kh2, jch, bms, tm, p1, ecc, oa, 1, jf )
         ! Conserve angular momentum
         ! FIXME: this does not really work properly for differentially rotating stars
         h(VAR_OMEGA, 1:kh2) = omega*vi / h(VAR_INERT, 1)
         hpr = h
         write (jb,*) 'Start ZAHB'
         call star12_loop ( jo, ksv, kt5, kp, it, dty )
      end if
   end if

   ! Output the last converged model, unless it's at the He flash
   call output ( kp, jop, jo, jf )
   ! Save last model to file.mod as well (suggested by Selma)
   if ( .not. (mod(jmod, ksv) == 1 .or. ksv == 1) )  call output ( kp, 15, jo, jf )
   if ( jo/=8 .and. jo/=13 ) then
      if (isb == 2) then
         write (jb, *) 'Switch to *', 3 - jb
      else
         write (jb, *) 'Evolution done'
         write (jb, *) 'jo = ', jo
      end if
   end if
   if ( jo == 8 ) write (jb, *) 'He flash'
   if ( jo == 13 ) write (jb,*) 'Start ZAHB'
   flush ( jb )
end subroutine star12



!> Called from star12
!! - jip:  unit for input structure model
!! - kh2:  desired number of grid points (remesh)
!! - jch:  remesh option
!! - it:   unit of 'init.dat', used as a switch to detect ZAHB construction
!! - jf:   (output) bitflags indicating what information is saved in the model files: inertia, changes DH, nucleosynthesis etc.
!<
subroutine beginn ( jip, kh2, jch, dty, kp, it, jo, jf )
   use real_kind
   use mesh
   use mesh_enc
   use constants
   use settings
   use model_initialiser
   use control
   use test_variables
   use current_model_properties
   use binary_history
   use accretion_abundances
   use radiative_acceleration
   use nucleosynthesis, only: nucleosynthesis_enabled, Hnuc, DHnuc, Hnucpr
   use filenames
   use file_exists_module
   use indices

   implicit none
   integer, intent(in) :: jip, kh2
   integer, intent(inout) :: jch
   integer, intent(out) :: jf
   integer :: kp,it,jo
   integer :: ii,kb,jin,ik,ir,i
   real(double) :: dty,bms,ecc,p1,tm,oa
   real(double) :: bm_bck, bper_bck, age_bck
   real(double) :: nuc(50, nm)
   real(double) :: h1(nvar, nm), dh1(nvar, nm)

   logical :: status, post_he_flash
   character(len=500) :: fname

   if(debug_level.ge.2) write(6,'(/,A20,I12,ES20.5,8I12)')'beginn: ',jip, dty, kp, it, jo



   h(:, 1:kh) = 0.0d0
   dh(:, 1:kh) = 0.0d0
   hspn = 0.0d0
   rlf = 0.0d0
   zet = 0.0d0
   xit = 0.0d0
   age = 0.0d0
   mc = 0.0d0
   om = 0.0d0

   bm_bck = bm
   bper_bck = bper
   bm = 0.0d0
   bper = 0.0d0

   sm = 0.0d0
   enc = 0.0d0

   tfr = 0.0d0
   t0 = 0.0d0
   m0 = 0.0d0
   mta = 0.0d0
   om0 = 0.0d0
   omta = 0.0d0
   a0 = 0.0d0
   ata = 0.0d0
   e0 = 0.0d0
   eta = 0.0d0
   cdd = 0.0d0

   bp = 0.0d0
   horb = 0.0d0
   ro = 0.0d0
   ra2 = 0.0d0
   secc = 0.0d0
   tn = 0.0d0
   wmh = 0.0d0
   wmhe = 0.0d0
   mh = 0.0d0
   mhe = 0.0d0
   mco = 0.0d0

   be = 0.0d0
   lh = 0.0d0
   lhe = 0.0d0
   lc = 0.0d0
   lnu = 0.0d0
   lth = 0.0d0
   mcb = 0.0d0
   msb = 0.0d0
   rcb = 0.0d0
   tct = 0.0d0

   sect(1:9999) = 1.0d20
   kb = 1

   !     Just constructed a ZAHB model; DON'T homogenise the model,
   !     whatever the input .dat file says.
   !     Without this flag, the post-he flash model will have its
   !     composition reset to ZAMS composition if we started from a
   !     homogeneous model - which is not what we want.
   !     We should also force the accretion mode to exponential.
   post_he_flash = .false.
   if (jo == 13) then
      jch = min(3, jch)
      jo = -1
      cmi_mode = 1
      post_he_flash = .true.
   end if

   ! Allocate data for nucleosynthesis calculations
   call allocate_nucleosynthesis_data(kh2)

   ! Abort if the requested number of meshpoints is larger than the
   ! size of the array
   if (kh2 > NM) then
      write (0, *) 'cannot rezone to ', kh2, 'meshpoints. maximum size is ', NM, 'meshpoints.'
      jo = -2
      return
   end if

   ! Read data for initial model (often last model of previous run)
   ! e.g. SM = stellar mass, solar units; DTY = next timestep, years
   do ii=1, ktw
      call load_star_model(jip,ii, h, dh, nuc, sm,dty,age,bper,bms,ecc,p1,enc,kh,kp,jmod,jb,jin,jf)

      ! Store the nucleosynthesis data structures:
      call set_initl_nucleosynt_abundances(ii)
      if (nucleosynthesis_enabled .and. iand(jf, 8) == 8) Hnuc(ii,:,:) = nuc(:, :)

      if ( jmod == 0 .and. kp>0 ) kp = kp + 1
      if ( jip.eq.13 .or. jip.eq.14 ) dty = ct3*dty
      if (post_he_flash .and. ktw.eq.1) then                ! Directly after He flash, restore orbital period and binary mass
         bper = bper_bck
         bms = bm_bck/cmsn
      end if

      if ( ktw == 2 ) jb = 1
      write (jb, '(1x, 1p, 2d14.6, d17.9, 5d14.6, 0p, 6i5)') sm, dty, age, bper, bms, ecc, p1, enc, kh, kp, jmod, jb, jin, jf

      if (mutate) p1 = 1.0e10
      tm = cmsn*sm

      ! Decrease timestep for mutation runs
      if ( mutate .and. jmod>0 ) dty = min(dty, 1.0d3)

      ! Make sure the timestep is small when starting ZAHB construction
      if (it == 24) then
         age = 0.0d0
         dty = min(dty, 1.0d3)
      end if

      if ( ii == 1 ) then
         ! Initialise timestep
         dt = csy*dty

         ! Initialise orbital angular momentum
         bm = cmsn*bms
         om = bm - tm
         oa = cg1*tm*om*(cg2*bper/bm)**c3rd*sqrt(1.0d0 - ecc*ecc)

         if ( jb /= 1 ) then
            ! For *2 of a binary, read in the mass-loss history of *1 from io12
            !The output file at unit 3 = io12 may not be open yet!
            !> \todo Want to do this in a more central place, e.g. filenames.f90?
            inquire(unit=3,opened=status)
            if(.not.status) then
               write (fname, '(A,i1)') 'fort.',3
               if (.not. file_exists(fname)) then
                  open(unit = 3, file=trim(basename)//'.io12')
               end if
            end if

            do ik = 1, 9999
               read (3, 102, end=8) ms(ik), ww, sdt(ik), ww, sect(ik), se(ik), wx, scm(ik), ww,  sang(ik), ww, ir
               if ( ik.gt.1 .and. sect(ik).gt.1.0d19 ) goto 8   ! FIXME - replace with exit? - may fail if more code is added below
            end do
8           continue
102         format (3(1p, d16.9, 5d10.3, 0p, f8.3, 7f8.5, 3f8.4, /),  &
                     (1p, d16.9, 5d10.3, 0p, f8.3, 7f8.3, 3f8.4, /),  &
                      1p, d16.9, 5d10.3, 0p, 11f8.3, i6)

         end if
      end if

      ! Setup initial timestep, store (and restore) the stellar age.
      age_bck = age
      call nextdt ( dty, jo, it )
      age = age_bck

      ! Sanity check
      if (jch < 3 .and. kh2 /= kh) then
         write(0, *) 'Warning: input mesh not equal to requested mesh, but no remesh specified.'
         write(0, *) 'Will use input mesh of ', kh, ' mesh points.'
         write(0, *) 'Change KH2 or JCH in init.dat to remove this warning.'
      end if

      ! Reset nucleosynthesis abundances to their default values if the star is to be homogeneous.
      if (nucleosynthesis_enabled .and. jch == 4) call set_initl_nucleosynt_abundances(ii)

      ! optionally rezone the model, e.g. for different no. of meshpoints.
      ! also initialise some variables that were not in input model.
      if (use_smooth_remesher) then
         call remesher ( kh2, jch, bm, tm, p1, ecc, oa, ii, jf )
      else
         call remesh ( kh2, jch, bm, tm, p1, ecc, oa, ii, jf )
      end if

      ! ZAHB construction, we need to take special care here since the ZAHB
      ! starting model will not have the correct envelope abundance. This is
      ! corrected by accreting material of the right composition, but we
      ! cannot do this immediately because we need to H burning shell to stay
      ! active and dumping pure He on top of it will make it less numerically
      ! stable.
      ! Composition accretion is switched on in printb when the stellar mass
      ! reaches the desired core mass.
      if (it == 24) ccac = 0.0

      ! Store model for star 1 in hpr
      if ( ii == 1 ) then
         h1 = h
         dh1 = dh
      end if
   end do  ! do ii=1,ktw

   ! Join the model for star 1 (in h1/dh1) and the model for star 2 (in h/dh) into a binary
   if ( ktw == 2 ) then
      ! Make star 2 the secondary
      h(index_secondary_start+1:index_secondary_start+nvstar, 1:kh) = h(1:nvstar, 1:kh)
      dh(index_secondary_start+1:index_secondary_start+nvstar, 1:kh) = dh(1:nvstar, 1:kh)

      ! Make star 1 the primary and restore orbital elements
      h(1:index_secondary_start, 1:kh) = h1(1:index_secondary_start, 1:kh)
      dh(1:index_secondary_start, 1:kh) = dh1(1:index_secondary_start, 1:kh)
   end if

   ! Store current model as "previous" model
   hpr = h

   if (nucleosynthesis_enabled) then
      DHnuc = 0.0d0
      Hnucpr = Hnuc
   end if

   ! store some numbers for possible restart with BACKUP
   jhold = 2
   prev(2:81) = (/hspn, rlf, zet, xit, age, bm, mc, om, bper, sm, enc,  &   ! zq(1:17)
       0.d0, 0.d0, tfr, t0, m0, mta, om0, omta, a0, ata, e0, eta, cdd,  &   ! zq(18:30)
       bp, horb, ro, ra2, 0.d0, secc, tn, wmh, wmhe, mh, mhe, mco,      &   ! zq(31:43)
       0.d0, be, lh, lhe, lc, lnu, lth, mcb, msb, rcb, tct/)


   pprev(1:81) = prev(1:81)
   jm1 = jmod
   ! Determine whether I and phi are computed or not, for OUTPUT
   jf = 0
   do i = 11, 50
      if ( id(i).eq.12 .or. id(i).eq.14 ) jf = jf + 1
   end do
   if ( jf.eq.1 ) jf = 0
   ! Store derivative information
   if (store_changes) jf = jf + 4
   ! Store nucleosynthesis data
   if (nucleosynthesis_enabled) jf = jf + 8
   return
end subroutine beginn

